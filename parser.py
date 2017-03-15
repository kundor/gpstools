#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 09:04:30 2017

@author: nima9589
"""
from collections import defaultdict
import sys
import os
import numpy as np
from binex import read_record, info
from ascii_read import gpstotsecgps, poslist
from satpos import mvec, myinterp
from coords import llh2xyz, get_ellipsoid_ht
import plot

SNR_REC = np.dtype([('prn', 'u1'), ('time', 'M8[us]'), ('el', 'f'), ('az', 'f'), ('snr', 'u2')])
HK_REC = np.dtype([('rxid', 'u1'), ('time', 'M8[us]'), ('mac', 'u8'), ('lon', 'i4'),
                   ('lat', 'i4'), ('alt', 'i4'), ('volt', 'u2'), ('temp', 'i2'),
                   ('msgct', 'u2'), ('err', 'u1')])
GPSepoch = np.datetime64('1980-01-06', 'us')
PDIR='/usr/local/adm/config/apache/htdocs/i/vapr/VB001'

def allcoeffs(poslist, epoch, n=5):
    """Return fitted coefficients for poslist an N*3 numpy array of ECEF locations
    for the desired satellite every 15 minutes, and epoch the start time as a numpy datetime64."""
    sec0 = gpstotsecgps(epoch)
    return [np.linalg.solve(np.array([mvec(sec0 + i*900, n) for i in range(idx-n, idx+n+1)]),
                            poslist[idx-n:idx+n+1])
            for idx in range(n, len(poslist)-n)]

def satcoeffs(pl, epoch, n=5):
    secn = gpstotsecgps(epoch) + n*900
    return [myinterp(secn, 900, allcoeffs(pl[:,prn,:], epoch, n)) for prn in range(32)]

def enutrans(lat, lon):
    sin, cos = np.sin, np.cos
    return np.array([[-sin(lon), -sin(lat)*cos(lon), cos(lat)*cos(lon)],
                     [ cos(lon), -sin(lat)*sin(lon), cos(lat)*sin(lon)],
                     [ 0,         cos(lat),          sin(lat)]])

def azel(sxloc, rxloc, trans):
    """Given satellite location sxloc, receiver location rxloc (ECEF, meters)
    and a ENU translation matrix trans, return azimuth and elevation (degrees)."""
    east, north, up = (sxloc - rxloc) @ trans
    az = np.arctan2(east, north) * 180 / np.pi
    if az < 0:
        az += 360
    el = np.arctan2(up, np.sqrt(east*east + north*north)) * 180 / np.pi
    return az, el

def weeksow_to_np(week, sow):
    """Convert GPS week and (fractional) second-of-week to numpy datetime64 in GPS time."""
    return GPSepoch + np.timedelta64(week, 'W') + np.timedelta64(round(sow*1e6), 'us')

def sod(time):
    """Convert numpy datetime64 to seconds of day."""
    return (time - time.astype('M8[D]')) / np.timedelta64(1, 's')

def assignrec(arr, ind, val):
    """Assign val to arr[ind], growing it if necessary."""
    try:
        arr[ind] = val
    except IndexError:
        newlen = int(len(arr)*1.2)
        arr.resize((newlen,), refcheck=False)
        arr[ind] = val

def addrecords(SNR, time, snrs, curi, cofns, rxloc, trans):
    """Add snr records to array SNR, starting at index curi."""
    totsec = gpstotsecgps(time)
    for prn, snr in snrs:
        sxloc = mvec(totsec) @ cofns[prn-1](totsec)
        az, el = azel(sxloc, rxloc, trans)
        assignrec(SNR, curi, (prn, time, el, az, snr))
        curi += 1
    return curi

def reader(fid):
    """A generator that yields all the available records from fid whenever the stream ends.

    After resuming, the next values are the arrays enlarged by subsequently received
    records (the entire arrays are returned each time, not just the new records).
    """
    SNRs = defaultdict(lambda: np.recarray((1000,), dtype=SNR_REC))
    curind = defaultdict(int) # how far we've filled the corresponding SNR array
    rxloc = {} # for each rxid, its location, filled in when we get a HK record
    trans = {} # for each rxid, a ENU translation matrix, filled in when location known
    HK = np.recarray((1000,), dtype=HK_REC)
    curh = 0
    cofns = None # compute as soon as we find a good time record
    numtot, numempty, numearly, numnoloc = [defaultdict(int) for _ in range(4)]
    while True:
        try:
            rid, vals = read_record(fid)
        except EOFError:
            yield {rxid: SNR[:curind[rxid]] for rxid, SNR in SNRs.items()}, HK[:curh]
        if rid == 192:
            rxid, weeksow, snrs = vals
            if not snrs or weeksow[0] < 1000 or weeksow[0] > 2222:
                numempty[rxid] += (not snrs)
                numearly[rxid] += (weeksow[0] < 1000)
                numtot[rxid] += 1
                continue
            time = weeksow_to_np(*weeksow)
            if numempty[rxid] or numearly[rxid]:
                info("Skipped {} records ({} empty, {} early) at {:%H:%M:%S}.".format(
                        numtot[rxid], numempty[rxid], numearly[rxid], time.tolist()))
                numtot[rxid] = numempty[rxid] = numearly[rxid] = 0
            if rxid not in rxloc:
                numnoloc[rxid] += 1
                continue
            if numnoloc[rxid]:
                info("Skipped", numnoloc[rxid], "records before receiver", rxid, "location was known.")
                numnoloc[rxid] = 0
            if not cofns or time > end_time:
                pl, epoch = poslist(time)
                cofns = satcoeffs(pl, epoch)
                end_time = epoch + np.timedelta64(15, 'm') * (len(pl) - 7) # 1.5 hours before last entry
            curind[rxid] = addrecords(SNRs[rxid], time, snrs, curind[rxid], cofns,
                                      rxloc[rxid], trans[rxid])
        elif rid == 193:
            rxid = vals[0]
            vals = rxid, weeksow_to_np(*vals[1]), *vals[2:]
            assignrec(HK, curh, vals)
            curh += 1
            if rxid not in rxloc:
                if vals[1] > np.datetime64('now') + np.timedelta64(1, 'W'):
                    info('Rx', rxid, ': Not using location from future time', vals[1])
                    continue
                lon, lat, alt = vals[3:6]
                if lon == lat == 0:
                    info('Rx', rxid, ': Not using location 0,0')
                    continue
                lon /= 1e7
                if lon > 0 and vals[1] < np.datetime64('2017-03-14'):
                    lon *= -1 # All our data before this date is in the western hemisphere,
                    # but some was reported with the wrong sign; fix it.
                    info('Forcing longitude to western hemisphere, to correct '
                         'bad data before 2017-03-14')
                lat /= 1e7
                if alt:
                    alt /= 1000
                else:
                    info('Obtaining terrain altitude from Google and Unavco.')
                    alt = get_ellipsoid_ht(lat, lon)
                info("Receiver", rxid, "reported at", lon, "°E, ", lat, "°N, ", alt, " m.")
                lon *= np.pi / 180
                lat *= np.pi / 180
                rxloc[rxid] = llh2xyz(lat, lon, alt)
                trans[rxid] = enutrans(lat, lon)

def readall(fid):
    """Return all the records currently in fid.
    
    Return a dictionary of rxid to SNR recarrays, and an HK recarray.
    """
    return next(reader(fid))

def out_snr88(SNR, filename):
    """Output the records in the numpy array SNR to filename in snr88 format."""
    with open(filename, 'wt') as fid:
        for rec in SNR:
            fid.write('{:3} {:9.4f} {:9.4f} {:9.1f}      0.      0. {:6.2f}     0.\n'.format(
                       rec.prn, rec.el, rec.az, sod(rec.time), rec.snr/10))

def write_all_snr88(SNRs, prefix='rx'):
    """Write out each receiver's snr records in a file in the current directory.

    Filename is rx##DOY0.YY.snr88. This assumes that all the data is from the same
    GPS day.
    """
    for rxid, SNR in SNRs.items():
        epoch = SNR[0].time
        fname = prefix + '{:02}'.format(rxid) + epoch.tolist().strftime('%j0.%y.snr88')
        info('Writing', fname, '...')
        out_snr88(SNR, fname)

def macformat(mac):
    hs = '{:016X}'.format(mac)
    return ':'.join(hs[i:i+4] for i in range(0, len(hs), 4))

def hkreport(HK, file=None):
    """Print out the contents of housekeeping messages in the array HK."""
    for rec in HK:
        print('Rx{0.rxid:02} {1:%x %X}  {2}  {3:10.5f}ºE {4:9.5f}ºN {5:8.3f}m  {6:.2f}V {0.temp:2}ºC msgct:{0.msgct:5}  flags: {0.err:02X}'.format(rec, rec.time.tolist(), macformat(rec.mac), rec.lon/1e7, rec.lat/1e7, rec.alt/1000, rec.volt/100), file=file)

def translate(fid):
    """Translate BINEX data to ASCII formats.

    Read in data from the provided file object, and write out snr88 files
    and a housekeeping report in the current directory.
    """
    SNRs, HK = readall(fid)
    write_all_snr88(SNRs)
    hkfile = HK[0].time.tolist().strftime('HK_%j.%y.txt') # we need a site identifier!
    info('Writing', hkfile, '...')
    with open(hkfile, 'wt') as fid:
        hkreport(HK, fid)

# To read from stdin:
# translate(open(0, 'rb'))

def _symlink(src, dest):
    if os.path.islink(dest):
        os.remove(dest)
    elif os.path.exists(dest):
        info('Could not create symlink', dest, '->', src, 'because', dest, 'exists.')
        return
    os.symlink(src, dest)

def makeplots(SNRs, HK, symlink=True, pdir=PDIR, hours=4, endtime=None):
    old = os.getcwd()
    os.chdir(pdir)
    plot.allrises(SNRs)
    plot.tempvolt(HK, hours, endtime)
    hkfile = HK[0].time.tolist().strftime('HK_%j.%y.txt')
    with open(hkfile, 'wt') as fid:
        hkreport(HK, fid)
    if symlink:
        _symlink(hkfile, 'HK.txt')
    for rxid, SNR in SNRs.items():
        allsnr = plot.prn_snr(SNR, rxid, hours, endtime)
        nsx = plot.numsats(SNR, rxid, minelev=10, hrs=hours, endtime=endtime)
        avg = plot.meansnr(SNR, rxid, hours, endtime)
        if symlink:
            suf = '-RX{:02}.png'.format(rxid)
            _symlink(allsnr, 'ALLSNR' + suf)
            _symlink(nsx, 'NS' + suf)
            _symlink(avg, 'AVG' + suf)
            _symlink('TV'+avg[3:], 'TV' + suf)
    os.chdir(old)


if __name__ == "__main__": # When this file is run as a script
    if len(sys.argv) != 2:
        print("Usage:", sys.argv[0], "<file>\n"
              "<file> should be the path to VAPR BINEX data. "
              "snr88 files and a housekeeping report will be written out "
              "in the current directory.")
        sys.exit(1)
    with open(sys.argv[1], 'rb') as fid:
        translate(fid)

