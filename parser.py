#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 09:04:30 2017

@author: nima9589
"""
from collections import defaultdict
import sys
import numpy as np
from binex import read_record, info
from ascii_read import gpstotsecgps, poslist
from satpos import mvec, myinterp
from coords import llh2xyz

SNR_REC = np.dtype([('prn', 'u1'), ('time', 'M8[us]'), ('el', 'f'), ('az', 'f'), ('snr', 'u2')])
HK_REC = np.dtype([('rxid', 'u1'), ('time', 'M8[us]'), ('mac', 'u8'), ('lon', 'i4'),
                   ('lat', 'i4'), ('alt', 'i4'), ('volt', 'u2'), ('temp', 'i2'),
                   ('msgct', 'u2'), ('err', 'u1')])
GPSepoch = np.datetime64('1980-01-06', 'us')

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

def readall(fid):
    SNRs = defaultdict(lambda: np.recarray((1000,), dtype=SNR_REC))
    curind = defaultdict(int) # how far we've filled the corresponding SNR array
    rxloc = {} # for each rxid, its location, filled in when we get a HK record
    trans = {} # for each rxid, a ENU translation matrix, filled in when location known
    HK = np.recarray((1000,), dtype=HK_REC)
    curh = 0
    cofns = None # compute as soon as we find a good time record
    numempty = numearly = numnoloc = 0
    while True:
        try:
            rid, vals = read_record(fid)
        except EOFError:
            break
        if rid == 192:
            rxid, weeksow, snrs = vals
            if not snrs:
                numempty += 1
                continue
            if weeksow[0] < 1000:
                numearly += 1
                continue
            if numempty or numearly:
                info("Skipped", numempty, "empty and", numearly, "early records.")
                numempty = numearly = 0
            if rxid not in rxloc:
                numnoloc += 1
                continue
            if numnoloc:
                info("Skipped", numnoloc, "records before receiver", rxid, "location was known.")
                numnoloc = 0
            time = weeksow_to_np(*weeksow)
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
                lon, lat, alt = vals[3:6]
                lon /= -1e7  # KLUDGE: have non-negated west latitude?
                lat /= 1e7
                alt /= 1000
                info("Receiver", rxid, "reported at", lon, "°E, ", lat, "°N, ", alt, " m.")
                lon *= np.pi / 180
                lat *= np.pi / 180
                rxloc[rxid] = llh2xyz(lat, lon, alt)
                trans[rxid] = enutrans(lat, lon)
    HK.resize((curh,))
    for rxid in curind:
        SNRs[rxid].resize((curind[rxid],))
    return SNRs, HK


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
        print('Rx{0.rxid:02} {1:%x %X}  {2}  {3:.4f}ºE {4:.4f}ºN {5:.3f}m  {6:.2f}V {0.temp:2}ºC msgct:{0.msgct:5}  flags: {0.err:02X}'.format(rec, rec.time.tolist(), macformat(rec.mac), rec.lon/1e7, rec.lat/1e7, rec.alt/1000, rec.volt/100), file=file)

def translate(fid):
    """Translate BINEX data to ASCII formats.

    Read in data from the provided file object, and write out snr88 files
    and a housekeeping report in the current directory.
    """
    SNRs, HK = readall(fid)
    write_all_snr88(SNRs)
    epoch = HK[0].time
    hkfile = epoch.tolist().strftime('HK_%j.%y.txt') # we need a site identifier!
    info('Writing', hkfile, '...')
    with open(hkfile, 'wt') as fid:
        hkreport(HK, fid)

# To read from stdin:
# translate(open(0, 'rb'))

if __name__ == "__main__": # When this file is run as a script
    if len(sys.argv) != 2:
        print("Usage:", sys.argv[0], "<file>\n"
              "<file> should be the path to VAPR BINEX data. "
              "snr88 files and a housekeeping report will be written out "
              "in the current directory.")
        sys.exit(1)
    with open(sys.argv[1], 'rb') as fid:
        translate(fid)

