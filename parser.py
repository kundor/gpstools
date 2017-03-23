#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 09:04:30 2017

@author: nima9589
"""
from collections import defaultdict
import sys
import numpy as np
from binex import read_record
from ascii_read import gpstotsecgps, poslist, gpsweekgps
from satpos import mvec, myinterp
from coords import llh2xyz, get_ellipsoid_ht
from utility import info, mode

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

def rx_locations(HK):
    """Return a dictionary of receiver locations determined from the housekeeping messages."""
    HK = cleanhk(HK)
    HK = HK[HK.lon != 0]
    HK = HK[HK.lat != 0]
    rxloc = {}
    trans = {}
    for rx in np.unique(HK.rxid):
        lon = np.median(HK[HK.rxid == rx].lon) / 1e7
        lat = np.median(HK[HK.rxid == rx].lat) / 1e7
        alt = np.median(HK[HK.rxid == rx].alt) / 1000
        if alt == 0:
            info('Obtaining terrain altitude for RX{:02} from Google and Unavco.'.format(rx))
            alt = get_ellipsoid_ht(lat, lon)
        info("RX{:02} found at {}°E, {}°N, {} m.".format(rx, lon, lat, alt))
        rxloc[rx] = llh2xyz(lat * np.pi / 180, lon * np.pi / 180, alt)
        trans[rx] = enutrans(lat, lon)
    return rxloc, trans

def reader(fid, preSNRs=None, preHK=None):
    """A generator that yields all the available records from fid whenever the stream ends.

    After resuming, the next values are the arrays enlarged by subsequently received
    records (the entire arrays are returned each time, not just the new records).
    """
    snralloc = 10**6
    hkalloc = 1000
    SNRs = defaultdict(lambda: np.recarray((snralloc,), dtype=SNR_REC))
    curind = defaultdict(int) # how far we've filled the corresponding SNR array
    if preSNRs:
        for rx, SNR in preSNRs.items():
            rxlen = len(SNR)
            if rxlen > snralloc:
                newlen = int(rxlen*1.2)
                SNRs[rx] = np.recarray((newlen,), dtype=SNR_REC)
            SNRs[rx][:rxlen] = SNR
            curind[rx] = rxlen
    if preHK is None:
        HK = np.recarray((hkalloc,), dtype=HK_REC)
        curh = 0
        rxloc = {} # for each rxid, its location, filled in when we get a HK record
        trans = {} # for each rxid, a ENU translation matrix, filled in when location known
    else:
        curh = len(preHK)
        rxloc, trans = rx_locations(preHK)
        newlen = max(hkalloc, int(curh*1.2))
        HK = np.recarray((newlen,), dtype=HK_REC)
        HK[:curh] = preHK # copy so that we own the memory (we'll resize it)
    cofns = None # compute as soon as we find a good time record
    numtot, numempty, numearly, numnoloc = [defaultdict(int) for _ in range(4)]
    thisweek = gpsweekgps(np.datetime64('now'))
    while True:
        try:
            rid, vals = read_record(fid)
        except EOFError:
            sliced = {rxid: SNR[:curind[rxid]] for rxid, SNR in SNRs.items()}, HK[:curh]
            fid = (yield sliced) or fid
            continue
        except ValueError as e:
            info(e)
            continue
        if rid == 192:
            rxid, weeksow, snrs = vals
            if not snrs or weeksow[0] < 1000 or weeksow[0] > thisweek + 1:
                numempty[rxid] += (not snrs)
                numearly[rxid] += (weeksow[0] < 1000)
                numtot[rxid] += 1
                continue
            time = weeksow_to_np(*weeksow)
            if numempty[rxid] or numearly[rxid]:
                info("Skipped {:2} records ({:2} empty, {:2} early) at {:%H:%M:%S}.".format(
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
                if vals[1] > np.datetime64('now') + np.timedelta64(1, 'h'):
                    info('Rx', rxid, ': Not using location from future time', vals[1])
                    continue
                lon, lat, alt = vals[3:6]
                if lon == 0 or lat == 0:
                    info('Rx', rxid, ': Not using location with 0')
                    continue
                lon /= 1e7
                if lon > 0 and vals[1] < np.datetime64('2017-03-14'):
                    lon *= -1 # All our data before this date is in the western hemisphere,
                    # but some was reported with the wrong sign; fix it.
                    info('Forcing longitude to western hemisphere, to correct '
                         'bad data before 2017-03-14')
                lat /= 1e7
                if not (-90 <= lat <= 90 and -180 <= lon <= 360):
                    info('Not using bad location {}°E, {}°N'.format(lon, lat))
                    continue
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
        else:
            info('Unknown record {}:'.format(rid), vals.hex())

HK_RANGE = ((0, 63), # rxid
            [np.datetime64('2001-01-01'), np.datetime64('now') + np.timedelta64(1, 'W')], # time
            (0, 2**64 - 1), # mac
            (-180e7, 360e7), # longitude * 1e7
            (-90e7, 90e7), # latitude * 1e7
            (-1e6, 9e6), # altitude (mm)
            (0, 500), # voltage * 100
            (-90, 60), # temp (°C)
            (0, 2**16 - 1), # msgct
            (0, 255)) # err

def cleanhk(HK):
    """Return HK masked to only sensible entries."""
    macs = {}
    for rx in np.unique(HK.rxid):
        macs[rx] = mode(HK[HK.rxid == rx].mac)
    macmask = np.array([h.mac == macs[h.rxid] for h in HK])
    HK_RANGE[1][1] = np.datetime64('now') + np.timedelta64(1, 'h')
    return HK[np.logical_and.reduce([a <= HK[f] for (a, b), f in zip(HK_RANGE, HK_REC.names)]
                                  + [HK[f] <= b for (a, b), f in zip(HK_RANGE, HK_REC.names)]
                                  + [macmask])]

def validhk(hkr):
    """Is the given HK record good?"""
    HK_RANGE[1][1] = np.datetime64('now') + np.timedelta64(1, 'h')
    return all(a <= h <= b for (a, b), h in zip(HK_RANGE, hkr))

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

def hkreport(HK, file=sys.stdout):
    """Print out the contents of housekeeping messages in the array HK."""
    for rec in HK:
        file.write('Rx{:02} '.format(rec.rxid))
        try:
            file.write('{:%x %X}'.format(rec.time.tolist()))
        except ValueError:
            file.write('{}'.format(rec.time))
        file.write('  {}  {:10.5f}°E {:9.5f}°N {:8.3f}m  {:.2f}V {:2}°C msgct:{:5}  '
                   'flags: {:02X}\n'.format(macformat(rec.mac),
                                            rec.lon/1e7,
                                            rec.lat/1e7,
                                            rec.alt/1000,
                                            rec.volt/100,
                                            rec.temp,
                                            rec.msgct,
                                            rec.err))

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

if __name__ == "__main__": # When this file is run as a script
    if len(sys.argv) != 2:
        print("Usage:", sys.argv[0], "<file>\n"
              "<file> should be the path to VAPR BINEX data. "
              "snr88 files and a housekeeping report will be written out "
              "in the current directory.")
        sys.exit(1)
    with open(sys.argv[1], 'rb') as fid:
        translate(fid)

