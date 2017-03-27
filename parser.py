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
from ascii_read import poslist
from gpstime import gpstotsecgps, gpsweekgps
from coords import llh2xyz, get_ellipsoid_ht
from utility import info, debug, mode

SNR_REC = np.dtype([('prn', 'u1'), ('time', 'M8[us]'), ('el', 'f'), ('az', 'f'), ('snr', 'u2')])
HK_REC = np.dtype([('rxid', 'u1'), ('time', 'M8[us]'), ('mac', 'u8'), ('lon', 'i4'),
                   ('lat', 'i4'), ('alt', 'i4'), ('volt', 'u2'), ('temp', 'i2'),
                   ('msgct', 'u2'), ('err', 'u1')])
GPSepoch = np.datetime64('1980-01-06', 'us')
SNR_ALLOC = 10**6 # how much to preallocate for these types
HK_ALLOC = 1000

class GrowArray:
    """Allows appending to a numpy array, reallocating as necessary."""
    def __init__(self, alloc, dtype, initarr=None):
        if initarr:
            self.curind = len(initarr)
            alloc = max(alloc, int(self.curind * 1.2))
            self.arr = np.empty((alloc,), dtype=dtype)
            self.arr[:self.curind] = initarr # copy so that we own the memory (we'll resize it)
        else:
            self.curind = 0
            self.arr = np.empty((alloc,), dtype=dtype)

    def append(self, val):
        try:
            self.arr[self.curind] = val
        except IndexError:
            newlen = int(len(self.arr) * 1.2)
            self.arr.resize((newlen,), refcheck=False)
            self.arr[self.curind] = val
        self.curind += 1

    def sliced(self):
        return self.arr[:self.curind]

def growSNR(arr=None):
    return GrowArray(SNR_ALLOC, SNR_REC, arr)

def growHK(arr=None):
    return GrowArray(HK_ALLOC, HK_REC, arr)

def enutrans(lat, lon):
    sin, cos = np.sin, np.cos
    return np.array([[-sin(lon), -sin(lat)*cos(lon), cos(lat)*cos(lon)],
                     [ cos(lon), -sin(lat)*sin(lon), cos(lat)*sin(lon)],
                     [ 0,         cos(lat),          sin(lat)]])

class SatPositions:
    """Class which precomputes satellite positions once a second, and gives them out
    when called. When the desired time is outside the range, a new set is precomputed."""
    @staticmethod
    def mvecs(times):
        P = 2 * np.pi / 86164.090530833 * times
        return np.hstack((np.ones((len(P), 1)),
                          np.sin(np.outer(P, np.arange(1, 6))),
                          np.cos(np.outer(P, np.arange(5, 0, -1)))))
        # Better yet, go with a constant 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1

    def __init__(self):
        self.endtime = np.datetime64('0')
        self.start = 0
        self.rxlocs = {}
        self.rxazel = {}

    def __call__(self, prn, time):
        if time >= self.endtime:
            self.update(time, 6)
        return self.pos[int(gpstotsecgps(time)) - self.start, prn - 1, :]

    def update(self, time, hours=None):
        """Compute an n*32*3 array of satellite positions for all PRNS, once per second,
        starting around *time* and lasting a maximum of the given *hours*."""
        pl, epoch = poslist(time)
        sec0 = gpstotsecgps(epoch)
        tim0 = gpstotsecgps(time)
        step0 = int((tim0 - sec0)/900)
        self.start = int(sec0 + step0*900)
        stop = int(sec0) + 900*(len(pl) - 6)
        endidx = len(pl) - 5
        if hours is not None:
            stop = min(stop, self.start + hours * 3600)
            endidx = min(endidx, 1 + step0 + hours * 4)
        COFS = np.array([[np.linalg.solve(self.mvecs(sec0 + np.arange(idx-5, idx+6)*900), pl[idx-5:idx+6, prn, :])
                      for prn in range(32)] for idx in range(step0, endidx)])
        M = self.mvecs(np.arange(self.start, stop))
        S = np.arange(0, 1, 1/900)
        S.shape = (900, 1, 1, 1)
        MS = 1 - S
        # K = np.kron(COFS[:-1], MS) + np.kron(COFS[1:], S)
        # This is equivalent to the kronecker product, but much faster:
        K = (COFS[:-1,np.newaxis,:,np.newaxis,:,np.newaxis,:,np.newaxis]
            * MS[np.newaxis, :, np.newaxis, :,np.newaxis, :,np.newaxis, :]
            + COFS[1:,np.newaxis,:,np.newaxis,:,np.newaxis,:,np.newaxis]
            * S[np.newaxis, :, np.newaxis, :,np.newaxis, :,np.newaxis, :])
        K.shape = (-1, *COFS.shape[1:])
        self.endtime = time + np.timedelta64(stop - int(tim0) - 1, 's')
        self.pos = np.sum(M.reshape((-1, 1, 11, 1)) * K, axis=2)
        for rxid in self.rxlocs:
            self.compazel(rxid)

    def rxadd(self, rxid, lat, lon, alt):
        lat *= np.pi / 180
        lon *= np.pi / 180
        self.rxlocs[rxid] = (llh2xyz(lat, lon, alt), enutrans(lat, lon))
        if self.start:
            self.compazel(rxid)

    def compazel(self, rxid):
        ENU = (self.pos - self.rxlocs[rxid][0]) @ self.rxlocs[rxid][1]
        az = np.arctan2(ENU[:, :, 0], ENU[:, :, 1]) * 180 / np.pi
        az[az < 0] += 360
        el = np.arctan2(ENU[:, :, 2], np.linalg.norm(ENU[:, :, :2], axis=2)) * 180 / np.pi
        self.rxazel[rxid] = np.stack((az, el), axis=-1)

def weeksow_to_np(week, sow):
    """Convert GPS week and (fractional) second-of-week to numpy datetime64 in GPS time."""
    return GPSepoch + np.timedelta64(week, 'W') + np.timedelta64(round(sow*1e6), 'us')

def sod(time):
    """Convert numpy datetime64 to seconds of day."""
    return (time - time.astype('M8[D]')) / np.timedelta64(1, 's')

def addrecords(SNR, time, snrs, satpos, rxid):
    """Add snr records to array SNR."""
    idx = int(gpstotsecgps(time)) - satpos.start
    for prn, snr in snrs:
        az, el = satpos.rxazel[rxid][idx, prn - 1, :]
        SNR.append((prn, time, el, az, snr))

def rx_locations(HK, satpos=None):
    """Return a dictionary of receiver locations determined from the housekeeping messages."""
    rxllh = {}
    HK = cleanhk(HK)
    HK = HK[HK['lon'] != 0]
    HK = HK[HK['lat'] != 0]
    for rx in np.unique(HK['rxid']):
        lon = np.median(HK[HK['rxid'] == rx]['lon']) / 1e7
        lat = np.median(HK[HK['rxid'] == rx]['lat']) / 1e7
        alt = np.median(HK[HK['rxid'] == rx]['alt']) / 1000
        if alt == 0:
            info('Obtaining terrain altitude for RX{:02} from Google and Unavco.'.format(rx))
            alt = get_ellipsoid_ht(lat, lon)
        info("RX{:02} found at {}°E, {}°N, {} m.".format(rx, lon, lat, alt))
        rxllh[rx] = (lon, lat, alt)
        if satpos:
            satpos.rxadd(lat, lon, alt)
    return rxllh

def reader(fid, preSNRs=None, preHK=None):
    """A generator that yields all the available records from fid whenever the stream ends.

    After resuming, the next values are the arrays enlarged by subsequently received
    records (the entire arrays are returned each time, not just the new records).
    """
    SNRs = defaultdict(growSNR)
    if preSNRs:
        for rx, SNR in preSNRs.items():
            SNRs[rx] = growSNR(SNR)
    HK = growHK(preHK)
    satpos = SatPositions()
    if preHK:
        rx_locations(preHK, satpos)
    numtot, numempty, numearly, numnoloc = [defaultdict(int) for _ in range(4)]
    thisweek = gpsweekgps(np.datetime64('now'))
    while True:
        try:
            rid, vals = read_record(fid)
        except EOFError:
            sliced = {rxid: SNR.sliced() for rxid, SNR in SNRs.items()}, HK.sliced()
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
            if rxid not in satpos.rxlocs:
                numnoloc[rxid] += 1
                continue
            if numnoloc[rxid]:
                info("Skipped", numnoloc[rxid], "records before receiver", rxid, "location was known.")
                numnoloc[rxid] = 0
            if time > satpos.endtime:
                debug('Updating satellite locations;', time, ' > ', satpos.endtime)
                satpos.update(time, 4)
            addrecords(SNRs[rxid], time, snrs, satpos, rxid)
        elif rid == 193:
            rxid = vals[0]
            vals = rxid, weeksow_to_np(*vals[1]), *vals[2:]
            HK.append(vals)
            if rxid not in satpos.rxlocs:
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
                satpos.rxadd(rxid, lat, lon, alt)
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
    for rx in np.unique(HK['rxid']):
        macs[rx] = mode(HK[HK['rxid'] == rx]['mac'])
    macmask = np.array([h['mac'] == macs[h['rxid']] for h in HK])
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

    Return a dictionary of rxid to SNR arrays, and an HK array.
    """
    return next(reader(fid))

def out_snr88(SNR, filename):
    """Output the records in the numpy array SNR to filename in snr88 format."""
    with open(filename, 'wt') as fid:
        for rec in SNR:
            fid.write('{:3} {:9.4f} {:9.4f} {:9.1f}      0.      0. {:6.2f}     0.\n'.format(
                       rec['prn'], rec['el'], rec['az'], sod(rec['time']), rec['snr']/10))

def write_all_snr88(SNRs, prefix='rx'):
    """Write out each receiver's snr records in a file in the current directory.

    Filename is rx##DOY0.YY.snr88. This assumes that all the data is from the same
    GPS day.
    """
    for rxid, SNR in SNRs.items():
        epoch = SNR[0]['time']
        fname = prefix + '{:02}'.format(rxid) + epoch.tolist().strftime('%j0.%y.snr88')
        info('Writing', fname, '...')
        out_snr88(SNR, fname)

def macformat(mac):
    hs = '{:016X}'.format(mac)
    return ':'.join(hs[i:i+4] for i in range(0, len(hs), 4))

def hkreport(HK, file=sys.stdout):
    """Print out the contents of housekeeping messages in the array HK."""
    for rec in HK:
        file.write('Rx{:02} '.format(rec['rxid']))
        try:
            file.write('{:%x %X}'.format(rec['time'].tolist()))
        except ValueError:
            file.write('{}'.format(rec['time']))
        file.write('  {}  {:10.5f}°E {:9.5f}°N {:8.3f}m  {:.2f}V {:2}°C msgct:{:5}  '
                   'flags: {:02X}\n'.format(macformat(rec['mac']),
                                            rec['lon']/1e7,
                                            rec['lat']/1e7,
                                            rec['alt']/1000,
                                            rec['volt']/100,
                                            rec['temp'],
                                            rec['msgct'],
                                            rec['err']))

def translate(fid):
    """Translate BINEX data to ASCII formats.

    Read in data from the provided file object, and write out snr88 files
    and a housekeeping report in the current directory.
    """
    SNRs, HK = readall(fid)
    write_all_snr88(SNRs)
    hkfile = HK[0]['time'].tolist().strftime('HK_%j.%y.txt') # we need a site identifier!
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

