#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 09:04:30 2017

@author: nima9589
"""
from collections import defaultdict
import numpy as np
from binex import *
from ascii_read import gpstotsecgps, poslist, mvec, satcoeffs

SNR_REC = np.dtype([('prn', 'u1'), ('time', 'M8[us]'), ('el', 'f'), ('az', 'f'), ('snr', 'u2')])
HK_REC = np.dtype([('rxid', 'u1'), ('time', 'M8[us]'), ('mac', 'u8'), ('lon', 'i4'),
                   ('lat', 'i4'), ('alt', 'i4'), ('volt', 'u2'), ('temp', 'i2'),
                   ('msgct', 'u2'), ('err', 'u1')])
GPSepoch = np.datetime64('1980-01-06', 'us')

def weeksow_to_np(week, sow):
    """Convert GPS week and (fractional) second-of-week to numpy datetime64 in GPS time."""
    return GPSepoch + np.timedelta64(week, 'W') + np.timedelta64(round(sow*1e6), 'us')

def assignrec(arr, ind, val):
    """Assign val to arr[ind], growing it if necessary."""
    try:
        arr[ind] = val
    except IndexError:
        newlen = int(len(arr)*1.2)
        arr.resize((newlen,))
        arr[ind] = val

def addrecords(SNR, weeksow, snrs, curi, cofns):
    """Add snr records to array SNR, starting at index curi."""
    time = weeksow_to_np(*weeksow)
    totsec = gpstotsecgps(time)
    for prn, snr in snrs:
        sxloc = mvec(totsec) @ cofns[prn](totsec)
        az, el = azel(sxloc)
        assignrec(SNR, curi, (prn, time, el, az, snr))
        curi += 1
    return curi

def readall(fid):
    SNRs = defaultdict(lambda: np.empty(1000, dtype=SNR_REC))
    curind = defaultdict(int)
    HK = np.empty(1000, dtype=HK_REC)
    curh = 0
    cofns = None
    numempty = 0
    while True:
        rid, vals = read_record(fid)
        if rid == 192:
            rxid, weeksow, snrs = vals
            if not snrs or weeksow[0] < 1000:
                numempty += 1
                continue
            if numempty:
                info("Skipped", numempty, "empty or early records.")
                numempty = 0
            if not cofns:
                time0 = weeksow_to_np(*weeksow)
                pl, epoch = poslist(time0)
                cofns = satcoeffs(pl)
            curind[rxid] = addrecords(SNRs[rxid], weeksow, snrs, curind[rxid], cofns)
        elif rid == 193:
            vals = vals[0], weeksow_to_np(*vals[1]), *vals[2:]
            assignrec(HK, curh, vals)
            curh += 1
    HK.resize((curh,))
    for rxid in curind:
        SNRs[rxid].resize((curind[rxid],))
    return SNRs, HK
