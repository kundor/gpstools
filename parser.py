#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 09:04:30 2017

@author: nima9589
"""
from collections import defaultdict
from binex import *
import numpy as np

IN = open(0, 'rb')

SNR_REC = np.dtype([('prn', 'u1'), ('time', 'M8[us]'), ('el', 'f'), ('az', 'f'), ('snr', 'u2')])
HK_REC = np.dtype([('rxid', 'u1'), ('time', 'M8[us]'), ('mac', 'u8'), ('lon', 'i4'),
                   ('lat', 'i4'), ('alt', 'i4'), ('volt', 'u2'), ('temp', 'i2'),
                   ('msgct', 'u2'), ('err', 'u1')])
GPSepoch = np.datetime64('1980-01-06', 'us')

def weeksow_to_np(week, sow):
    """Convert GPS week and (fractional) second-of-week to numpy datetime64 in GPS time."""
    return GPSepoch + np.timedelta64(week, 'W') + np.timedelta64(sow*1e6, 'us')

def gpstotsecgps(ndt):
    """GPS seconds since epoch for a numpy datetime64 in GPS time."""
    return (ndt - GPSepoch) / np.timedelta64(1, 's')

def addrecords(SNR, weeksow, snrs, curi, cofns):
    """Add snr records to array SNR, starting at index curi."""
    time = weeksow_to_np(*weeksow)
    totsec = gpstotsecgps(time)
    for prn, snr in snrs:
        sxloc = mvec(totsec) @ cofns[prn](totsec)
        az, el = azel(sxloc)
        try:
            SNR[curi] = (prn, time, el, az, snr)
        except IndexError:
            newlen = int(len(SNR)*1.2)
            SNR.resize((newlen,))
            SNR[curi] = (prn, time, el, az, snr)
        curi += 1
    return curi

def readall():
    SNRs = defaultdict(lambda: np.empty(1000, dtype=SNR_REC))
    curind = defaultdict(int)
    HKs = np.empty(1000, dtype=HK_REC)
    while True:
        rid, vals = read_record(IN)
        if rid == 192:
            rxid, weeksow, snrs = vals
            curind[rxid] = addrecords(SNRs[rxid], weeksow, snrs, curind[rxid], cofns)

