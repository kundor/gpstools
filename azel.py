#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  1 09:10:01 2017

@author: nima9589
"""

import numpy as np
from ascii_read import readsp3
from fetchnav import getsp3file
from parser import enutrans
from coords import llh2xyz, deg2rad
from gpstime import gpstotsecgps

SNR89rec = np.dtype([('prn', 'u1'), ('el', 'f'), ('az', 'f'), ('sod', 'f'), ('snr', 'f')])

RXLOC = {
        'vpr2' : (39.94939333, -105.194425, 1727.794) # lat, lon, alt
        }

def yeardoy2np(year, doy):
    ndt = np.datetime64(str(year), 'us')
    ndt += np.timedelta64(doy - 1, 'D')
    return ndt

def poslist(day):
    """Return an array of sp3 satellite positions for the three days around the given day."""
    min15 = np.timedelta64(15, 'm')
    pl, epoch = zip(readsp3(getsp3file(day - np.timedelta64(12, 'h'))),
                    readsp3(getsp3file(day + np.timedelta64(4, 'h'))),
                    readsp3(getsp3file(day - np.timedelta64(28, 'h'))))
    if epoch[1] != epoch[0] + len(pl[0])*min15 or epoch[2] != epoch[1] + len(pl[1])*min15:
        print("Difference between sp3 files is not 900 seconds.")
# TODO: when dealing with ultrarapid files, we may have problems with this
    pl = np.concatenate(pl)
    return pl, epoch[0]

def mvecs(times):
    P = 2 * np.pi / 86164.090530833 * times
    return np.hstack((np.ones((len(P), 1)),
                      np.sin(np.outer(P, np.arange(1, 6))),
                      np.cos(np.outer(P, np.arange(5, 0, -1)))))
    # Better yet, go with a constant 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1

def satpos(pl, epoch, gts):
    """Compute an n*32*3 array of satellite positions for all PRNS
    at given times."""
    sec0 = gpstotsecgps(epoch)
    step0 = int((gts[0] - sec0)/900)
    endidx = len(pl) - 5
    COFS = np.array([[np.linalg.solve(mvecs(sec0 + np.arange(idx-5, idx+6)*900),
                                      pl[idx-5:idx+6, prn, :])
                      for prn in range(32)]
                     for idx in range(step0, endidx)])
    M = mvecs(gts)
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
    return np.sum(M.reshape((-1, 1, 11, 1)) * K, axis=2)

def addazel(fname, oname):
    """Add azimuths and elevations to the snr89 file fname, writing to oname.

    Fname is parsed as SSSSJJJ0.YY.snr89, where SSSS is the site, JJJ is the day of year,
    and YY is taken as the two-digit year from 2000 to 2099.
    """
    site = fname[:4]
    doy = int(fname[4:7])
    year = 2000 + int(fname[9:11])
    try:
        rxllh = RXLOC[site]
    except KeyError:
        print('First four characters of filename, ' + site
              + ', are not a recognized site name. Receiver location unknown.')
        return
    rxloc = llh2xyz(deg2rad(rxllh))
    etrans = enutrans(*deg2rad(rxllh[:2]))
    day = yeardoy2np(year, doy)
    pl, epoch = poslist(day)
    gts0 = gpstotsecgps(day)
    snr = np.loadtxt(fname, dtype=SNR89rec, usecols=(0,1,2,3,6))
    satxyz = satpos(pl, epoch, snr['sod'] + gts0)
    ENU = (satxyz - rxloc) @ etrans
    az = np.arctan2(ENU[:, :, 0], ENU[:, :, 1]) * 180 / np.pi
    az[az < 0] += 360
    snr['az'] = az
    snr['el'] = np.arctan2(ENU[:, :, 2], np.linalg.norm(ENU[:, :, :2], axis=2)) * 180 / np.pi
    np.savetxt(oname, snr, fmt='%d %.4f %.4f %.1f 0 0 %.4f')
