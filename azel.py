#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  1 09:10:01 2017

@author: nima9589
"""

import sys
import os
import numpy as np
from parser import SatPositions
from gpstime import GPSepoch, gpstotsecgps
from utility import debug

RXLOC = {
        'vpr2' : (39.94939333, -105.194425, 1727.794) # lat, lon, alt
        }

def snr89iter(fname):
    """An iterator returning prn, time, and SNR from each line in an snr89 file
    with integer az/el (to be ignored.)
    """
    with open(fname) as fid:
        for line in fid:
            words = line.split()
            if len(words) != 7:
                print('"' + line.rstrip() + '" is not a valid snr89 record (' + fname + ')')
                continue
            try:
                yield int(words[0]), float(words[3]), float(words[6])
            except ValueError:
                print('"' + line.rstrip() + '" is not a valid snr89 record (' + fname + ')')
                continue


def gpstotsecyeardoy(year, doy):
    ndt = np.datetime64(str(year), 'us')
    ndt += np.timedelta64(doy - 1, 'D')
    return gpstotsecgps(ndt)

def snrazel(snriter, year, doy, rxloc):
    """A generator that yields tuples of prn, el, az, sod, snr from an iterator that gives prn, sec of day, and SNR.

    Requires lat, lon in degrees, alt in meters.
    """
    satpos = SatPositions()
    satpos.rxadd(0, *rxloc)
    gts0 = int(gpstotsecyeardoy(year, doy))
    for prn, sod, snr in snriter:
        time = gts0 + int(sod)
        if time >= satpos.endgts:
            ntime = GPSepoch + np.timedelta64(time, 's')
            debug('Updating satellite locations;', ntime, ' > ', satpos.endtime)
            satpos.update(ntime, 4)
        idx = time - satpos.start
        try:
            az0, el0 = satpos.rxazel[0][idx, satpos.prndict[prn], :]
            az1, el1 = satpos.rxazel[0][idx+1, satpos.prndict[prn], :]
        except KeyError:
            continue
        frac = sod % 1
        az = (1 - frac)*az0 + frac*az1
        el = (1 - frac)*el0 + frac*el1
        yield prn, el, az, sod, snr

def addazel89(fname, oname):
    """Add azimuths and elevations to the snr89 file fname, writing to oname.

    Fname is parsed as SSSSJJJ#.YY.snr89, where SSSS is the site, JJJ is the day of year,
    and YY is taken as the two-digit year from 2000 to 2099.
    """
    site = fname[:4]
    doy = int(fname[4:7])
    year = 2000 + int(fname[9:11])
    try:
        rxloc = RXLOC[site]
    except KeyError:
        print('First four characters of filename, ' + site
              + ', are not a recognized site name. Receiver location unknown.')
        return
    with open(oname, 'wt') as fid:
        for snrec in snrazel(snr89iter(fname), year, doy, rxloc):
            fid.write('{:3} {:9.4f} {:9.4f} {:9.1f}      0.      0. {:6.2f}\n'.format(*snrec))


if __name__ == "__main__": # When this file is run as a script
    if len(sys.argv) != 3:
        print("Usage:", sys.argv[0], "<file> <newfile>\n"
              "<file> should be the filename of VAPR SNR89 data. "
              "Azimuths and elevations will be added, and <newfile> will be "
              "written out (in the same directory as <file>, unless an absolute "
              "path is given).")
        sys.exit(1)
    fname = sys.argv[1]
    ndir = os.path.dirname(fname)
    if ndir:
        os.chdir(ndir)
        fname = os.path.basename(fname)
    addazel89(fname, sys.argv[2])
