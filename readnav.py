#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 15:30:54 2017

@author: nima9589
"""
import subprocess
from urllib.parse import urlparse
from urllib.request import urlopen
from contextlib import closing
import pandas as pd
import numpy as np
from utility import info

CDDIS_HRLY = 'ftp://cddis.gsfc.nasa.gov/gnss/data/hourly/%Y/%j/hour%j0.%yn.Z'
"""Where to fetch current hourly broadcast nav message; strftime format."""

BKG_15 = 'ftp://igs.bkg.bund.de/NTRIP/BRDC/%Y/%j/brdc%j0.%yn.Z'
"""Where to fetch current broadcast nav message, for a 24-hour sliding window
updated every 15 minutes from NTRIP streams, at the German BKG GNSS Data Center."""


'http://celestrak.com/GPS/SOF/latest.sof'
#: satellite outage listing; xml format; latest

'ftp://ftp.agi.com/pub/Catalog/SOF/Current_SOF.txt'
#: contains filename of the latest SOF file, which can be fetched under:
'ftp://ftp.agi.com/pub/Catalog/SOF/'
#...or:
'http://celestrak.com/GPS/SOF/'


def epoc2dt64(s):
    sfmt = "20{:>02}-{:>02}-{:>02}T{:>02}:{:>02}:{:>04}"
    return np.datetime64(sfmt.format(*(s.decode().split())))

def isurl(url):
    try:
        return urlparse(url).scheme != ""
    except Exception:
        return False

def getbytes(file):
    if isurl(file):
        with urlopen(file) as req:
            return req.read()
    try: # pathlib.Path
        return file.read_bytes()
    except Exception:
        pass
    try: # string
        file = open(file, 'rb')
    except TypeError:
        pass # already an open file?
    with closing(file):
        return file.read()

def getlines(file):
    byts = getbytes(file)
    if byts.startswith(b'\x1f\x9d'): # LZW-compressed
        byts = subprocess.check_output('gunzip', input=byts)
    return byts.splitlines()

def read_rinex_nav(file):
    """Parse RINEX nav message file to Pandas DataFrame.

    file: URL, file path, or open (binary-mode) file object.
    """
    lines = getlines(file)
    try:
        idx, = [i for i, l in enumerate(lines) if l.strip() == b'END OF HEADER']
    except ValueError:
        info('Incorrectly formatted RINEX 2 nav file (END OF HEADER does not appear once).')
        return
    lines = [b''.join(lines[i:i+8]).replace(b'D', b'E')
             for i in range(idx + 1, len(lines), 8)]
    names = ['PRN', 'ToC', 'Bias', 'Drift', 'Rate',
             'IODE', 'Crs', 'Δn', 'M0',
             'Cuc', 'Ecc', 'Cus', 'SqrtA',
             'ToE', 'Cic', 'Ω', 'CIS',
             'i0', 'Crc', 'ω', 'dΩ',
             'di', 'CL2', 'Week', 'L2P',
             'Acc', 'Health', 'TGD', 'IODC',
             'ToM', 'Fit']
    widths = np.array([2, 20] + [19]*29)
    widths[5::4] += 3
    dtype = [np.uint8, 'M8[ms]'] + [np.float]*29
    nnav = np.genfromtxt(lines, names=names, delimiter=widths, comments=None, converters={1: epoc2dt64}, dtype=dtype)
    return pd.DataFrame(nnav).set_index(['ToC', 'PRN'])


# ToC: Time of Clock (year, month, day, hour, minute, second)
#      20.3.3.3.3.1
# IODE: Issue of Data, Ephemeris -- changes with any change in the ephemeris parameters
#       whenever the terms of iode in subframes 2 and 3 and the 8 LSBs of the IODC in subframe 1
#       do not match, a data set cutover has occurred
# The IODE is an 8-bit number equal to the 8 LSBs of the 10-bit IODC of the same data set.
#       Will be different than any value transmitted by the same sat in the preceding 6 hours
# ToE: Time of Ephemeris (gps second of week)
#   ToE which is offset from an hour boundary is the first data set after a new upload cutover
#   20.3.3.4.3 ?
# IODC: Issue of Data, clock. The issue number of the data set, to detect any change
#       The transmitted IODC will be different from any value transmitted by the same SV
#       in the preceding 7 days
# ToM: Time of Message (gps second of week, derived from Z-count in Hand Over Word ?)
#      20.3.3.5.2.4 ?
