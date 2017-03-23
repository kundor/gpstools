#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read VAPR ascii format into a numpy array.

VAPR ascii files have the format:
time,date,prn,snr,prn,snr,prn,snr,prn,snr

where each line has up to four prn,snr pairs (or as few as 0.)
prn is 2 digits, and snr is 2.1 (four characters), so each pair(with two commas)
is 8 characters.
time is HHMMSS.FFF, where the fractional part FFF is milliseconds.
date is DDMMYY (two digit year).
Thus each line is 10 + 1 + 6 + 8*4 = 49 characters.

Lines with no prn,snr pairs are discarded. Lines with the year between 70 and 99
inclusive are discarded (the receiver starts up producing spurious dates at the
GPS epoch, 060180).

Every once in a while there is a blank line, followed by a line like
HK:: BATT:4.13V TTR:13694s U2:UART2_0.TXT U3:/UART
(50 characters) which continues with a normal data line.
At this time, switch to reading the next file, cyclically.

Created on Mon Nov 21 10:14:06 2016

@author: nima9589
"""

import numpy as np
from numpy import sin, cos
#import glob
import re
from coords import xyz2llh
from gpstime import leapsecsnp, gpstotsecgps
from gpsazel import satcoeffs, mvec
from fetchnav import getsp3file
from utility import info

#%% Stuff to configure
udir = '/inge/scratch/UART/UART3'

outfile = 'uart_cat.txt'

LOC = (-1283649.0796, -4726431.0920, 4074789.6026)  # vpr3 site at Marshall

thresh = np.timedelta64(1, 'h')
"""A gap larger than this means to go back to the first file."""

gapthresh = np.timedelta64(5, 'm')
"""A gap larger than this is commented on."""

approxnumrec = 115920000
"""The number of snr records to preallocate"""

#%% Constants
rec = np.dtype([('prn', 'u1'), ('time', 'M8[us]'), ('el', 'f'), ('az', 'f'), ('snr', 'u2')])

# I don't know if SNRs below 10 are reported with leading zeros or not
goodline = re.compile(rb'([01][0-9]|2[0-3])' # hour
                      rb'([0-5][0-9]){2}\.[0-9]{3},' # minute second
                      rb'([0-2][0-9]|3[01])' # day
                      rb'(0[0-9]|1[012])' # month
                      rb'[0-9]{2}' # year
                      rb'(,(0[1-9]|[12][0-9]|3[012]),' # prn
                      rb'(0|[1-9][0-9]\.[0-9])){,4}[\r\n]+') # snr

hkline = re.compile(rb'HK:: BATT:[0-9]\.[0-9]{2}V TTR:[-0-9][0-9]{4}s U2:/?UART2[_/][0-9]+\.TXT U3:/?UART?3?')

#%% coordinate data for az/el computations
lat, lon, ht = xyz2llh(LOC)
trans = np.array([[-sin(lon), -sin(lat)*cos(lon), cos(lat)*cos(lon)],
                  [ cos(lon), -sin(lat)*sin(lon), cos(lat)*sin(lon)],
                  [ 0,         cos(lat),          sin(lat)]])

#%% Time stuff
def nptime(lin):
    """Convert first two fields into a numpy datetime64[us]

    Assumes two-digit year is in 21st century. Microsecond resolution.
    """
    return np.datetime64(b'20%b-%b-%bT%b:%b:%b' % (lin[15:17], lin[13:15], lin[11:13], lin[:2], lin[2:4], lin[4:10]), 'us')

#%%#-----------------------------------------------------------------------#%%#

def nexttime(file):
    """Returns the time in the next line of file."""
    pos = file.tell()
    line = next(file)
    file.seek(pos)
    return nptime(line)

def bigjump(file, curtime):
    """Return true if the time difference between the next line of file and curtime is big."""
    jump = nexttime(file) - curtime
    if jump > thresh:
        info('Jump of', jump.astype('m8[m]'), 'seen in file', file.name)
        return True
    return False

def endfile(file):
    """Return true if no data remains."""
    pos = file.tell()
    try:
        next(file)
    except StopIteration:
        return True
    file.seek(pos)
    return False

def readhk(file):
    """Read housekeeping (HK) line. Check for proper format, then discard."""
    hk = file.read(50)
    if not hk:
        info("Housekeeping line expected; instead found end of file?", file.name)
        return
    if not hkline.fullmatch(hk):
        info('The 50 characters following a blank line do not match expected HK format:\n', hk.decode())
        #raise OSError('Gave up parsing the files')

def windpastshortlines(file):
    pos = file.tell()
    for line in file:
        if len(line) > 20:
            break
        pos = file.tell()
    file.seek(pos)
    if b'HK' in line:
        readhk(file)

def checkline(line, name):
    if b'HK' in line:
        info('HK found in line not following a blank line in', name, '-- What?')
        #raise OSError('Gave up parsing these files')
        return False
    elif len(line) < 20:
        return False
    elif not goodline.fullmatch(line):
        info('Bad line in', name, line.decode().strip())
        return False
    elif int(line[15:17]) >= 70:
        return False
    return True

def closewhenempty(file):
    if endfile(file):
        info('Closing file', file.name)
        file.close()
        return True
    return False

def getnewtime(line, curtime, name):
    thetime = nptime(line)
    if curtime is not None:
        thegap = thetime - curtime
        if thegap > gapthresh:
            info('Gap of', thegap.astype('m8[m]'), 'in file', name)
    return thetime

def nextfilenum(i, files, curtime):
    """Increment i, cyclically.

    Reset to 0 when the next file doesn't have data for the current time.
    """
    i = (i + 1) % len(files)
    if files[i].closed:
        return 0
    if i and bigjump(files[i], curtime):
        return 0
    return i

def checkallclosed(files):
    if all(f.closed for f in files):
        info('All files closed')
    else:
        info('Not all files closed?!?')

def prepfiles(filenames):
    files = [open(f, 'rb') for f in filenames]
     # Use binary mode so that seek and iteration can mix.
    for f in files:
        windpastshortlines(f)
    return files

#%%
def concat(filenames):
    """A generator yielding all the good data lines, in chronological order, from the given files."""
    files = prepfiles(filenames)
    i = 0
    curtime = None
    while True:
        for line in files[i]:
            if line.isspace():
                break
            if checkline(line, files[i].name):
                curtime = getnewtime(line, curtime, files[i].name)
                yield line
        readhk(files[i])
        if closewhenempty(files[i]) and i == 0:
            break
        i = nextfilenum(i, files, curtime)
    checkallclosed(files)

def vfilter(filename):
    """A generator yielding all the good data lines from the given file."""
    with open(filename, 'rb') as file:
        for line in file:
            if line.isspace():
                readhk(file)
            elif checkline(line, filename):
                yield line

#%% Satellite positions
def sp3headtime(line):
    fields = ((3, 7), (8, 10), (11, 13), (14, 16), (17, 19), (20, 29))
    return np.datetime64('{}-{}-{}T{}:{}:{}'.format(*(line[b:e].strip().zfill(e - b) for b, e in fields)), 'us')

def readsp3(sp3file):
    with open(sp3file) as fid:
        line = fid.readline()
        time0 = sp3headtime(line)
        line = fid.readline().split()
        assert line[3] == '900.00000000'
        line = fid.readline()
        numsat =  int(line.split()[1])
        if numsat != 32:
            end = min(len(line), 10 + 3*numsat)
            prnlist = [int(line[i:i+2])-1 for i in range(10, end, 3)]
            if numsat > 17:
                line = fid.readline()
                end = min(len(line), 10 + 3*(numsat - 17))
                prnlist += [int(line[i:i+2])-1 for i in range(10, end, 3)]
    pl = np.loadtxt(sp3file, skiprows=22, usecols=(1, 2, 3), comments=['*', 'E', 'V'])
    pl.shape = (-1, numsat, 3) # reshape to 96xNx3 (-1 means 'fit the data')
    pl *= 1000
    if numsat != 32:
        B = np.zeros((len(pl), 32, 3))
        B[:, prnlist, :] = pl
        pl = B
    return pl, time0

def poslist(time):
    """Return an array of sp3 satellite positions for the day(s) around the given time."""
    min15 = np.timedelta64(15, 'm')
    pl, epoch = readsp3(getsp3file(time))
    endtime = epoch + min15*len(pl)
    if endtime - time < np.timedelta64(2, 'h'):
        gtime = min(endtime + np.timedelta64(4, 'h'),
                    np.datetime64('now', 'us') + np.timedelta64(8, 'h'))
        pl, epoch = readsp3(getsp3file(gtime))
    if time - epoch < np.timedelta64(2, 'h'):
        pl0, epoch0 = readsp3(getsp3file(time - np.timedelta64(12, 'h')))
        if epoch != epoch0 + len(pl0)*min15:
            info("Difference between sp3 files is not 900 seconds.")
# TODO: when dealing with ultrarapid files, we may have problems with this
        pl = np.concatenate((pl0, pl))
        epoch = epoch0
    return pl, epoch


#%%
def azel(sxloc):
    east, north, up = (sxloc - LOC) @ trans
    az = np.arctan2(east, north) * 180 / np.pi
    if az < 0:
        az += 360
    el = np.arctan2(up, np.sqrt(east*east + north*north)) * 180 / np.pi
    return az, el

def addrecords(SNR, line, curi, leaps, cofns):
    """Add snr records from uart-format line to array SNR, starting at index curi."""
    time = nptime(line) + leaps
    totsec = gpstotsecgps(time)
    words = line.strip().split(b',')
    for prn, snr in zip(words[2::2], words[3::2]):
        prn = int(prn)
        snr = round(float(snr)*10)
        sxloc = mvec(totsec) @ cofns[prn](totsec)
        az, el = azel(sxloc)
        try:
            SNR[curi] = (prn, time, el, az, snr)
        except IndexError:
            newlen = int(len(SNR)*1.2)
            info('Resizing SNR to', newlen)
            SNR.resize((newlen,))
            SNR[curi] = (prn, time, el, az, snr)
        curi += 1
    return curi

#%%#-----------------------------------------------------------------------#%%#

#files = glob.glob(os.path.join(udir, '*.TXT'))
"""
os.chdir(udir)
files = [str(n) + '.TXT' for n in range(1, 92)]
"""
#%% To just concatenate
def concatenate_files(files, outfile=outfile):
    with open(outfile, 'wb') as out:
        for line in concat(files):
            out.write(line)

#%% To make a numpy array from an iterable
def makearray(lineiter):
    SNR = np.empty(approxnumrec, dtype=rec) # an estimate of the number of records...
    line = next(lineiter)
    time0 = nptime(line)
    leaps = leapsecsnp(time0)
    pl, epoch = poslist(time0 + leaps)
    cofns = satcoeffs(pl)
    curi = addrecords(SNR, line, 0, leaps, cofns)
    for line in lineiter:
        curi = addrecords(SNR, line, curi, leaps, cofns)
    SNR.resize((curi,))
    return SNR

#%%
def process_files(files):
    """Make a numpy array of the records in a list of UART filenames."""
    return makearray(concat(files))

def process_catfile(filename=outfile):
    """Make a numpy array of the records in a concatenated, cleaned file."""
    with open(filename, 'rb') as fin:
        return makearray(fin)
