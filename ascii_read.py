#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read VAPR ascii format.

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

import os
import re
import numpy as np
try:
    from utility import info
except ImportError:
    info = print

#%% Stuff to configure
outfile = 'uart_cat.txt'

thresh = np.timedelta64(1, 'h')
"""A gap larger than this means to go back to the first file."""

gapthresh = np.timedelta64(5, 'm')
"""A gap larger than this is commented on."""

#%% Constants

# I don't know if SNRs below 10 are reported with leading zeros or not
goodline = re.compile(rb'([01][0-9]|2[0-3])' # hour
                      rb'([0-5][0-9]){2}\.[0-9]{3},' # minute second
                      rb'([0-2][0-9]|3[01])' # day
                      rb'(0[0-9]|1[012])' # month
                      rb'[0-9]{2}' # year
                      rb'(,(0[1-9]|[12][0-9]|3[012]),' # prn
                      rb'(0|[1-9][0-9]\.[0-9])){,4}[\r\n]+') # snr

hkline = re.compile(rb'HK:: BATT:[0-9]\.[0-9]{2}V TTR:[-0-9][0-9]{4}s U2:/?UART2[_/][0-9]+\.TXT U3:/?UART?3?')

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

#%% To just concatenate
def concatenate_files(files, outfile=outfile):
    with open(outfile, 'wb') as out:
        for line in concat(files):
            out.write(line)

def concatenate_dir(udir, outfile=outfile):
    """Look for files in the given *udir* named UART3_0.TXT and UART3/#.TXT, where # is a number,
    and concatenate them into *outfile*.
    """
    old = os.getcwd()
    os.chdir(udir)
    try:
        if not os.path.isfile('UART3_0.TXT') or not os.path.isdir('UART3'):
            info('Could not find UART3_0.TXT and directory UART3')
            return
        files = ['UART3_0.TXT'] + [os.path.join('UART3', str(n) + '.TXT') for n in range(1, 999)]
        files = [f for f in files if os.path.isfile(f) and os.path.getsize(f)]
        #or, files = glob.glob(os.path.join(udir, 'UART3', '*.TXT')); then sort on
        # int(os.path.splitext(os.path.basename(file))[0])
        concatenate_files(files, outfile)
    finally:
        os.chdir(old)

if __name__ == "__main__": # When this file is run as a script
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description="Combine VAPR ASCII files",
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("directory", default='.', nargs='?',
                        help='should contain VAPR ASCII data, in files UART3_0.TXT '
                             'and UART3/#.TXT where # is a number.')
    parser.add_argument("outfile", default=outfile, nargs='?',
                        help='Where to write out ordered, cleaned data '
                             '(within <directory>, unless an absolute path is given).')
    args = parser.parse_args()
    concatenate_dir(args.directory, args.outfile)

