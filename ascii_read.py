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
At this time, switch to reading the next file.

Created on Mon Nov 21 10:14:06 2016

@author: nima9589
"""

import numpy as np
#import glob
import os
import re

udir = '/inge/scratch/UART/UART3'

thresh = np.timedelta64(1, 'h')
"""A gap larger than this means to go back to the first file."""

gap = np.timedelta64(5, 'm')
"""A gap larger than this is commented on."""


rec = np.dtype(('prn', 'u1'), ('time', 'M8[us]'), ('el', 'f'), ('az', 'f'), ('snr', 'u2'))

# I don't know if SNRs below 10 are reported with leading zeros or not
goodline = re.compile(rb'([01][0-9]|2[0-3])([0-5][0-9]){2}\.[0-9]00,([0-2][0-9]|3[01])(0[0-9]|1[012])[0-9]{2}'
                      '(,(0[1-9]|[12][0-9]|3[012]),(0|[1-9][0-9]\.[0-9])){,4}[\r\n]+')

hkline = re.compile(rb'HK:: BATT:[0-9]\.[0-9]{2}V TTR:[0-9]{5}s U2:/?UART2[_/][0-9]+\.TXT U3:/UART?')

def nptime(lin):
    """Convert first two fields into a numpy datetime64[us]

    Assumes two-digit year is in 21st century. Microsecond resolution.
    """
    return np.datetime64(b'20%b-%b-%bT%b:%b:%b' % (lin[15:17], lin[13:15], lin[11:13], lin[:2], lin[2:4], lin[4:10]))

def bigjump(file, curtime):
    pos = file.tell()
    try:
        line = next(file)
    except StopIteration:
        return True # don't use this file, go back to start
    file.seek(pos)
    if nptime(line) - curtime > thresh:
        return True

def concat(files, out):
    i = 0
    while True:
        for line in files[i]:
            if line.isspace():
                break
            elif b'HK' in line:
                print('HK found in line not following a blank line--What?')
                print('File ', files[i].name)
                raise OSError('Gave up parsing these files')
            if len(line) < 20:
                continue
            if int(line[15:17]) >= 70:
                continue
            if not goodline.fullmatch(line):
                print('Bad line ', line.decode())
                continue
            thetime = nptime(line)
            if thetime - curtime > gap:
                print('Gap of ', thetime - curtime, ' in file ', files[i].name)
            out.write(line)
            curtime = thetime
        hk = files[i].read(50)
        if not hk:
            if i:
                print("I'm confused, what is going on")
            break
        if not hkline.fullmatch(hk):
            print('The 50 characters following a blank line do not match expected HK format:\n', hk.decode())
            raise OSError('Gave up parsing the files')
        i = (i + 1) % len(files)
        if i and bigjump(files[i], curtime):
            i = 0

#files = glob.glob(os.path.join(udir, '*.TXT'))
os.chdir(udir)
files = [open(str(n) + '.TXT', 'rb') for n in range(1, 92)]
# Using binary mode so that seek and iteration can mix.
out = open('uart_cat.txt', 'wb')