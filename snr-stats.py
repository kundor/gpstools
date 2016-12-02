#!/usr/bin/env python3
from math import sqrt
import re

def gensnrs(lineiter):
    for line in lineiter:
        words = line.split(b',')
        for snr in words[3::2]:
            yield float(snr)

def gensnr89(filename):
    int89line = re.compile(r'(0[1-9]|[12][0-9]|3[012]) '
                           r'([0-8][0-9]|90) '
                           r'([012][0-9][0-9]|3[0-5][0-9]) '
                           r'[0-9]{1,5}\.[0-9] '
                           r'0 0 ([0-9]{2}\.[0-9])\n')
    with open(filename) as fid:
        for line in fid:
            mch = int89line.fullmatch(line)
            if not mch:
                print('Bad line', line.strip())
                continue
            yield float(mch.group(4))

def updatespans(snr, spanstarts, spannum):
    for i, beg in enumerate(spanstarts):
        if beg <= snr < beg + 12.8:
            spannum[i] += 1

def reportspans(spanstarts, spannum, pnum):
    for ct, beg in zip(spannum, spanstarts):
        print('{:.2f}% in [{}, {:.1f})'.format(ct/pnum, beg, beg + 12.8))

def paren(num):
    return '(' + str(num) + ')'

def calcsnrstat(snriter, spanstarts):
    snrmin = 100.
    snrminp = 100.
    snrmax = 0.
    snrmean = 0.
    snrmeanp = 0.
    num = 0
    nump = 0
    spannum = [0] * len(spanstarts)
    for snr in snriter:
        snrmean = num/(num+1) * snrmean + snr/(num + 1)
        num += 1
        if snr > snrmax:
            snrmax = snr
        if snr < snrmin:
            snrmin = snr
        if snr > 0:
            snrmeanp = nump/(nump+1) * snrmeanp + snr/(nump + 1)
            nump += 1
            if snr < snrminp:
                snrminp = snr
        updatespans(snr, spanstarts, spannum)
    print('Min:', snrmin, paren(snrminp))
    print('Max:', snrmax)
    print('Mean:', snrmean, paren(snrmeanp))
    print('Num:', num, paren(nump))
    reportspans(spanstarts, spannum, num/100)
    return snrmean, num

def calcsnrstd(snriter, spanstarts, snrmean, snrmeanp, num, nump):
    spannum = [0] * len(spanstarts)
    sqrsum = 0.
    sqrsump = 0.
    for snr in snriter:
        sqrsum += (snr - snrmean)**2
        if snr > 0:
            sqrsump += (snr - snrmeanp)**2
        updatespans(snr, spanstarts, spannum)
    snrstd = sqrt(sqrsum / num)
    snrstdp = sqrt(sqrsump / nump)
    print('Std:', snrstd, paren(snrstdp))
    reportspans(spanstarts, spannum, num/100)
    return snrstd

def calcsnrspans(snriter, spanstarts):
    num = 0
    spannum = [0] * len(spanstarts)
    for snr in snriter:
        num += 1
        updatespans(snr, spanstarts, spannum)
    reportspans(spanstarts, spannum, num/100)

"""
spanstarts = [36, 38, 40]

with open('/inge/scratch/UART/UART3/uart_cat.txt', 'rb') as fid:
    snrmean, num = calcsnrstat(gensnrs(fid), spanstarts)

spanstarts = [36.4, 37, 37.6, 37.8, 38.2, 38.6, 39]

with open('/inge/scratch/UART/UART3/uart_cat.txt', 'rb') as fid:
    calcsnrstd(gensnrs(fid), spanstarts, snrmean, num)

from ascii_read import filter
snrmean, num = calcsnrstat(gensnrs(filter('uart3.txt')), spanstarts)
"""
