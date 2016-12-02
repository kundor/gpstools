#!/usr/bin/env python3
from math import sqrt
import re

def gensnrs(lineiter):
    """Get the SNR values from each line in VAPR ASCII format.

    Lines should all be valid, so lineiter should either be a prefiltered
    file object, or run through filter(filename).
    """
    for line in lineiter:
        words = line.split(b',')
        for snr in words[3::2]:
            yield float(snr)

def gensnr89(filenames):
    """Get the SNR values from a list of files in snr89 format."""
    int89line = re.compile(r'(0[1-9]|[12][0-9]|3[012]) '
                           r'([0-8][0-9]|90) '
                           r'([012][0-9][0-9]|3[0-5][0-9]) '
                           r'[0-9]{1,5}\.[0-9]{1,3} '
                           r'0 0 ([0-9]{1,2}\.[0-9])\n')
    for f in filenames:
        with open(f) as fid:
            for line in fid:
                mch = int89line.fullmatch(line)
                if not mch:
                    print('Bad line', line.strip())
                    continue
                yield float(mch.group(4))

def reportspans(spanstarts, bins):
    pnum = sum(bins) / 100
    for beg in spanstarts:
        begind = round(beg*10)
        ct = sum(bins[begind:begind+128])
        print('{:.2f}% in [{}, {:.1f})'.format(ct/pnum, beg, beg + 12.8))

def bestspan(bins):
    maxct = 0
    maxi = 0
    for i in range(0, len(bins) - 128):
        ct = sum(bins[i:i+128])
        if ct > maxct:
            maxct = ct
            maxi = i
    pnum = sum(bins) / 100
    print('Best span [{}, {}) with {:.2f}%'.format(maxi/10, (maxi + 128)/10, maxct / pnum))

def paren(num):
    return '(' + str(num) + ')'

def calcsnrstat(snriter):
    snrmin = 100.
    snrminp = 100.
    snrmax = 0.
    snrmean = 0.
    snrmeanp = 0.
    num = 0
    nump = 0
    bins = [0.] * 600
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
        bins[round(snr*10)] += 1
    print('Min:', snrmin, paren(snrminp))
    print('Max:', snrmax)
    print('Mean:', snrmean, paren(snrmeanp))
    print('Num:', num, paren(nump))
    return snrmean, num, snrmeanp, nump, bins

def calcsnrstd(snriter, snrmean, snrmeanp, num, nump):
    sqrsum = 0.
    sqrsump = 0.
    for snr in snriter:
        sqrsum += (snr - snrmean)**2
        if snr > 0:
            sqrsump += (snr - snrmeanp)**2
    snrstd = sqrt(sqrsum / num)
    snrstdp = sqrt(sqrsump / nump)
    print('Std:', snrstd, paren(snrstdp))
    return snrstd, snrstdp

"""
with open('/inge/scratch/UART/UART3/uart_cat.txt', 'rb') as fid:
    snrmean, num, snrmeanp, nump, bins = calcsnrstat(gensnrs(fid))

reportspans([36, 38, 40], bins)
spanstarts = [36.4, 37, 37.6, 37.8, 38.2, 38.6, 39]
reportspans(spanstarts, bins)

with open('/inge/scratch/UART/UART3/uart_cat.txt', 'rb') as fid:
    calcsnrstd(gensnrs(fid), snrmean, snrmeanp, num, nump)

from ascii_read import filter
snrmean, num, snrmeanp, nump, bins = calcsnrstat(gensnrs(filter('uart3.txt')))

filenames = ['vpr3%d0.16.snr89' % d for d in range(188, 240)]
snrmean, num, snrmeanp, nump, bins = calcsnrstat(gensnr89(filenames))
"""
