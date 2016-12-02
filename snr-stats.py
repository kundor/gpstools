#!/usr/bin/env python3
from math import sqrt

spanstarts = [36, 38, 40]

def calcsnrstat(lineiter):
    snrmin = 100.
    snrminp = 100.
    snrmax = 0.
    snrmean = 0.
    num = 0
    spannum = [0] * len(spanstarts)
    for line in lineiter:
        words = line.split(b',')
        for snr in words[3::2]:
            snr = float(snr)
            snrmean = num/(num+1) * snrmean + snr/(num + 1)
            num += 1
            if snr > snrmax:
                snrmax = snr
            if snr < snrmin:
                snrmin = snr
            if snr > 0 and snr < snrminp:
                snrminp = snr
            for i, beg in enumerate(spanstarts):
                if beg <= snr < beg + 12.8:
                    spannum[i] += 1
    print('Min:', snrmin, '(' + str(snrminp) + ')')
    print('Max:', snrmax)
    print('Mean:', snrmean)
    print('Num:', num)
    pnum = num / 100
    for ct, beg in zip(spannum, spanstarts):
        print(ct/pnum, '% in [', beg, ', ', beg + 12.8, ')', sep='')
    return snrmean

def calcsnrstd(lineiter, snrmean):
    spannum = [0] * len(spanstarts)
    sqrsum = 0.
    for line in lineiter:
        words = line.split(b',')
        for snr in words[3::2]:
            snr = float(snr)
            sqrsum += (snr - snrmean)**2
            for i, beg in enumerate(spanstarts):
                if beg <= snr < beg + 12.8:
                    spannum[i] += 1
    snrstd = sqrt(sqrsum / num)
    print('Std:', snrstd)
    for ct, beg in zip(spannum, spanstarts):
        print(ct/pnum, '% in [', beg, ', ', beg + 12.8, ')', sep='')
    return snrstd

"""
with open('/inge/scratch/UART/UART3/uart_cat.txt', 'rb') as fid:
    snrmean = calcsnrstat(fid),

spanstarts = [36.4, 37, 37.6, 37.8, 38.2, 38.6, 39]

with open('/inge/scratch/UART/UART3/uart_cat.txt', 'rb') as fid:
    calcsnrstd(fid, snrmean)

from ascii_read import filter
snrmean = calcsnrstat(filter('uart3.txt'))
"""
