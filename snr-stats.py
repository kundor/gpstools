from math import sqrt
import re

def gensnrs(lineiter):
    """Get the SNR values from each line in VAPR ASCII format.

    Lines should all be valid, so lineiter should either be a prefiltered
    file object, or run through vfilter(filename).
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

def gensnrnp(arr):
    """Get the SNR values from a numpy recarray with SNR_REC datatype."""
    for rec in arr:
        yield rec.snr / 10

def reportspans(spanstarts, bins):
    pnum = sum(bins) / 100
    pnump = sum(bins[1:]) / 100
    for beg in spanstarts:
        begind = round(beg*10)
        ct = sum(bins[begind:begind+128])
        print('{:.2f}% ({:.2f}%) in [{}, {:.1f})'.format(ct/pnum, ct/pnump, beg, beg + 12.8))

def bestspan(bins):
    maxct = 0
    maxi = 0
    for i in range(0, len(bins) - 128):
        ct = sum(bins[i:i+128])
        if ct > maxct:
            maxct = ct
            maxi = i
    pnum = sum(bins) / 100
    pnump = sum(bins[1:]) / 100
    print('Best span [{}, {}) with {:.2f}% ({:.2f}%)'.format(maxi/10, (maxi + 128)/10, maxct / pnum, maxct / pnump))

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
    M2 = 0.
    M2p = 0.
    bins = [0.] * 600
    for snr in snriter:
        num += 1
        delta = snr - snrmean
        snrmean += delta / num
        delta2 = snr - snrmean
        M2 += delta * delta2
        if snr > snrmax:
            snrmax = snr
        if snr < snrmin:
            snrmin = snr
        if snr > 0:
            nump += 1
            delta = snr - snrmeanp
            snrmeanp += delta / nump
            delta2 = snr - snrmeanp
            M2p += delta * delta2
            if snr < snrminp:
                snrminp = snr
        bins[round(snr*10)] += 1
    snrstd = sqrt(M2 / (num - 1))
    snrstdp = sqrt(M2p / (nump - 1))
    print('Min:', snrmin, paren(snrminp))
    print('Max:', snrmax)
    print('Mean: {:.2f} ({:.2f})'.format(snrmean, snrmeanp))
    print('Num:', num, paren(nump))
    print('Std: {:.2f} ({:.2f})'.format(snrstd, snrstdp))
    return bins

"""
with open('/inge/scratch/UART/UART3/uart_cat.txt', 'rb') as fid:
    bins = calcsnrstat(gensnrs(fid))

reportspans([36, 38, 40], bins)
spanstarts = [36.4, 37, 37.6, 37.8, 38.2, 38.6, 39]
reportspans(spanstarts, bins)

from ascii_read import vfilter
bins = calcsnrstat(gensnrs(vfilter('uart3.txt')))

filenames = ['vpr3%d0.16.snr89' % d for d in range(188, 240)]
bins = calcsnrstat(gensnr89(filenames))
"""
