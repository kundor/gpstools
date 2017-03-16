#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 11:54:05 2017

@author: nima9589
"""
import os
import time
import sys
from collections import defaultdict
import numpy as np
from snrstats import calcsnrstat, gensnrnp
import plot
from parser import reader, readall, hkreport
from utility import info, debug

PLOTDIR='/usr/local/adm/config/apache/htdocs/i/vapr/VB001'

PLOT_IVAL = np.timedelta64(5, 'm')

def _symlink(src, dest):
    if src is None:
        return
    if os.path.islink(dest):
        os.remove(dest)
    elif os.path.exists(dest):
        info('Could not create symlink', dest, '->', src, 'because', dest, 'exists.')
        return
    os.symlink(src, dest)

def format_stats(rxid, stat, statp):
    """Format statistics as from calcsnrstat into an HTML table row."""
    return """
        <tr>
          <td>{:02}</td>
          <td>{:.1f}&ndash;{:.1f}</td>
          <td>{:.2f}</td>
          <td>({:.2f})</td>
          <td>{:.2f}</td>
          <td>({:.2f})</td>
        </tr>""".format(rxid, statp.min, stat.max, statp.mean, stat.mean, statp.std, stat.std)

def makeplots(SNRs, HK, symlink=True, pdir=PLOTDIR, hours=4, endtime=None):
    old = os.getcwd()
    os.chdir(pdir)
    #plot.allrises(SNRs) # skip for now
    tvs = plot.tempvolt(HK, hours, endtime)
    if symlink:
        for tv in tvs:
            _symlink(tv, tv[:7] + '.png')
    hkfile = HK[0].time.tolist().strftime('HK_%j.%y.txt')
    with open(hkfile, 'wt') as fid:
        hkreport(HK, fid)
    if symlink:
        _symlink(hkfile, 'HK.txt')
    snrtab = open('snrtab.html', 'wt')
    for rxid, SNR in SNRs.items():
        snrtab.write(format_stats(rxid, *calcsnrstat(gensnrnp(SNR))))
        allsnr = plot.prn_snr(SNR, rxid, hours, endtime)
        nsx = plot.numsats(SNR, rxid, minelev=10, hrs=hours, endtime=endtime)
        avg = plot.meansnr(SNR, rxid, hours, endtime)
        if symlink:
            suf = '-RX{:02}.png'.format(rxid)
            _symlink(allsnr, 'ALLSNR' + suf)
            _symlink(nsx, 'NS' + suf)
            _symlink(avg, 'AVG' + suf)
    snrtab.close()
    os.chdir(old)

def plotupdate(fname):
    olens = defaultdict(int)
    otic = np.datetime64('2000-01-01', 'ms')
    attempt = 0
    with open(fname, 'rb') as fid:
        for SNRs, HK in reader(fid):
            nlens = defaultdict(int, {rx: len(SNRs[rx]) for rx in SNRs}, hk=len(HK))
            tic = np.datetime64('now')
            if nlens == olens:
                attempt += 1
                if attempt > 5:
                    info('Now new records at', tic, '. Giving up.')
                    break
                info('No new records at', tic, '. Sleeping', attempt*5)
                time.sleep(attempt*5)
                continue
            attempt = 0
            olens = nlens
            debug('{:2} new records {} at'.format(sum(nlens.values()) - sum(olens.values()),
                                                 [nlens(rx) - olens[rx] for rx in nlens]),
                  tic, 'timestamped', max(max(s[-1].time for s in SNRs.values()), HK[-1].time))
            if tic - otic > PLOT_IVAL:
                info('Starting plotting at', tic)
                makeplots(SNRs, HK)
                info('Done at', np.datetime64('now'))
                otic = tic
            else:
                time.sleep(2)

def usage():
    print("Usage:", sys.argv[0], "<file> [hours] [endtime]\n"
          "  <file> should be the path to VAPR BINEX data. \n"
          "  hours (optional): how many hours to plot, default 4\n"
          "  endtime (optional): plots are made preceding this time (default now)\n"
          "                      Format 2017-03-15T17:06:59\n")
    sys.exit(1)

if __name__ == "__main__": # When this file is run as a script
    if len(sys.argv) < 2 or len(sys.argv) > 4:
        usage()
    hours = 4
    endtime = np.datetime64('now')
    if len(sys.argv) > 2:
        try:
            hours = float(sys.argv[2])
        except ValueError:
            print('Hours parameter must be a number')
            usage()
    if len(sys.argv) > 3:
        try:
            endtime = np.datetime64(sys.argv[3])
        except ValueError:
            print('Endtime must be given in the format YYYY-MM-DDTHH:MM:SS')
            usage()
    with open(sys.argv[1], 'rb') as fid:
        SNRs, HK = readall(fid)
    makeplots(SNRs, HK, hours=hours, endtime=endtime)
