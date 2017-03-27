#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 11:54:05 2017

@author: nima9589
"""
import os
import time
import sys
import shutil
from collections import defaultdict
import numpy as np
from snrstats import calcsnrstat
import plot
from parser import reader, readall, hkreport, cleanhk
from utility import info, debug, pushdir
import config

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

def rxline(rxid, HK):
    """A list entry for the web page for the given receiver."""
    lon = np.median(HK[HK['rxid'] == rxid]['lon']) / 1e7
    lat = np.median(HK[HK['rxid'] == rxid]['lat']) / 1e7
    alt = np.median(HK[HK['rxid'] == rxid]['alt']) / 1000
    return """    <li>
        RX{0:02}: <a href="#RX{0:02}-HK">HK</a> <a href="#RX{0:02}-SNR">SNR</a>
        {1:.5f}&deg;W, {2:.5f}&deg;N, {3:.3f}m.
        </li>""".format(rxid, lon, lat, alt)

def makeplots(SNRs, HK, symlink=True, pdir=None, snrhours=None, hkhours=None, endtime=None, **snrargs):
    if pdir is None:
        pdir = config.PLOTDIR
    if snrhours is None:
        snrhours = config.PLOT_SNR_HOURS
    if hkhours is None:
        hkhours = config.PLOT_HK_HOURS
    with pushdir(pdir):
        with open('updatetime.txt', 'at', encoding='utf-8') as fid:
            fid.write('. Updating...')
        #plot.allrises(SNRs) # skip for now
        tvs = plot.tempvolts(cleanhk(HK), hkhours, endtime)
        if symlink:
            for tv in tvs:
                _symlink(tv, tv[:7] + '.png')
        day = endtime or np.datetime64('now')
        hkfile = day.tolist().strftime('HK_%j.%y.txt')
        with open(hkfile, 'wt', encoding='utf-8') as fid:
            hkreport(HK, fid)
        if symlink:
            _symlink(hkfile, 'HK.txt')
        snrtab = open('snrtab-new.html', 'wt', encoding='utf-8')
        for rxid, SNR in SNRs.items():
            snrtab.write(format_stats(rxid, *calcsnrstat(SNR['snr'] / 10)))
            allsnr = plot.prn_snr(SNR, rxid, snrhours, endtime, **snrargs)
            nsx = plot.numsats(SNR, rxid, minelev=10, hrs=hkhours, endtime=endtime)
            avg = plot.meansnr(SNR, rxid, hkhours, endtime, minelev=10)
            if symlink:
                suf = '-RX{:02}.png'.format(rxid)
                _symlink(allsnr, 'ALLSNR' + suf)
                _symlink(nsx, 'NS' + suf)
                _symlink(avg, 'AVG' + suf)
        snrtab.close()
        os.replace('snrtab-new.html', 'snrtab.html')
        with open('rxlist.html', 'wt', encoding='utf-8') as fid:
            for rxid in SNRs:
                fid.write(rxline(rxid, HK))
        updatetime()

def midnightplots(SNRs, HK, day=None, ddir=None):
    """Make 24-hour archive page in *ddir* (default config.DAILYDIR).

    *day* is a numpy datetime64 in the desired day (default: the immediately
    preceding day.)
    *SNRs* and *HK* should include data for the desired 24 hours.
    """
    if ddir is None:
        ddir = config.DAILYDIR
    if day is None:
        day = np.datetime64('now') - np.timedelta64(6, 'h')
    day = day.astype('M8[D]')
    ddir = day.tolist().strftime(ddir)
    try:
        os.makedirs(ddir, exist_ok=True)
    except FileExistsError as e:
        info('Could not create directory', ddir, '\n    ', e)
        return
    etime = day + np.timedelta64(1, 'D')
    makeplots(SNRs, HK, pdir=ddir, snrhours=24, hkhours=24, endtime=etime,
              figlength=16, doazel=False, minelev=15)
    index = os.path.join(config.PLOTDIR, 'index.html')
    shutil.copy(index, ddir)

def updatetime(append=''):
    updstr = time.strftime('%b %d %Y %H:%M:%S UTC', time.gmtime()) + time.strftime(' (%H:%M:%S %Z)')
    with open('updatetime.txt', 'wt') as fid:
        fid.write(updstr + append)

def report_failure(rectic, dattic, pdir=None):
    """Record data failure in updatetime.txt"""
    if pdir is None:
        pdir = config.PLOTDIR
    failstr = ' (No data received since {}, timestamped {})'.format(rectic, dattic)
    with pushdir(pdir):
        updatetime(failstr)

def current_binex(time=None):
    if time is None:
        time = np.datetime64('now')
    return time.tolist().strftime(config.BINEX_FILES)

def after_midnight():
    """True if current time is within 15 minutes after UTC midnight."""
    now = np.datetime64('now')
    return np.timedelta64(0) < (now - now.astype('M8[D]')) < np.timedelta64(15, 'm')

def plotupdate(fname=None, handover=None, oldstate=None):
    """Follow stream and update web plot periodically.

    With no arguments, follow the current file and attempt handover at UTC midnight.
    With handover=False, stop when the file is exhausted.
    If fname is given, follow that file, updating plots until it's exhausted;
    if handover=True, then attempt to switch to the current file (this only makes
    sense if fname is yesterday's data.)
    """
    yesterday = False
    if oldstate:
        recgen, fid = oldstate
    elif fname is None and handover is not False:
        yesterday = np.datetime64('now') - np.timedelta64(1, 'D')
        fid = open(current_binex(yesterday), 'rb')
    elif fname is None:
        fid = open(current_binex(), 'rb')
    else:
        fid = open(fname, 'rb')
    ofile = os.path.abspath(fid.name)
    olens = defaultdict(int)
    oldtic = rectic = dattic = np.datetime64('2000-01-01', 'ms')
    attempt = 0
    if oldstate is None:
        recgen = reader(fid)
    try:
        for SNRs, HK in recgen:
            if yesterday:
                info("Done prepopulating. Switching to today's file.")
                fid.close()
                fid = open(current_binex(), 'rb')
                ofile = os.path.abspath(current_binex())
                recgen.send(fid)
                yesterday = False
                continue
            nlens = defaultdict(int, {rx: len(SNRs[rx]) for rx in SNRs}, hk=len(HK))
            tic = np.datetime64('now')
            if nlens == olens:
                attempt += 1
                if attempt > 6:
                    if handover or (fname is None and handover is not False and after_midnight()):
                        if ofile != os.path.abspath(current_binex()):
                            info('No new records at', tic, '. Attempting handover.')
                            fid.close()
                            info('Starting midnight plotting at', tic)
                            midnightplots(SNRs, HK)
                            info('Done at', np.datetime64('now'))
                            fid = open(current_binex(), 'rb')
                            ofile = os.path.abspath(current_binex())
                            recgen.send(fid)
                            attempt = 0
                            continue
                    info('No new records at', tic, '. Reporting.')
                    report_failure(rectic, dattic)
                    time.sleep(config.PLOT_IVAL / np.timedelta64(1, 's'))
                info('No new records at', tic, '. Sleeping', attempt*5)
                time.sleep(attempt*5)
                continue
            attempt = 0
            rectic = tic
            dattic = max(max(s[-1]['time'] for s in SNRs.values()), HK[-1]['time'])
            debug('{:2} new records {} at'.format(sum(nlens.values()) - sum(olens.values()),
                                                 [nlens[rx] - olens[rx] for rx in nlens]),
                  tic, 'timestamped', dattic)
            olens = nlens
            if tic - oldtic > config.PLOT_IVAL:
                info('Starting plotting at', tic)
                makeplots(SNRs, HK, endtime=tic, minelev=10)
                info('Done at', np.datetime64('now'))
                oldtic = tic
            else:
                time.sleep(2)
    except KeyboardInterrupt:
        return recgen, fid

def usage():
    print("Usage:", sys.argv[0], "<file> [snr_hours] [hk_hours] [endtime]\n"
          "  <file> should be the path to VAPR BINEX data. \n"
          "  snr_hours (optional): how many hours to plot SNR values, default", config.PLOT_SNR_HOURS, "\n"
          "  hk_hours (optional): how many hours to plot housekeeping values, default", config.PLOT_HK_HOURS, "\n"
          "  endtime (optional): plots are made preceding this time (default now)\n"
          "                      Format 2017-03-15T17:06:59\n")
    sys.exit(1)

if __name__ == "__main__": # When this file is run as a script
    if len(sys.argv) < 2 or len(sys.argv) > 5:
        usage()
    snrhours = config.PLOT_SNR_HOURS
    hkhours = config.PLOT_HK_HOURS
    endtime = np.datetime64('now')
    if len(sys.argv) > 2:
        try:
            snrhours = float(sys.argv[2])
        except ValueError:
            print('Hours parameter must be a number')
            usage()
    if len(sys.argv) > 3:
        try:
            hkhours = float(sys.argv[3])
        except ValueError:
            print('Hours parameter must be a number')
            usage()
    if len(sys.argv) > 4:
        try:
            endtime = np.datetime64(sys.argv[4])
        except ValueError:
            print('Endtime must be given in the format YYYY-MM-DDTHH:MM:SS')
            usage()
    with open(sys.argv[1], 'rb') as fid:
        SNRs, HK = readall(fid)
    makeplots(SNRs, HK, snrhours=snrhours, hkhours=hkhours, endtime=endtime)
