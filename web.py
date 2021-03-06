#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 11:54:05 2017

@author: nima9589
"""
import os
import sys
import time
import shutil
import json
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import defaultdict
import numpy as np
import config
import plot
from snrstats import calcsnrstat
from parser import reader, current_binex, readtimes, hkreport, cleanhk, snr88str
from utility import info, debug, pushdir, static_vars, email_exc, fmt_timesite
from findfirst import findfirstclosed
from vents import VENTS

def guard_and_time(msg, fn, tic, *args, **kwargs):
    info('Starting', msg, 'at', tic)
    try:
        fn(*args, **kwargs)
    except Exception:
        email_exc()
    info('Done at', np.datetime64('now'))


class Task:
    """An intermittent task, to be run every *ival* seconds in the message loop."""
    def __init__(self, ival, descrip, func):
        self.ival = np.timedelta64(ival, 's') # if ival is already timedelta64, this just changes its units
        self.descrip = descrip
        self.func = func
        self.lastrun = np.datetime64(0, 'ms')

    def __call__(self, tic, *args, **kwargs):
        if (tic - self.lastrun) > self.ival:
            guard_and_time(self.descrip, self.func, tic, *args, **kwargs)
            self.lastrun = tic

def _move(srcs, suf):
    """Move each filename in *srcs* to the name obtained by removing *suf*."""
    for src in srcs:
        if not src:
            continue
        if src.count(suf) != 1:
            info('Suffix', suf, 'not found uniquely in filename', src)
            continue
        os.replace(src, src.replace(suf, '', 1))

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

def rxloc(rxid, HK):
    """Return longitude, latitude, altitude for the given receiver."""
    RHK = HK[HK['rxid'] == rxid]
    return (np.median(RHK['lon']) / 1e7,
            np.median(RHK['lat']) / 1e7,
            np.median(RHK['alt']) / 1000)

def rxline(rxid, HK):
    """A list entry for the web page for the given receiver."""
    lng, lat, alt = rxloc(rxid, HK)
    return """    <li>
        RX{id:02}: <a href="#RX{id:02}-HK">HK</a> <a href="#RX{id:02}-SNR">SNR</a>
        {lng:.5f}&deg;{ew}, {lat:.5f}&deg;{ns}, {alt:.3f}m.
        </li>
        """.format(id=rxid, lng=abs(lng), ew='EW'[lng < 0], lat=abs(lat), ns='NS'[lat < 0], alt=alt)

def rxhkplots(rxid):
    """HTML to include the housekeeping plots for RX *rxid*."""
    rxstr = 'RX{:02}'.format(rxid)
    return """<h3 class="rxhead"><a name="{0}-HK">{0}</a></h3>
        <img src="TV-{0}.png"><br>
        <img src="NS-{0}.png"><br>
        <img src="AVG-{0}.png"><br>
        <hr>
        """.format(rxstr)

def snrplots(rxid):
    """HTML to include the SNR plots for RX *rxid*."""
    rxstr = 'RX{:02}'.format(rxid)
    return """<h3 class="rxhead"><a name="{0}-SNR">{0}</a></h3>
        <img src="ALLSNR-{0}.png">
        """.format(rxstr)

def sitedatajson(HK, rxlist=None, site=None):
    if site is None:
        site = config.SITE
    if rxlist is None:
        rxlist = np.unique(HK['rxid'])
    updated = str(np.datetime64('now', 'ms')) + 'Z'
    with open('sitedata.json', 'wt') as fid:
        json.dump({'rx': [dict(zip(('lng', 'lat', 'alt'), rxloc(rx, HK)), rxid=int(rx)) for rx in rxlist],
                   'vents': VENTS[site], 'updated': updated},
                  fid, indent=1)

@static_vars(rxlist=[])
def plotspage(HK, endtime):
    """Write out the receiver list and image lists for the web page."""
    since = endtime - min(config.PLOT_HK_HOURS, config.PLOT_SNR_HOURS) * np.timedelta64(3600, 's')
    HK = cleanhk(HK, since=since)
    rxlist = sorted(np.unique(HK['rxid']))
    if rxlist == plotspage.rxlist and all(os.path.exists(f) for f in ('rxlist.html', 'plots.html', 'sitedata.json')):
        return
    with open('rxlist.html', 'wt', encoding='utf-8') as fid:
         for rxid in rxlist:
             fid.write(rxline(rxid, HK))
    with open('plots.html', 'wt', encoding='utf-8') as fid:
        for rxid in rxlist:
            fid.write(rxhkplots(rxid))
        for rxid in rxlist:
            fid.write(snrplots(rxid))
    sitedatajson(HK, rxlist)
    plotspage.rxlist = rxlist

def makeplots(SNRs, HK, domove=True, plotdir=None, **plotargs):
    suf = '-new' * domove
    plotargs = {'hk_hours': config.PLOT_HK_HOURS,
                'snr_hours': config.PLOT_SNR_HOURS,
                'endtime': np.datetime64('now'),
                'suffix': suf,
                'minelev': config.MINELEV,
                **plotargs}
    # Any values already in plotargs will override these values
    hkfile = 'HK' + suf + '.txt'
    tabfile = 'snrtab' + suf + '.html'
    day = plotargs['endtime']
    hkstart = findfirstclosed(HK['time'], day - np.timedelta64(1, 'D'), day)
    if hkstart == -1:
        hkstart = 0 # If there are no good dates, put everything in the report
    files = [hkfile, tabfile]
    if plotdir is None:
        plotdir = fmt_timesite(config.PLOTDIR, day)
    with pushdir(plotdir):
        noteupdate()
        files += plot.tempvolts(cleanhk(HK), pos='top', **plotargs)
        with open(hkfile, 'wt', encoding='utf-8') as fid:
            hkreport(HK[hkstart:], fid)
        snrtab = open(tabfile, 'wt', encoding='utf-8')
        for rxid, SNR in SNRs.items():
            snrtab.write(format_stats(rxid, *calcsnrstat(SNR['snr'] / 10)))
            files += [plot.prn_snr(SNR, rxid, **plotargs),
                      plot.numsats(SNR, rxid, pos='mid', **plotargs),
                      plot.meansnr(SNR, rxid, pos='bot', **plotargs)]
        snrtab.close()
        if domove:
            _move(files, suf)
        plotspage(HK, day)
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
    ddir = fmt_timesite(ddir, day)
    try:
        os.makedirs(ddir, exist_ok=True)
    except FileExistsError as e:
        info('Could not create directory', ddir, '\n    ', e)
        return
    etime = day + np.timedelta64(1, 'D')
    makeplots(SNRs, HK, domove=False, plotdir=ddir, snr_hours=24, hk_hours=24, endtime=etime,
              figlength=14, doazel=False, snrminel=15)
    index = os.path.join(fmt_timesite(config.PLOTDIR, day), 'index.html')
    try:
        shutil.copy(index, ddir)
    except OSError:
        info('Could not copy index.html from plot dir to daily dir.')
    ylink = os.path.join(ddir, 'yesterlink.html')
    yday = day - np.timedelta64(1, 'D')
    with open(ylink, 'wt') as fid:
        fid.write(yday.tolist().strftime('../../%Y/%j'))

def noteupdate():
    updstr = time.strftime('.<br>Updating (%H:%M:%S %Z)...')
    with open('updatetime.txt', 'at', encoding='utf-8') as fid:
        fid.write(updstr)

def updatetime(append=''):
    updstr = time.strftime('%b %d %Y %H:%M:%S UTC', time.gmtime()) + time.strftime(' (%H:%M:%S %Z)')
    with open('updatetime.txt', 'wt', encoding='utf-8') as fid:
        fid.write(updstr + append)

def report_failure(rectic, dattic, pdir=None):
    """Record data failure in updatetime.txt"""
    if pdir is None:
        pdir = fmt_timesite(config.PLOTDIR)
    stamp = str(dattic)[:21] # chop off after deciseconds
    failstr = '''<br>
        <span class="fail">No data received since {}, timestamped {}</span>
        '''.format(rectic, stamp)
    with pushdir(pdir):
        updatetime(failstr)

def after_midnight():
    """True if current time is within 15 minutes after UTC midnight."""
    now = np.datetime64('now')
    return np.timedelta64(0) < (now - now.astype('M8[D]')) < np.timedelta64(15, 'm')

def write_snr88(SNRs, rxfids, olens, ofile):
    """Write out new SNR records to snr88 files."""
    for rx in SNRs:
        if rx not in rxfids:
            rtimes = SNRs[rx]['time']
            rtimes = rtimes[rtimes < np.datetime64('now') + np.timedelta64(1, 'h')]
            rtimes = rtimes[rtimes > np.datetime64('2000-01-01')]
            epoch = rtimes[-1]
            rxfids[rx] = open(os.path.join(os.path.dirname(ofile),
                              'rx{:02}'.format(rx) + epoch.tolist().strftime('%j0.%y.snr88')), 'w')
        for rec in SNRs[rx][olens[rx]:]:
            if rec['snr'] and rec['el'] > 10:
                rxfids[rx].write(snr88str(rec))

def plotupdate(fname=None, handover=None, oldstate=None, write88=False):
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
        try:
            fid = open(current_binex(yesterday), 'rb')
        except FileNotFoundError:
            yesterday = False
            fid = open(current_binex(), 'rb')
    elif fname is None:
        fid = open(current_binex(), 'rb')
    else:
        fid = open(fname, 'rb')
    ofile = os.path.abspath(fid.name)
    olens = defaultdict(int)
    rxfids = {}
    oldtic = rectic = dattic = np.datetime64('2000-01-01', 'ms')
    attempt = 0
    if oldstate is None:
        recgen = reader(fid)
    try:
        for SNRs, HK in recgen:
            if yesterday:
                if write88:
                    write_snr88(SNRs, rxfids, olens, ofile)
                    for rx in list(rxfids):
                        rxfids[rx].close()
                        del rxfids[rx]
                    olens = defaultdict(int, {rx: len(SNRs[rx]) for rx in SNRs}, hk=len(HK))
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
                            guard_and_time('midnight plotting', midnightplots, tic, SNRs, HK)
                            fid = open(current_binex(), 'rb')
                            ofile = os.path.abspath(current_binex())
                            recgen.send(fid)
                            for rx in list(rxfids):
                                rxfids[rx].close()
                                del rxfids[rx]
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
            if SNRs:
                dattic = max(s[-1]['time'] for s in SNRs.values())
            if HK.size and HK[-1]['time'] > dattic:
                dattic = HK[-1]['time']
            debug('{:2} new records {} at'.format(sum(nlens.values()) - sum(olens.values()),
                                                 [nlens[rx] - olens[rx] for rx in nlens]),
                  tic, 'timestamped', dattic)
            if write88:
                write_snr88(SNRs, rxfids, olens, ofile)
            olens = nlens
            if tic - oldtic > config.PLOT_IVAL:
                guard_and_time('plotting', makeplots, tic, SNRs, HK, endtime=tic)
                oldtic = tic
            else:
                time.sleep(2)
    except KeyboardInterrupt:
        return recgen, fid

if __name__ == "__main__": # When this file is run as a script
    parser = ArgumentParser(description='Create plots for the VAPR web site.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('site', metavar='SITE',
                        help='which site to plot for')
    parser.add_argument('-s', dest='snr_hours', default=config.PLOT_SNR_HOURS, type=float,
                        help='how many hours to plot SNR values')
    parser.add_argument('-k', dest='hk_hours', default=config.PLOT_HK_HOURS, type=float,
                        help='how many hours to plot housekeeping values')
    parser.add_argument('-e', dest='endtime', default='now', type=np.datetime64,
                        help='plots are made preceding this time. '
                             'Use ISO 8601-type format, e.g. 2017-03-15T17:06:59')
    parser.add_argument('-m', dest='minelev', default=config.MINELEV, type=float,
                        help='Only use values from above this elevation')
    parser.add_argument('-d', dest='plotdir', default=config.PLOTDIR,
                        help='the directory in which to put plotted images')
    parser.add_argument('-y', dest='daily', action='store_true',
                        help='Create daily overview plots: set defaults SNR_HOURS=HK_HOURS=24, '
                             'PLOTDIR=' + config.DAILYDIR.replace('%','%%') + ', '
                             'ENDTIME=last UTC midnight. '
                             'Make SNR plots longer, with minimum elevation at least 15, '
                             'and remove az/el maps')
    parser.add_argument('-u', dest='update', action='store_true',
                        help='Run as a daemon, updating plots every'+str(config.PLOT_IVAL)
                            +'. Any other arguments, except SITE, are ignored')
    args = parser.parse_args()
    config.SITE = args.site
    if args.update:
        config.LOGFILE = open(fmt_timesite(config.LOGFILELOC), 'w')
        plotupdate(write88=True)
        sys.exit()
    if args.daily:
        parser.set_defaults(snr_hours=24, hk_hours=24, minelev=15, plotdir=config.DAILYDIR,
                            endtime=np.datetime64('now', 'D'))
        args = parser.parse_args()
        args.figlength = 14
        args.doazel = False
        args.snrminel = max(15, args.minelev)
        args.domove = False
    del args.daily

    dhour = max(args.hk_hours, args.snr_hours)
    start = args.endtime - dhour * np.timedelta64(3600, 's')
    args.plotdir = fmt_timesite(args.plotdir, start)

    dhour = max(args.hk_hours, args.snr_hours + 24)
    start = args.endtime - dhour * np.timedelta64(3600, 's')

    SNRs, HK = readtimes(start, args.endtime)
    makeplots(SNRs, HK, **vars(args))
