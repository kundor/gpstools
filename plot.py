from math import floor, ceil, pi
import os
from contextlib import suppress
import time
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
import matplotlib as mp
mp.use('Agg', warn=False) # non-interactive backend, allows importing without crashing X-less ipython
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D #analysis:ignore
from gpstime import gpssowgps
from utility import info, mode
from findfirst import findfirstgt
import config

def posneg(arr):
    """Check that input consists of some positive entries followed by negative entries,
    with no mixing. Return index of first negative entry, or None."""
    neg = False
    ind = None
    for i, k in enumerate(arr):
        if not neg and k < 0:
            ind = i
            neg = True
        if neg and k > 0:
            info("{} (entry {}) is positive, after negative entry at {}!".format(k, i, ind))
            return ind
    return ind

def dorises(snrdata, prn):
    """Plot snr vs. elevation for all periods of rising elevation.

    Finds all ascensions of satellite prn in snrdata, and plots them
    in subplots of a single figure.
    """
    _, el, az, sod, snr = zip(*(r for r in snrdata if r.prn == prn))
    riz = rises(el, sod)
    for i, (rb, re, _) in enumerate(riz):
        plt.subplot(len(riz), 1, i+1)
        plt.scatter(el[rb:re+1], snr[rb:re+1], s=1)
        plt.xlim(el[rb], el[re])
        plt.xlabel('Elevation')
        plt.ylabel('SNR (dB-Hz)')
        plt.title('PRN {}, DoY {}, {:02d}:{:02d}--{:02d}:{:02d}, Az {}--{}'.format(
            prn, snrdata.doy, int(sod[rb]//3600), int(sod[rb] % 3600)//60,
            int(sod[re]//3600), int(sod[re] % 3600)//60,
            min(az[rb:re+1]), max(az[rb:re+1])))
    plt.tight_layout()

def _setsnrlim(ax, snrs, restrict=False):
    y0, y1 = config.SNR_RANGE
    if not np.count_nonzero(snrs > 0):
        ax.set_ylim(y0, y1)
        return
    s0, s1 = min(snrs[snrs > 0]), max(snrs)
    if restrict:
        y0 = max(y0, s0)
        y1 = min(y1, s1)
    ax.set_ylim(y0, y1)
    if s0 < y0 or s1 > y1:
        s0, s1, y0, y1 = [round(x, 1) for x in (s0, s1, y0, y1)]
        info('Warning: SNRs are off the plot: {} in [{}, {}), {} in ({}, {}]'.format(
            np.count_nonzero(snrs[snrs > 0] < y0), s0, y0,
            np.count_nonzero(snrs > y1), s1, y1))

def dorises2(snrdata, prn):
    """Plot snr vs. elevation for all periods of rising-then-falling elevation.

    Plot each ascension in a different figure, and save to a file instead of
    displaying. Input is a numpy array of SNR records. One image is generated
    per rising/falling arc, named as PRN-DOY-HR-Q#.png (HR is the time the arc
    begins, and # is the quadrant wherein the satellite rose.)
    """
    plt.ioff()
    snrdata = snrdata[snrdata['prn'] == prn]
    if not snrdata.size:
        info('PRN {} not found'.format(prn))
        return
    eles, azis, snr = snrdata['el'], snrdata['az'], snrdata['snr'] / 10
    sow = gpssowgps(snrdata['time'])
    riz = rises(eles, sow, prn)
    for beg, peak, end in riz:
        if not np.count_nonzero(snr[beg:end+1]):
            info('Only 0 snr, PRN {}, {} to {} ({}--{})'.format(
                prn, beg, end, sowdhrmin(sow[beg]), sowhrmin(sow[end])))
            continue
        fig = plt.figure(figsize=(12, 4))
        gs = mp.gridspec.GridSpec(1, 2, width_ratios=[4, 1])
        ax0 = plt.subplot(gs[0])
        pel = eles[peak]
        twoel = 2*pel
        ax0.plot(eles[beg:peak+1], snr[beg:peak+1], 'k.', ms=1)
        ax0.plot(twoel - eles[peak+1:end+1], snr[peak+1:end+1], 'k.', ms=1)
        ax0.axvline(eles[peak])
        xt = np.arange(pel/5, twoel-1, pel/5)
        xlab = ['{:.0f}°'.format(x if x <= pel else twoel - x) for x in xt]
        plt.xticks(xt, xlab)
        plt.xlim(floor(eles[beg]), ceil(twoel - eles[end]))
        _setsnrlim(ax0, snr[beg:end+1])
        plt.xlabel('Elevation')
        plt.ylabel('SNR (dB-Hz)')
        doy = snrdata[beg]['time'].tolist().strftime('%j')
        daz = np.diff(azis[beg:end+1])
        if (daz >= 0).all():
            azinc = r'$\uparrow$'
        elif (daz <= 0).all():
            azinc = r'$\downarrow$'
        elif np.count_nonzero(daz < 0) == 1 and daz[daz < 0] < -350:
            azinc = r'$\uparrow$'
        elif np.count_nonzero(daz > 0) == 1 and daz[daz > 0] > 350:
            azinc = r'$\downarrow$'
        else:
            info('Azimuths are not monotone.')
            azinc = ''
        plt.title('PRN {}, DoY {}, {}--{}, Az {:.0f}--{:.0f} {}'.format(
            prn, doy, sowhrmin(sow[beg]), sowhrmin(sow[peak]), azis[beg], azis[end], azinc))
        ax1 = plt.subplot(gs[1], projection='polar')
        polarazel(azis[beg:end+1], eles[beg:end+1], ax1)
        plt.tight_layout()
        quart = int(5 - azis[beg] // 90)  # rising quarter
        if quart == 5:
            quart = 1
        hr = (int(sow[beg]) % 86400) // 3600
        plt.savefig('{:02d}-{}-{:02d}-Q{}.png'.format(prn, doy, hr, quart))
        plt.close(fig)

def allrises(SNR):
    """Plot all rising/falling arcs for each PRN found for each rxid in SNR.

    Plots are saved to images in subdirectory named rx## under the current directory.
    """
    for rxid in SNR:
        rdir = 'rx{:02}'.format(rxid)
        with suppress(FileExistsError):
            os.mkdir(rdir)
        os.chdir(rdir)
        for prn in np.unique(SNR[rxid]['prn']):
            dorises2(SNR[rxid], prn)
        os.chdir(os.pardir)

def sodhrmin(sod):
    return '{:02d}:{:02d}'.format(sod//3600, (sod % 3600)//60)

def sowhrmin(sow):
    return sodhrmin(int(sow) % 86400)

def sowdhrmin(sow):
    return str(int(sow) // 86400) + ';' + sowhrmin(sow)

def rises(el, sod, prn=None):
    difel = np.diff(el)
    starts = [-1] + np.argwhere(np.diff(sod) > 1000).ravel().tolist() + [len(difel)]
    riz = []
    for beg, end in zip(starts, starts[1:]):
        if sod[end] - sod[beg+1] < 600: # less than 10 minute arc
            info('Less than 10 minutes, PRN {}, {} to {} ({}--{})'.format(
                prn, beg+1, end, sowdhrmin(sod[beg+1]), sowhrmin(sod[end])))
            continue
        peak = posneg(difel[beg+1:end])
        if peak == 0: # only falling elevations I guess
            info('Only falling elevations? PRN {}, {} to {} ({}--{})'.format(
                prn, beg+1, end, sowdhrmin(sod[beg+1]), sowhrmin(sod[end])))
            continue
        if peak is not None:
            peak += beg + 1
        else:
            peak = end
        if el[peak] < 20:
            info('Max elevation {}, PRN {}, {} to {} ({}--{})'.format(
                el[peak], prn, beg+1, end, sowdhrmin(sod[beg+1]), sowhrmin(sod[end])))
            continue
        riz.append([beg+1, peak, end])
    return riz

class UTCLocalFormat:
    """Format datetime tick labels as two lines, UTC HH:MM on top, local HH:MM on bottom."""
    def __init__(self, ax):
        """The last tick, determined by ax (the axis we are ticking), has timezones appended."""
        self.ax = ax
    def __call__(self, x, pos):
        dt = mp.dates.num2date(x)
        if pos is None:
            return dt.strftime('%H:%M')
        if pos == len(self.ax.major.locator()) - 1:
            return dt.strftime('%H:%M (UTC)\n')  + dt.astimezone().strftime('%H:%M (%Z)')
        return dt.strftime('%H:%M\n') + dt.astimezone().strftime('%H:%M')

def _gethouraxes(figsize, pos, shareax=None, **kwargs):
    """Return figure and axes set up for plotting against time, labeled HH:MM."""
    if 'top' not in pos:
        figsize[1] -= 0.3
    if 'bot' not in pos:
        figsize[1] -= 0.4
    fig = plt.figure(figsize=figsize)
    if 'top' not in pos and 'title' in kwargs:
        del kwargs['title']
    if shareax is None:
        ax = fig.add_subplot(1, 1, 1, **kwargs)
        ax.xaxis_date()
    else:
        ax = fig.add_subplot(1, 1, 1, sharex=shareax, sharey=shareax, **kwargs)
    #ax.set_xlabel('Time (UTC)') # To save space on stacked plots: we'll only label the bottom x-axis
    bot = 0.06 / figsize[1]
    height = 1 - 0.15 / figsize[1]
    if 'bot' in pos:
        bot = 0.42 / figsize[1]
        height -= bot
        ax.xaxis.set_major_formatter(mp.ticker.FuncFormatter(UTCLocalFormat(ax.xaxis)))
    else:
        ax.xaxis.set_major_formatter(mp.ticker.NullFormatter())
    if 'top' in pos:
        height -= 0.18 / figsize[1]
    ax.set_position([0.055, bot, 0.89, height]) # ensure the plots line up with eachother
    return fig, ax

def _thresh(hrs, endtime):
    if endtime is None:
        endtime = np.datetime64('now')
    if hrs is not None:
        return endtime - hrs * np.timedelta64(3600, 's'), endtime
    return None, None

def prn_snr(SNR, rxid=None, snr_hours=None, endtime=None, omit_zero=True, minelev=0, figlength=12,
            doazel=True, suffix='', snrminel=None, **unused):
    """Plot SNR for each tracked satellite over the time period given.

    If rxid is given, an image is written out in the current directory, named
    ALLSNR-RX##-DOY.png. Otherwise the figure is returned.
    If snr_hours is set to 'all', then all data in the provided array are plotted.
    If snr_hours is not given, then config.PLOT_SNR_HOURS is used.
    """
    diftime = np.timedelta64(1, 'D') - np.timedelta64(246, 's')
    cur, prev = 'Current', 'Previous day'
    if snr_hours is None:
        snr_hours = config.PLOT_SNR_HOURS
    if snr_hours == 'all':
        thresh = None
    else:
        thresh, endtime = _thresh(snr_hours, endtime)
        if endtime > np.datetime64('now') - np.timedelta64(1, 'h'):
            cur, prev = 'Today', 'Yesterday'
    if snrminel is not None:
        minelev = snrminel
    if minelev:
        SNR = SNR[SNR['el'] > minelev]
    if omit_zero:
        SNR = SNR[SNR['snr'] > 0]
    if thresh:
        OSNR = SNR[findfirstgt(SNR['time'], thresh - diftime) :
                   findfirstgt(SNR['time'], endtime - diftime)]
        SNR = SNR[findfirstgt(SNR['time'], thresh) : findfirstgt(SNR['time'], endtime)]
    uprns, ct = np.unique(SNR['prn'], return_counts=True)
    prns = uprns[ct > 1200] # at least 20 minutes data
    numsat = len(prns)
    if not numsat:
        # relax our standards and try again
        prns = uprns
        numsat = len(prns)
    if not numsat:
        info('No records', 'for RX{:02}'.format(rxid) if rxid else '',
             'in the given time period', thresh, 'to', endtime)
        return
    fig = plt.figure(figsize=(figlength, 3*numsat))
    gs = mp.gridspec.GridSpec(numsat, 1 + doazel, width_ratios=[figlength - 2, 1.8*doazel],
                              left=0.53/figlength, right=1-(0.5-0.272*doazel)/figlength,
                              bottom=0.08/numsat, top=1-0.02/numsat,
                              wspace=0.12, hspace=0.09)
    axes = [0]*numsat
    dax = None
    minsnr = 100.
    maxsnr = 0.
    for i, prn in enumerate(prns):
        psnr = SNR[SNR['prn'] == prn]
        gs2 = mp.gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[i, 0], hspace=0., height_ratios=[2, 1])
        if i:
            ax = axes[i] = fig.add_subplot(gs2[0], sharex=axes[0], sharey=axes[0])
        else:
            ax = axes[0] = fig.add_subplot(gs2[0])
            ax.xaxis_date()
            ax.xaxis.set_major_formatter(mp.dates.DateFormatter('%H:%M'))
        ax.plot(psnr['time'].tolist(), psnr['snr'] / 10, 'k.', ms=2, label=cur)
        rxlab = 'RX{:02}: '.format(rxid) if rxid else ''
        ax.set_ylabel(rxlab + 'PRN {:02}'.format(prn))
        minsnr = min(minsnr, np.percentile(psnr['snr'], 0.02) / 10)
        maxsnr = max(maxsnr, np.percentile(psnr['snr'], 99.99) / 10)
        if thresh:
            opsnr = OSNR[OSNR['prn'] == prn]
            dax = plot_old_snr(fig, gs2, ax, psnr, opsnr, diftime, prev, dax)
        if doazel:
            ax1 = fig.add_subplot(gs[i, 1], projection='polar')
            polarazel(psnr['az'], psnr['el'], ax1, label_el=False)
    axes[0].set_ylim(minsnr, maxsnr)
    if thresh:
        axes[0].set_xlim(thresh.tolist(), endtime.tolist())
        if len(OSNR):
            axes[0].legend(loc='best', fontsize='small', numpoints=1, markerscale=4,
                           handletextpad=0.4, handlelength=0.8)
    else:
        axes[0].set_xlim(min(SNR['time'].tolist()), max(SNR['time'].tolist()))
    #fig.tight_layout(pad=0.1, h_pad=0.01, w_pad=0.01)
    if rxid:
        fname = 'ALLSNR-RX{:02}{}.png'.format(rxid, suffix)
        fig.savefig(fname)
        plt.close(fig)
        return fname
    return fig

def _twinax(ax, **kwargs):
    """Like ax.twinx, but you can specify keyword arguments for the new axis."""
    ax2 = ax._make_twin_axes(sharex=ax, **kwargs)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position('right')
    ax2.yaxis.set_offset_position('right')
    ax.yaxis.tick_left()
    ax2.xaxis.set_visible(False)
    ax2.patch.set_visible(False)
    ax2.set_position(ax.get_position())
    return ax2

def plot_old_snr(fig, gs, ax, psnr, opsnr, diftime, prev, dax):
    window = 31
    if opsnr.size < window or psnr.size < window:
        return dax
    tims = opsnr['time'] + diftime
    today = pd.Series(psnr['snr'] / 10, psnr['time'])
    yestr = pd.Series(opsnr['snr'] / 10, tims)
    #mn_dif = today.resample('5s').mean() - yestr.resample('5s').mean()
    sm_td = pd.Series(savgol_filter(today, window, 4), psnr['time'])
    sm_yd = pd.Series(savgol_filter(yestr, window, 4), tims)
    sm_dif = sm_td.resample('1s').first() - sm_yd.resample('1s').first()
    els = pd.Series(psnr['el'], psnr['time'])
    ax.plot(yestr, 'g.', ms=2, zorder=0, label=prev)
    plt.setp(ax.get_xticklabels(), visible=False)
    m0 = np.floor(np.nanpercentile(sm_dif, 0.05) * 10) / 10
    m1 = np.ceil(np.nanpercentile(sm_dif, 99.95) * 10) / 10
    if dax:
        y0, y1 = dax.get_ylim()
        m0 = min(m0, y0)
        m1 = max(m1, y1)
    m0 = max(m0, -10)
    m1 = min(m1, 10)
    ax2 = fig.add_subplot(gs[1], sharex=ax, sharey=dax)
    ax2.yaxis.set_major_locator(mp.ticker.MaxNLocator(nbins=4, prune='both', integer=True))
    ax2.fill_between(sm_dif.index, 0, sm_dif, color='m', label='Difference')
    ax2.set_ylabel('Difference', labelpad=2, color='m')
    ax2.tick_params('y', colors='m')
    ax2.set_ylim(m0, m1)
    ax3 = _twinax(ax2)
    ax3.plot(els, 'k', label='Elevation')
    ax3.set_ylabel('Elevation', labelpad=2)
    ax3.yaxis.set_major_locator(mp.ticker.MultipleLocator(20))
    ax3.set_ylim(10, 90)
    return ax2

def _expandlim(minmax, ax):
    """Expand y-range to the given (y0, y1) (but don't shrink it.)"""
    y0, y1 = minmax
    a0, a1 = ax.get_ylim()
    if a0 <= y0 and a1 >= y1:
        return
    ax.set_ylim(min(a0, y0), max(a1, y1))

def _localmidnight(ax):
    """Put a line at local midnight on the x-axis.

    X-axis should be labeled with matplotlib-style UTC dates, and is assumed
    to span about 24 hours."""
    x0, x1 = ax.get_xlim()
    utc_hr = -(time.localtime().tm_gmtoff / 60 / 60) # local midnight in UTC hours
    dfrac = utc_hr / 24 # fractional day when local midnight occurs
    if dfrac > (x0 % 1):
        midnt = floor(x0) + dfrac
    elif dfrac < (x1 % 1):
        midnt = floor(x1) + dfrac
    elif floor(x1) - floor(x0) > 1:
        midnt = floor(x0) + 1 + dfrac
    else:
        info('Local midnight did not occur in axis')
        return
    ax.axvline(midnt, color='0.3', zorder=0)

def tempvolt(hk, pos, shareax=None):
    """Plot temperature and voltage from the array of HK records.

    Return a figure.
    """
    times = hk['time'].tolist()
    volt = hk['volt'] / 100
    temp = hk['temp']
    fig, ax = _gethouraxes([10, 3], pos, shareax)
    ax.plot(times, volt, 'b.', ms=6)
    ax.set_ylabel('Volts', labelpad=3, color='b')
    ax.tick_params('y', colors='b')
    ax.yaxis.set_major_formatter(mp.ticker.FormatStrFormatter('%.1f'))
    _expandlim(config.VOLT_RANGE, ax)
    if min(np.diff(ax.yaxis.get_major_locator()())) < 0.1:
        ax.yaxis.set_major_locator(mp.ticker.MultipleLocator(0.1))
    ax.grid(True)
    ax2 = _twinax(ax)
    ax2.plot(times, temp, c='r')
    ax.set_xlim(min(times), max(times))
    _localmidnight(ax)
    ax2.set_ylabel('Temperature (°C)', labelpad=4, color='r')
    ax2.tick_params('y', colors='r')
    _expandlim(config.TEMP_RANGE, ax2)
    return fig, ax

def tempvolts(hk, hk_hours=None, endtime=None, pos='topbot', suffix='', **unused):
    """Plot temperature and voltage for each receiver in a array of housekeeping records.

    The two plots are on the same axes, against time. One image per receiver is
    written out in the current directory, named TV-RX##-DOY.png.
    If hk_hours is given, hk_hours hours preceding the datetime64 endtime (default now)
    are plotted (otherwise, all data is plotted.)
    """
    plt.ioff()
    thresh, endtime = _thresh(hk_hours, endtime)
    fnames = []
    figs = []
    ax = None
    if hk_hours is not None:
        hk = hk[findfirstgt(hk['time'], thresh) : findfirstgt(hk['time'], endtime)]
    for rx in np.unique(hk['rxid']):
        mask = hk['rxid'] == rx
        if not mask.any():
            info('No records for RX{:02} in the given time period'.format(rx), thresh, 'to', endtime)
            continue
        fig, ax = tempvolt(hk[mask], pos, ax)
        doy = mode(hk[mask]['time'].astype('M8[D]')).tolist()
        title = 'Rx{:02} {:%Y-%m-%d}'.format(rx, doy)
        if pos != 'top':
            title += ': Voltage and Temperature'
        if 'top' in pos:
            ax.set_title(title)
        if hk_hours is not None:
            ax.set_xlim(thresh.tolist(), endtime.tolist())
        figs += [fig]
        fnames += ['TV-RX{:02}{}.png'.format(rx, suffix)]
    # Wait to save them, in case the shared axes are resized
    for fig, fname in zip(figs, fnames):
        fig.savefig(fname)
        plt.close(fig)
    return fnames

def numsats(snr, rxid=None, hk_hours=None, endtime=None, minelev=0, pos='topbot', suffix='', **unused):
    """Plot number of tracked satellites vs. time from a array of SNR records.

    If rxid is given, an image is written out in the current directory, named
    NS-RX##-DOY.png. Otherwise the figure is returned.
    """
    if hk_hours is not None:
        thresh, endtime = _thresh(hk_hours, endtime)
        snr = snr[findfirstgt(snr['time'], thresh) : findfirstgt(snr['time'], endtime)]
    if minelev:
        snr = snr[snr['el'] > minelev]
    if not len(snr):
        info('No records', 'for RX{:02}'.format(rxid) if rxid else '',
             'above {}° from'.format(minelev), thresh, 'to', endtime)
        return
    time, idx = np.unique(snr['time'], return_index=True)
    idx.sort() # if the times aren't in order for some reason
    numsats = [len(set(snr[a:b]['prn'])) for a, b in zip(idx, np.append(idx[1:], len(snr)))]
#   time, idx, numsats = np.unique(snr['time'], return_index=True, return_counts=True)
#    for n, a, b in zip(numsats, idx, np.append(idx[1:], len(snr))):
#        if len(set(snr[a:b]['prn'])) != n:
#            info('Same PRN in multiple records?')
#            debug(snr[a:b])
    doy = mode(time.astype('M8[D]')).tolist() # most common day
    time = snr['time'][idx].tolist() # get datetime.datetime, which matplotlib can handle
    if rxid:
        title = 'Rx{:02} {:%Y-%m-%d}'.format(rxid, doy)
    else:
        title = 'VAPR {:%Y-%m-%d}'.format(doy)
    if pos != 'top':
        title += ': Number of satellites over {}°'.format(minelev)
    if 'top' in pos:
        ylabl = 'Tracked satellites'
    else:
        ylabl = 'No. satellites over {}°'.format(minelev)
    fig, ax = _gethouraxes([10, 2.7], pos, title=title, ylabel=ylabl)
    ax.plot(time, numsats, '-o', ms=2)
    if hk_hours is not None:
        ax.set_xlim(thresh.tolist(), endtime.tolist())
    else:
        ax.set_xlim(min(time), max(time))
    ax.yaxis.set_major_locator(mp.ticker.MaxNLocator(integer=True))
    ax.set_ylim(min(numsats) - 1, max(numsats)+1)
    ax.tick_params(axis='y', labelleft='on', labelright='on')
    ax.grid(True)
    if rxid:
        fname = 'NS-RX{:02}{}.png'.format(rxid, suffix)
        fig.savefig(fname)
        plt.close(fig)
        return fname
    return fig

def meansnr(snr, rxid=None, hk_hours=None, endtime=None, minelev=0, pos='topbot', suffix='', **unused):
    """Plot mean snr vs. time from a array of SNR records.

    If rxid is given, an image is written out in the current directory, named
    AVG-RX##-DOY.png. Otherwise the figure is returned.
    """
    if hk_hours is not None:
        thresh, endtime = _thresh(hk_hours, endtime)
        snr = snr[findfirstgt(snr['time'], thresh) : findfirstgt(snr['time'], endtime)]
    if minelev:
        snr = snr[snr['el'] > minelev]
    if not len(snr):
        info('No records', 'for RX{:02}'.format(rxid) if rxid else '',
             'in the given time period', thresh, 'to', endtime)
        return
    time, idx = np.unique(snr['time'], return_index=True)
    idx.sort() # if the times aren't in order for some reason
    time = snr['time'][idx]
    doy = mode(time.astype('M8[D]')).tolist() # most common day
    idx = np.append(idx, len(snr))
    means = []
    pmeans = []
    for a, b in zip(idx, idx[1:]):
        tsnr = snr[a:b]['snr'] / 10
        means += [np.mean(tsnr)]
        if np.count_nonzero(tsnr):
            pmeans += [np.mean(tsnr[tsnr > 0])]
        else:
            pmeans += [0]
    if rxid:
        title = 'Rx{:02} {:%Y-%m-%d}: Mean SNR'.format(rxid, doy)
    else:
        title = 'VAPR {:%Y-%m-%d}: Mean SNR'.format(doy)
    if minelev:
        title += ' over {}°'.format(minelev)
    ylabl = 'Mean SNR (dB-Hz)'
    if 'top' not in pos and minelev:
        ylabl += ' over {}°'.format(minelev)
    fig, ax = _gethouraxes([10, 3], pos, title=title, ylabel=ylabl)
    pmeans = np.array(pmeans)
    ax.plot(time.tolist(), means, 'b.', ms=2, label='All values')
    ax.plot(time[pmeans > 0].tolist(), pmeans[pmeans > 0], 'r.', ms=2, label='Positive only')
    if hk_hours is not None:
        ax.set_xlim(thresh.tolist(), endtime.tolist())
    else:
        ax.set_xlim(min(time.tolist()), max(time.tolist()))
    _setsnrlim(ax, pmeans) # add argument True to restrict y-range to observed values
    ax.tick_params(axis='y', labelleft='on', labelright='on')
    ax.legend(loc='best', fontsize='small', numpoints=1, markerscale=4,
              handletextpad=0.4, handlelength=0.8)
    if rxid:
        fname = 'AVG-RX{:02}{}.png'.format(rxid, suffix)
        fig.savefig(fname)
        plt.close(fig)
        return fname
    return fig

def add_arrow(line, start_ind=None, direction='right', size=15, color=None):
    """Add an arrow to a line.

    line:       Line2D object
    position:   x-position of the arrow. If None, mean of xdata is taken
    direction:  'left' or 'right'
    size:       size of the arrow in fontsize points
    color:      if None, line color is taken.
    """
    # from http://stackoverflow.com/a/34018322/2132213
    if color is None:
        color = line.get_color()

    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if start_ind is None:
        start_ind = len(xdata) // 2
    if direction == 'right':
        end_ind = start_ind + 1
    else:
        end_ind = start_ind - 1
    if end_ind >= len(xdata):
        return

    line.axes.annotate('',
                       xytext=(xdata[start_ind], ydata[start_ind]),
                       xy=(xdata[end_ind], ydata[end_ind]),
                       arrowprops={'arrowstyle': '->', 'color': color},
                       size=size)

def _getaxis(ax, projection, figsize=None):
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 1, 1, projection=projection)
    return ax

def polarazel(azis, eles, ax=None, label_el=True):
    """Plot azimuth and elevations (in degrees) as a curve on a polar plot.

    Input should be numpy arrays."""
    ax = _getaxis(ax, 'polar')
    line = ax.plot((90 - azis)/180 * pi, 90 - eles)
    xloc, _ = plt.xticks()
    xlab = (90 - xloc * 180 / pi) % 360
    plt.xticks(xloc, ['{:.0f}°'.format(x) for x in xlab])
    yloc = range(20, 90, 20)
    ax.set_ylim(0, 90)
    if label_el:
        plt.yticks(yloc, (str(90 - y) for y in yloc))
    else:
        plt.yticks(yloc, [])
    add_arrow(line[0])

def snrVSel(snrdata, prn, secstart=0, secend=86400, color=None):
    snr = [r.snr for r in snrdata if r.prn == prn and secstart < r.sod < secend]
    el = [r.el for r in snrdata if r.prn == prn and secstart < r.sod < secend]
    plt.scatter(el, snr, s=1, color=color)
#    plt.xlim(min(el), max(el))
#    plt.xlabel('Elevation')
#    plt.ylabel('SNR')
#    plt.title('PRN {}, DoY {}'.format(prn, snrdata.doy))
    plt.tight_layout()

def iterSNRs(SNRs):
    for snr in SNRs:
        for rec in snr:
            yield rec.az, rec.el, rec.snr

def itergdo(gdo):
    for r in gdo.iterlist(obscode=('az', 'el', 'S1'), skip=True):
        for rec in r:
            if rec[0] is not None:
                yield rec

def azelbin(iterfn, dat, scale=2):
    snravg = np.zeros((360*scale, 90*scale))
    snrnum = np.zeros((360*scale, 90*scale))
    snrstd = np.zeros((360*scale, 90*scale))
    for az, el, snr in iterfn(dat):
        azi = floor(az*scale)
        eli = floor(el*scale)
        n = snrnum[azi, eli]
        snravg[azi, eli] = n/(n+1)*snravg[azi, eli] + 1/(n+1)*snr
        snrnum[azi, eli] += 1
    for az, el, snr in iterfn(dat):
        azi = floor(az*scale)
        eli = floor(el*scale)
        snrstd[azi, eli] += (snr - snravg[azi, eli])**2
    snrstd = np.sqrt(snrstd / np.where(snrnum > 0, snrnum, 1))
    snravg = np.ma.masked_where(snrnum == 0, snravg)
    snrstd = np.ma.masked_where(snrnum == 0, snrstd)
    plotazelbin(snravg, 'Mean SNR', scale)
    plotazelbin(snrstd, 'SNR standard deviation', scale)
    return snravg, snrnum, snrstd

def plotazelbin(dat, title, scale=2):
    plt.figure()
    vmin = np.percentile(dat.compressed(), 1)
    vmin = max(floor(vmin), np.min(dat))
    vmax = np.percentile(dat.compressed(), 99)
    vmax = min(ceil(vmax), np.max(dat))
    plt.pcolormesh(dat.T, vmin=vmin, vmax=vmax)
    axx = plt.gca()
    box = axx.get_position()
    xloc, xlab = plt.xticks()
    plt.xticks(xloc, xloc/scale)
    yloc, ylab = plt.yticks()
    plt.yticks(yloc, yloc/scale)
    plt.xlim(0, 360*scale)
    plt.ylim(10*scale, 90*scale)
    plt.title(title)
    plt.xlabel('Azimuth')
    plt.ylabel('Elevation')
    axc = plt.axes([box.x0 + box.width * 1.05, box.y0, 0.01, box.height])
    plt.colorbar(cax=axc)

def perpvecs(uvec):
    """Return two unit vectors, perpendicular to the given unit vector and eachother."""
    if uvec[0] == 0:
        n1 = np.array([1, 0, 0])
    else:
        n1 = np.array([-uvec[1], uvec[0], 0])
        n1 /= np.linalg.norm(n1)
    n2 = np.cross(uvec, n1)
    return n1, n2

def plotradial(base, uaxis, length, radfn, ax=None, **kwargs):
    """Plot a radial surface.

    The surface is at the radius determined by the function radfn
    around the line segment of the given length, starting from base, in the
    direction of the unit vector uaxis.
    radfn is given arrays of t and theta values.
    Plot onto the Axes3D object ax if given.
    """
    ax = _getaxis(ax, '3d')
    n1, n2 = perpvecs(uaxis)
    t = np.linspace(0, length, 100)
    theta = np.linspace(0, 2 * np.pi, 100)
    t, theta = np.meshgrid(t, theta)
    mouter = np.multiply.outer
    base = np.array(base).reshape(3, 1, 1)
    XYZ = base + mouter(uaxis, t) + radfn(t, theta) * mouter(n1, np.sin(theta)) + radfn(t, theta) * mouter(n2, np.cos(theta))
    ax.plot_surface(*XYZ, **kwargs)

def plotcylinder(base, uaxis, length, radius, ax=None, alpha=0.2, **kwargs):
    """Plot (the surface of) an open finite cylinder.

    The cylinder is at the given radius around the line segment of the given
    length, starting from base, in the direction of the unit vector uaxis.
    Plot onto the Axes3D object ax if given.
    """
    def radfn(t, theta):
        return radius
    plotradial(base, uaxis, length, radfn, ax, alpha=alpha, **kwargs)

def plotcone(base, uaxis, length, radius, ax=None, alpha=0.2, **kwargs):
    """Plot (the surface of) a finite section of a cone.

    The cone has its apex at base, and opens in the direction of the unit vector
    uaxis, for the given length. It opens to the given radius after 1 unit.
    Plot onto the Axes3D object ax if given.
    """
    def radfn(t, theta):
        return radius * t
    plotradial(base, uaxis, length, radfn, ax, alpha=alpha, **kwargs)

def draw_plume(gloc, cyl_inf, ax=None, stride=10):
    """Given an array of ECEF ground locations, and cylinder info (as from
    skytracks.cylinfo()), draw the ground surface and the cylinder.

    Stride defaults to 10, a coarse map; stride of 3 gives a finer map,
    at expense of responsiveness.
    """
    from coords import earthradius, xyz2lat, lengthlat, lengthlon
    ax = _getaxis(ax, '3d', (12, 12))
    base, uax, length, radius = cyl_inf
    N = gloc @ uax # distance upward from center of the earth, assuming that's how the cylinder is oriented
    rlat = xyz2lat(*base)
    msl = earthradius(rlat)
    N -= msl
    N /= N.max() # normalize to 0..1
    ax.plot_surface(gloc[:, :, 0], gloc[:, :, 1], gloc[:, :, 2],
                    facecolors=mp.cm.terrain(N), rstride=stride, cstride=stride)
    vmark = np.array(base) * 1.00005
    ax.scatter(*vmark, marker='^', s=100, c='r')
    # the above draws a triangle at the vent location, but the surface usually obscures it.
    # vmark is multiplied by 1.00005 to lift it above the surface a little,
    # but it still doesn't show up reliably.
    plotcylinder(base, uax, length, radius, ax, color='r')
    NR = np.array([[gloc[:, :, i].min(), gloc[:, :, i].max()] for i in [0, 1, 2]])
    ax.auto_scale_xyz(*NR)
    ax.set_aspect(lengthlon(rlat * 180 / pi) / lengthlat(rlat * 180 / pi))


def equalize3d(ax, maxi=True):
    """Set aspect ratio so all three axes have the same scale."""
    RR = np.stack((ax.get_xlim(), ax.get_ylim(), ax.get_zlim()))
    difr = np.diff(RR)
    if maxi:
        maxr = np.max(difr)
    else:
        maxr = np.min(difr)
    for i in range(3):
        ext = maxr - difr[i]
        RR[i, 0] -= ext/2
        RR[i, 1] += ext/2
    ax.auto_scale_xyz(*RR)
    ax.set_aspect('equal')
