from math import floor, ceil, pi
import numpy as np
import matplotlib as mp
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D #analysis:ignore
from ascii_read import gpssowgps

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
            print("{} (entry {}) is positive, after negative entry at {}!".format(k, i, ind))
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
        plt.ylabel('SNR')
        plt.title('PRN {}, DoY {}, {:02d}:{:02d}--{:02d}:{:02d}, Az {}--{}'.format(
            prn, snrdata.doy, int(sod[rb]//3600), int(sod[rb] % 3600)//60,
            int(sod[re]//3600), int(sod[re] % 3600)//60,
            min(az[rb:re+1]), max(az[rb:re+1])))
    plt.tight_layout()

def dorises2(snrdata, prn):
    """Plot snr vs. elevation for all periods of rising-then-falling elevation.

    Plot each ascension in a different figure, and save to a file instead of
    displaying. Input is a numpy recarray of SNR records. One image is generated
    per rising/falling arc, named as PRN-DOY-HR-Q#.png (HR is the time the arc
    begins, and # is the quadrant wherein the satellite rose.)
    """
    plt.ioff()
    snrdata = snrdata[snrdata.prn == prn]
    if not snrdata.size:
        print('PRN {} not found'.format(prn))
        return
    eles, azis, snr = snrdata.el, snrdata.az, snrdata.snr
    sow = gpssowgps(snrdata.time)
    riz = rises(eles, sow, prn)
    for beg, peak, end in riz:
        fig = plt.figure(figsize=(16,6))
        gs = mp.gridspec.GridSpec(1, 2, width_ratios=[4,1])
        ax0 = plt.subplot(gs[0])
        pel = eles[peak]
        twoel = 2*pel
        ax0.scatter(eles[beg:peak+1], snr[beg:peak+1], s=1)
        ax0.scatter(twoel - eles[peak+1:end+1], snr[peak+1:end+1], s=1)
        xt = np.arange(pel/5, twoel-1, pel/5)
        xlab = ['{:.0f}°'.format(x if x <= pel else twoel - x) for x in xt]
        plt.xticks(xt, xlab)
        plt.xlim(floor(eles[beg]), ceil(twoel - eles[end]))
        plt.xlabel('Elevation')
        plt.ylabel('SNR')
        doy = snrdata[beg].time.tolist().strftime('%j')
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
            print('Azimuths are not monotone.')
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
            print('Less than 10 minutes, PRN {}, {} to {} ({}--{})'.format(
                prn, beg+1, end, sowdhrmin(sod[beg+1]), sowhrmin(sod[end])))
            continue
        peak = posneg(difel[beg+1:end])
        if peak == 0: # only falling elevations I guess
            print('Only falling elevations? PRN {}, {} to {} ({}--{})'.format(
                prn, beg+1, end, sowdhrmin(sod[beg+1]), sowhrmin(sod[end])))
            continue
        if peak is not None:
            peak += beg + 1
        else:
            peak = end
        if el[peak] < 20:
            print('Max elevation {}, PRN {}, {} to {} ({}--{})'.format(
                el[peak], prn, beg+1, end, sowdhrmin(sod[beg+1]), sowhrmin(sod[end])))
            continue
        riz.append([beg+1, peak, end])
    return riz

def _mode(arr):
    """Return the value occuring most often in arr."""
    vals, ct = np.unique(arr, return_counts=True)
    return vals[np.argmax(ct)]

def _gethouraxes(figsize, **kwargs):
    """Return figure and axes set up for plotting against time, labeled HH:MM."""
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(1, 1, 1, **kwargs)
    ax.xaxis_date()
    ax.xaxis.set_major_formatter(mp.dates.DateFormatter('%H:%M'))
    ax.set_xlabel('Time (UTC)')
    return fig, ax

def _thresh(hrs, endtime):
    if endtime is None:
        endtime = np.datetime64('now')
    if hrs is not None:
        return endtime - hrs * np.timedelta64(3600, 's')


def tempvolt(hk, hrs=None, endtime=None):
    """Plot temperature and voltage for each receiver in a recarray of housekeeping records.

    The two plots are on the same axes, against time. One image per receiver is
    written out in the current directory, named TV-RX##-DOY.png.
    If hrs is given, hrs hours preceding the datetime64 endtime (default now)
    are plotted (otherwise, all data is plotted.)
    """
    plt.ioff()
    thresh = _thresh(hrs, endtime)
    for rx in np.unique(hk.rxid):
        mask = hk.rxid == rx
        if hrs is not None:
            mask = np.logical_and(mask, hk.time > thresh)
        times = hk[mask].time.tolist()
        volt = hk[mask].volt / 100
        temp = hk[mask].temp
        doy = _mode(hk[mask].time.astype('M8[D]')).tolist()
        title = 'Rx{:02} {:%Y-%m-%d}: Voltage and Temperature'.format(rx, doy)
        fig, ax = _gethouraxes((6, 3), title=title)
        ax.scatter(times, volt, c='b')
        ax.set_ylabel('Volts', color='b')
        ax2 = ax.twinx()
        ax2.scatter(times, temp, c='r')
        ax.set_xlim(min(times), max(times))
        ax2.set_ylabel('Temperature (°C)', color='r')
        fig.tight_layout()
        fig.savefig('TV-RX{:02}-{:%j}.png'.format(rx, doy))
        plt.close(fig)

def numsats(snr, rxid=None, hrs=None, endtime=None):
    """Plot number of tracked satellites vs. time from a recarray of SNR records.

    If rxid is given, an image is written out in the current directory, named
    NS-RX##-DOY.png. Otherwise the figure is returned.
    """
    if hrs is not None:
        thresh = _thresh(hrs, endtime)
        snr = snr[snr.time > thresh]
    time, numsats = np.unique(snr.time.astype('M8[s]'), return_counts=True) # round to seconds
    doy = _mode(time.astype('M8[D]')).tolist() # most common day
    time = time.tolist() # get datetime.datetime, which matplotlib can handle
    if rxid:
        plt.ioff()
        title = 'Rx{:02} {:%Y-%m-%d}: Number of satellites'.format(rxid, doy)
    else:
        title = 'VAPR {:%Y-%m-%d}: Number of satellites'.format(doy)
    fig, ax = _gethouraxes((6, 3), title=title, ylabel='Tracked satellites')
    ax.plot(time, numsats, '.')
    ax.set_xlim(min(time), max(time))
    ax.set_ylim(min(numsats) - 1, max(numsats)+1)
    fig.tight_layout()
    if rxid:
        fig.savefig('NS-RX{:02}-{:%j}.png'.format(rxid, doy))
        plt.close(fig)
    else:
        return fig

def meansnr(snr, rxid=None, hrs=None, endtime=None):
    """Plot mean snr vs. time from a recarray of SNR records.

    If rxid is given, an image is written out in the current directory, named
    AVG-RX##-DOY.png. Otherwise the figure is returned.
    """
    if hrs is not None:
        thresh = _thresh(hrs, endtime)
        snr = snr[snr.time > thresh]
    time, idx = np.unique(snr.time.astype('M8[s]'), return_index=True) # round to seconds
    doy = _mode(time.astype('M8[D]')).tolist() # most common day
    time = time.tolist() # get datetime.datetime, which matplotlib can handle
    idx = np.append(idx, len(snr))
    means = []
    for a, b in zip(idx, idx[1:]):
        means += [np.mean(snr[a:b].snr)/100]
    if rxid:
        plt.ioff()
        title = 'Rx{:02} {:%Y-%m-%d}: Mean SNR'.format(rxid, doy)
    else:
        title = 'VAPR {:%Y-%m-%d}: Mean SNR'.format(doy)
    fig, ax = _gethouraxes((6, 3), title=title, ylabel='Mean SNR')
    ax.scatter(time, means)
    ax.set_xlim(min(time), max(time))
    fig.tight_layout()
    if rxid:
        fig.savefig('AVG-RX{:02}-{:%j}.png'.format(rxid, doy))
        plt.close(fig)
    else:
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

def polarazel(azis, eles, ax=None):
    """Plot azimuth and elevations (in degrees) as a curve on a polar plot.

    Input should be numpy arrays."""
    ax = _getaxis(ax, 'polar')
    line = ax.plot((90 - azis)/180 * pi, 90 - eles)
    xloc, _ = plt.xticks()
    xlab = (90 - xloc * 180 / pi) % 360
    plt.xticks(xloc, ['{:.0f}°'.format(x) for x in xlab])
    yloc = range(20, 90, 20)
    plt.yticks(yloc, (str(90 - y) for y in yloc))
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
    base = np.array(base).reshape(3,1,1)
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
    ax.plot_surface(gloc[:,:,0], gloc[:,:,1], gloc[:,:,2],
                           facecolors=mp.cm.terrain(N), rstride=stride, cstride=stride)
    vmark = np.array(base) * 1.00005
    ax.scatter(*vmark, marker='^', s=100, c='r')
    # the above draws a triangle at the vent location, but the surface usually obscures it.
    # vmark is multiplied by 1.00005 to lift it above the surface a little,
    # but it still doesn't show up reliably.
    plotcylinder(base, uax, length, radius, ax, color='r')
    NR = np.array([[gloc[:,:,i].min(), gloc[:,:,i].max()] for i in [0,1,2]])
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
        RR[i,0] -= ext/2
        RR[i,1] += ext/2
    ax.auto_scale_xyz(*RR)
    ax.set_aspect('equal')
