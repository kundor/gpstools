"""Functions to make azimuth/elevation skytracks of satellites passing over a
given point, or a heatmap of ground locations with many lines-of-sight to satellites
passing through a given cylinder.
Currently hardcoded w.r.t. times, file locations.
"""
# -*- coding: utf-8 -*-
from math import sqrt
import os
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos
from plot import polarazel
from coords import xyz2llh, earthnormal_xyz, lengthlat, lengthlon, WGS84
from satpos import myinterp, mvec, coeffs
import findfirst

ETNA0 = (4878146, 1309017, 3886374) # ECEF in meters
REDBT = (-2802350, -1444879, 5527495)

def quickloadsp3(sp3file):
    pl = np.loadtxt(sp3file, skiprows=22, usecols=(1, 2, 3), comments=['*', 'E', 'V'])
    # This assumes all 32 PRNs are present, in order, for each epoch
    # Also assuming interval is 900 seconds (15 minutes)

    pl.shape = (-1, 32, 3) # reshape to 96x32x3 (-1 means 'fit the data')
    pl *= 1000
    return pl

def skytracks(sp3file, LOC):
    pl = quickloadsp3(sp3file)
    lat, lon, ht = xyz2llh(LOC)
    trans = np.array([[-sin(lon), -sin(lat)*cos(lon), cos(lat)*cos(lon)],
                      [ cos(lon), -sin(lat)*sin(lon), cos(lat)*sin(lon)],
                      [ 0,         cos(lat),          sin(lat)]])

    plt.figure()
    ax = plt.subplot(1, 1, 1, projection='polar')

    for prn in range(32): # prn is one less than the actual PRN :\
        ENU = (pl[:, prn, :] - LOC) @ trans
        az = np.arctan2(ENU[:, 0], ENU[:, 1]) * 180 / np.pi
        az[az < 0] += 360
        el = np.arctan2(ENU[:, 2], np.linalg.norm(ENU[:, :2], axis=1)) * 180 / np.pi
        ind = np.argwhere(el > 5).ravel()
        starts = [-1] + np.argwhere(np.diff(ind) > 1).ravel().tolist() + [len(ind) - 1]
        starts = np.array(starts) + 1
        for beg, end in zip(starts, starts[1:]):
            if end - beg < 3:
                continue
            polarazel(az[ind[beg:end]], el[ind[beg:end]], ax)

# skytracks('/scratch/sp3/igs19162.sp3', REDBT)

def quadformula(a, b, c):
    """Solve a x^2 + b x + c == 0 (real roots)."""
    if a == b == c == 0:
        raise ValueError("Infinite solutions")
    elif a == b == 0:
        return []
    elif a == 0:
        return [-c/b]
    a *= 2
    d = b*b - 2*a*c
    if d < 0:
        return []
    bb = -b / a
    if d > 0:
        dd = sqrt(d) / a
        return [bb - dd, bb + dd]
    return [bb]

def _clip(vals, t0, t1):
    """Ensure max(0, t0) <= vals[0] < vals[1] <= min(1, t1)"""
    if len(vals) != 2:
        return False
    assert vals[1] > vals[0]
    t0 = max(0, t0)
    t1 = min(1, t1)
    if vals[1] <= t0 or vals[0] >= t1:
        return False
    if vals[0] < t0 < vals[1]:
        vals[0] = t0
    if vals[0] < t1 < vals[1]:
        vals[1] = t1
    return True

def _lci(eproj, eperp, C, dif, uax, length):
    # Parameterize the line by L(t) = end0 + t*dif
    dproj = dif @ uax
    dperp = dif - dproj*uax
    A = dperp @ dperp
    B = 2 * eperp @ dperp
    slns = quadformula(A, B, C)
    t0 = -eproj / dproj # where line meets bottom of cylinder
    t1 = (length - eproj) / dproj # where line meets top of cylinder
    if not _clip(slns, t0, t1):
        return np.zeros((0, 3)) # empty array
    return np.outer(slns, dif)

def linecylinderintersect(end0, end1, cyl_inf):
    """Find intersections of a line segment with a finite cylinder.

    Return points where the line segment from end0 to end1 passes through the
    cylinder around the given axis, taken as a vector starting from base. The
    height of the cylinder is the length of axis. The radius of the cylinder
    is given. If no intersections exist, return nothing. If one or both endpoints
    is within the cylinder, return the endpoint(s).

    Assume that the line from end0 to end1 and axis go in the same direction,
    i.e. the dot product of axis and end1 - end0 is positive; this should be the
    case if end1 is a satellite above the horizon.
    """
    _, uax, length, _ = cyl_inf
    _, eproj, eperp, C = rxinfo(end0, cyl_inf)
    return end0 + _lci(eproj, eperp, C, end1 - end0, uax, length)

def cylinfo(base, length, radius):
    axis = earthnormal_xyz(base)
    return (base, axis, length, radius)

def rxinfo(rxloc, cyl_inf):
    base, uax, _, radius = cyl_inf
    end0 = rxloc - base
    eproj = end0 @ uax
    eperp = end0 - eproj*uax
    C = eperp @ eperp - radius**2
    return (rxloc, eproj, eperp, C)

def subrows(A, B):
    """Subtract each row in B from each row in A.

    If A is n*k and B is m*k, result is n*m*k.
    """
    n, k = A.shape
    m = B.shape[0]
    C = np.repeat(A, m, axis=0)
    C -= np.tile(B, (n, 1))
    C.shape = (n, m, k)
    return C

def heatmap(gloc, sxpos, cyl_inf, elthresh=None):
    np.seterr(invalid='ignore')
    if elthresh is None:
        elthresh = [5, 10, 20] # degrees above horizon (roughly)
    elthresh = np.array([cos((90 - el)*np.pi/180) for el in elthresh])
    base, uax, length, radius = cyl_inf
    heat = np.zeros(gloc.shape[:2] + (len(elthresh),), dtype=np.int32)
    sxpos = sxpos.reshape(-1, 3)
    end0 = gloc - base
    eproj = end0 @ uax
    eperp = end0 - np.multiply.outer(eproj, uax)
    C = np.sum(eperp**2, 2) - radius**2
    for i in range(len(gloc)):
        difs = subrows(sxpos, gloc[i])
        dproj = difs @ uax
        angdo = np.greater.outer(dproj / np.linalg.norm(difs, axis=2), elthresh)
        dperp = difs - np.multiply.outer(dproj, uax)
        A = 2 * np.sum(dperp**2, 2)
        B = 2 * np.sum(eperp[i] * dperp, -1)
        D = B**2 - 2*A*C[i]
        DD = np.sqrt(D) / A
        BB = -B / A
        SLNS = np.stack((BB - DD, BB + DD), -1)
        T0 = np.maximum(-eproj[i] / dproj, 0)
        T1 = (length - eproj[i]) / dproj
        GD = np.logical_and(SLNS[:, :, 0] < T1, SLNS[:, :, 1] > T0)
        heat[i] = np.sum(GD[:, :, np.newaxis] * angdo, 0)
    return heat

def drawheatmap(lons, lats, heat, hgts, vent, iskm=False, dolog=False):
    from matplotlib import colors
    plt.rcParams['axes.formatter.useoffset'] = False
    if iskm:
        aspect = 1
    else:
        aspect = lengthlat(vent[0]) / lengthlon(vent[0])
    fig = plt.figure(figsize=(8, 7.1875))
    ax = fig.add_axes([0.08, 0.07, 0.88, 0.88], aspect=aspect, anchor='SW')
                     #^left, bot, width, height
    ax.contour(lons, lats, hgts, colors='k') # this resizes the axes
    if dolog:
        surf = ax.pcolormesh(lons, lats, heat,
                             norm=colors.LogNorm(vmin=heat.min(), vmax=heat.max()))
    else:
        vmin = np.percentile(heat, 1)
        vmax = np.percentile(heat, 88)
        surf = ax.pcolormesh(lons, lats, heat, vmin=vmin, vmax=vmax)
    ax.autoscale(enable=False) # otherwise the scatter command messes up the axes
    ax.scatter(vent[1], vent[0], marker='^', s=100, c='r')
    box = ax.get_position()
    cax = fig.add_axes([0.9, box.y0, 0.04, box.height])
    fig.colorbar(surf, cax=cax)
    if iskm:
        ax.set_ylabel('North of peak (km)')
        ax.set_xlabel('East of peak (km)')
    else:
        ax.set_ylabel('Latitude (°)')
        ax.set_xlabel('Longitude (°)')
    ax.set_title('Satellites visible through cylinder')
    return ax

# Etna location info
peak_llh = (37.755278, 14.995277, 3340)
peak_xyz = (4879720.957403, 1307086.174806, 3886048.858672)

vent1_llh = (37.7466, 15.0020, 3020) # New South East Crater
vent1_xyz = (4879893.036667, 1307745.970325, 3885090.992257)

vent2_llh = (37.7507, 14.9922, 3200)
vent2_xyz = (4879984.887685, 1306875.998159, 3885561.186548)

def get_srtm(sdir='/homes/nima9589/SRTM'):
    files = ['N37E014.hgt', 'N37E015.hgt']
    files = [os.path.join(sdir, f) for f in files]
    # Loading SRTM data
    hgts1 = np.fromfile(files[0], dtype='>i2').reshape((3601, 3601))
    hgts2 = np.fromfile(files[1], dtype='>i2').reshape((3601, 3601))
    hgts = np.hstack((hgts1[:, :-1], hgts2))

    lons = np.arange(14, 16.0001, 1/3600)
    lats = np.arange(38, 36.9999, -1/3600)
    return hgts, lons, lats

def groundlocs(llh, dist):
    """Get ECEF ground locations for dist meters around llh coordinates."""
    hgts, lons, lats = get_srtm()
    londel = dist / lengthlon(llh[0])
    latdel = dist / lengthlat(llh[0])
    lonind = lons.searchsorted([llh[1] - londel, llh[1] + londel])
    latind = [findfirst.findfirstlt(lats, llh[0] + latdel),
              findfirst.findfirstlt(lats, llh[0] - latdel)]
    latrange = slice(*latind)
    lonrange = slice(*lonind)

    eht = hgts[latrange, lonrange]
    rlat = lats[latrange] * np.pi / 180
    rlon = lons[lonrange] * np.pi / 180
    slat = np.sin(rlat)
    enlat = WGS84.a / np.sqrt(1 - WGS84.e2 * slat**2)
    enlat.shape = (len(enlat), 1)
    slat.shape = (len(slat), 1)
    X = (enlat + eht) * np.outer(np.cos(rlat), np.cos(rlon))
    Y = (enlat + eht) * np.outer(np.cos(rlat), np.sin(rlon))
    Z = ((1 - WGS84.e2) * enlat + eht) * slat
    return np.stack((X, Y, Z), axis=-1), eht, lats[latrange], lons[lonrange]

def satpositions():
    pl = quickloadsp3('/homes/nima9589/igs19233.sp3')
    start = ((((1923 * 7) + 3) * 24 + 1) * 60 + 15) * 60
    cofns = [myinterp(start,
                      900,
                      [coeffs(range(start + 900*idx, start + 900*(idx+11), 900),
                              pl[idx:idx+11, prn, :]) for idx in range(len(pl) - 10)])
             for prn in range(32)]
    stop = start + 900*(len(pl) - 11)
    return np.array([[mvec(t) @ cofns[prn](t) for t in range(start, stop, 300)] for prn in range(32)])

"""
cyl_inf = cylinfo(vent1_xyz, 4000, 1414)
gloc, hgts, lats, lons = groundlocs(peak_llh, 6000)
sxpos = satpositions()
heat = heatmap(gloc, sxpos, cyl_inf)

nskm = (lats - peak_llh[0]) * lengthlat(peak_llh[0]) / 1000
ewkm = (lons - peak_llh[1]) * lengthlon(peak_llh[0]) / 1000
ventkm = (lengthlat(peak_llh[0], vent1_llh[0]) / 1000,
          (vent1_llh[1] - peak_llh[1]) * lengthlon(peak_llh[0]) / 1000)

drawheatmap(lons, lats, heat[:,:,0], hgts, vent1_llh) # or
drawheatmap(ewkm, nskm, heat[:,:,0], hgts, ventkm, True)
"""