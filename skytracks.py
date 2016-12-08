import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos
from math import sqrt
from plot import polarazel, plotcylinder
from coords import xyz2llh, earthnormal_xyz, lengthlat, lengthlon, WGS84
from satpos import myinterp, mvec, coeffs

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
    d = b*b - 4*a*c
    if d < 0:
        return []
    if d > 0:
        return [(-b - sqrt(d))/2/a, (-b + sqrt(d))/2/a]
    return [-b/2/a]

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

def _lci(end0, dif, uax, length, radius):
    # Parameterize the line by L(t) = end0 + t*dif
    eproj = end0 @ uax
    dproj = dif @ uax
    eperp = end0 - eproj*uax
    dperp = dif - dproj*uax
    A = dperp @ dperp
    B = 2 * eperp @ dperp
    C = eperp @ eperp - radius**2
    slns = quadformula(A, B, C)
    t0 = -eproj / dproj # where line meets bottom of cylinder
    t1 = (length - eproj) / dproj # where line meets top of cylinder
    if not _clip(slns, t0, t1):
        return np.zeros((0,3)) # empty array
    return end0 + np.outer(slns, dif)

def linecylinderintersect(end0, end1, base, axis, radius):
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
    length = np.linalg.norm(axis)
    uax = axis / length
    return base + _lci(end0 - base, end1 - end0, uax, length, radius)
#    blens = np.dot(ipts, uax) # make sure 0 < blen < length

def isectp(rxloc, sxloc, ventloc, radius, height):
    """Does the line from rxloc to sxloc intersect the cylinder above ventloc with given radius and height?"""
    if linecylinderintersect(rxloc, sxloc, ventloc, earthnormal_xyz(ventloc)*height, radius).size:
        return True
    return False

def cylinfo(base, length, radius):
    axis = earthnormal_xyz(base)
    return (base, axis, length, radius)

def misectp(rxloc, sxloc, cyl_inf):
    base, uax, length, radius = cyl_inf
    if _lci(rxloc - base, sxloc - rxloc, uax, length, radius).size:
        return True
    return False

# Loading SRTM data
hgts1 = np.fromfile('N37E014.hgt', dtype='>i2').reshape((3601, 3601))
hgts2 = np.fromfile('N37E015.hgt', dtype='>i2').reshape((3601, 3601))
hgts = np.hstack((hgts1[:,:-1], hgts2))

LONS = np.arange(14, 16.0001, 1/3600)
LATS = np.arange(38, 36.9999, -1/3600)

peak_llh = (37.755278, 14.995277, 3340)
peak_xyz = (4879720.957403, 1307086.174806, 3886048.858672)

vent1_llh = (37.7466, 15.0020, 3020) # New South East Crater
vent1_xyz = (4879893.036667, 1307745.970325, 3885090.992257)

vent2_llh = (37.7507, 14.9922, 3200)
vent2_xyz = (4879984.887685, 1306875.998159, 3885561.186548)

# We look at lats 37.71 to 37.8, lons 14.94 to 15.05 (about 5 km around the peak)
latrange = slice(720, 1045) # .2*3600 -- .29 * 3600
lonrange = slice(3384, 3781) # .94*3600 -- 1.05*3600
ax = plt.subplot(1, 1, 1)
ax.contour(LONS[lonrange], LATS[latrange], hgts[latrange, lonrange], colors='k')
ax.set_aspect(lengthlat(peak_llh[0]) / lengthlon(peak_llh[0]))

eht = hgts[latrange, lonrange]
rlat = LATS[latrange] * np.pi / 180
rlon = LONS[lonrange] * np.pi / 180
slat = np.sin(rlat)
enlat = WGS84.a / np.sqrt(1 - WGS84.e2 * slat * slat)
enlat.shape = (len(enlat), 1)
slat.shape = (len(slat), 1)
X = (enlat + eht) * np.outer(np.cos(rlat), np.cos(rlon))
Y = (enlat + eht) * np.outer(np.cos(rlat), np.sin(rlon))
Z = ((1 - WGS84.e2) * enlat + eht) * slat
gloc = np.stack((X, Y, Z), axis=-1)

pl = quickloadsp3('igs19233.sp3')
start = ((((1923 * 7) + 3) * 24 + 1) * 60 + 15) * 60
cofns = [myinterp(start,
                  900,
                  [coeffs(range(start + 900*idx, start + 900*(idx+11), 900),
                          pl[idx:idx+11, prn, :]) for idx in range(len(pl) - 10)])
         for prn in range(32)]
stop = start + 900*(len(pl) - 11)
sxpos = np.array([[mvec(t) @ cofns[prn](t) for t in range(start, stop, 60)] for prn in range(32)])

heat = np.zeros_like(eht)

cyl_inf = cylinfo(vent1_xyz, 6000, 2000)

for ri in np.ndindex(heat.shape):
    for si in np.ndindex(sxpos.shape[:2]):
        heat[ri] += misectp(gloc[ri], sxpos[si], cyl_inf)
