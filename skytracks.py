import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos
from plot import polarazel
from coords import xyz2llh, earthnormal_xyz, lengthlat, lengthlon, WGS84

ETNA0 = (4878146, 1309017, 3886374) # ECEF in meters
REDBT = (-2802350, -1444879, 5527495)

def skytracks(sp3file, LOC):
    pl = np.loadtxt(sp3file, skiprows=22, usecols=(1, 2, 3), comments=['*', 'E', 'V'])
    # This assumes all 32 PRNs are present, in order, for each epoch
    # Also assuming interval is 900 seconds (15 minutes)

    pl.shape = (-1, 32, 3) # reshape to 96x32x3 (-1 means 'fit the data')
    pl *= 1000

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

skytracks('/scratch/sp3/igs19162.sp3', REDBT)

def linecylinderintersect(end0, end1, axis, radius):
    """Intersections of a line segment from end0 to end1 with a cylinder centered on axis.

    Axis should be a unit vector.
    """
    theta = np.arctan(axis[1]/axis[0])
    phi = np.arcsin(axis[2])
    sithe2 = np.sin(theta)**2
    cothe2 = np.cos(theta)**2
    siphi2 = axis[2]**2
    cophi2 = np.cos(phi)**2
    D = sithe2 + cothe2 * siphi2
    E = cothe2 + sithe2 * cophi2
    # Equation of the cylinder is D x^2 + E y^2 + cophi2 z^2 = radius^2
    dif = end1 - end0
    A = np.outer((D, E, cophi2), dif*dif)
    B = 2*


    # Parameterize the line by L(t) = end0 + t (end1 - end0)
    # this intersects the cylinder where D t^2 + E t



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

vent1_axis = earthnormal_xyz(vent1_xyz)
vent2_axis = earthnormal_xyz(vent2_xyz)
