import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos
from plot import polarazel
from coords import xyz2llh

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
