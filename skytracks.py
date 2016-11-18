import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos
from plot import polarazel
from coords import xyz2llh

ETNA0 = (4878146, 1309017, 3886374) # ECEF in meters
REDBT = (-2802350, -1444879, 5527495)

sp3file = '/scratch/sp3/igs19162.sp3'
LOC = REDBT

#from satpos import _procheader
#from utility import fileread
#pl = np.full((96, 32, 3), np.NaN)

#with fileread(sp3file) as fid:
#    _procheader(fid)
#    epoch = 0
#    for line in fid:
#        if line[0] in ('E', 'V'):
#            continue
#        elif line[0] == '*':
#            epoch += 1

pl = np.loadtxt(sp3file, skiprows=22, usecols=(1, 2, 3), comments=['*', 'E', 'V'])
# This assumes all 32 PRNs are present, in order, for each epoch
# Also assuming interval is 900 seconds (15 minutes)

pl.shape = (-1, 32, 3) # reshape to 96x32x3 (-1 means 'fit the data')
pl *= 1000

lat, lon, ht = xyz2llh(LOC)
trans = np.array([[-sin(lon),           cos(lon),          0       ],
                  [-sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat)],
                  [ cos(lat)*cos(lon),  cos(lat)*sin(lon), sin(lat)]])
trans = trans.T

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
