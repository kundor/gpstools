"""
Functions to convert between coordinate systems, including
Earth-centered, Earth-fixed Cartesian coordinates (in meters);
latitude, longitude, and height above the ellipsoid (in radians and meters);
East-North-Up from a base point (in meters);
and azimuth and elevation from a base point (in radians).
"""

from math import atan2, cos, sin, sqrt, pi
from urllib.request import urlopen
from urllib.parse import urlencode
from urllib.error import HTTPError
import json
import re
import scipy.integrate
import numpy as np
import config
from utility import info

class WGS84:
    """Parameters defining the WGS84 ellipsoid."""
    a = 6378137.
    f = 1./298.257223563
    e2 = 2*f - f*f
    b = a * (1 - f)

    @classmethod
    def ellnormal(cls, lat):
        """Distance from the surface of the ellipsoid at latitude lat to the
        z-axis, along the ellipsoid normal.
        """
        sl = sin(lat)
        return cls.a / sqrt(1 - cls.e2 * sl * sl)
# equivalent: a^2 / sqrt(a^2 cos^2(lat) + b^2 sin^2(lat))

def lengthlat(lat1, lat2=None):
    """Distance in meters between two given latitudes in degrees.

    If only one argument is given, return length of one degree of latitude
    centered at lat1.
    """
    if lat2 is None:
        lat2 = lat1 + 0.5
        lat1 -= 0.5
    lat1 *= pi / 180
    lat2 *= pi / 180
    def integrand(x):
        sx = sin(x)
        return pow(1 - WGS84.e2 * sx * sx, -3/2)
    reserr = scipy.integrate.quad(integrand, lat1, lat2)
    return WGS84.a * (1 - WGS84.e2) * reserr[0]

def lengthlon(lat):
    """Length (m) of one degree of longitude at given latitude (degrees)."""
    lat *= pi / 180
    return pi * WGS84.ellnormal(lat) * cos(lat) / 180

def earthradius(lat):
    """Distance from the earth's center to the WGS84 ellipsoid at the given latitude (radians)."""
    cl2 = cos(lat)**2
    sl2 = sin(lat)**2
    return sqrt((WGS84.a**4 * cl2 + WGS84.b**4 * sl2)
                / (WGS84.a**2 * cl2 + WGS84.b**2 * sl2))

def jsonfetch(url, **kwargs):
    """Make a GET request to url with parameters from the keyword arguments; parse result as JSON."""
    with urlopen(url + '?' + urlencode(kwargs)) as resp:
        enc = resp.info().get_content_charset('utf-8')
        return json.loads(resp.read().decode(enc))

def google_ht(lat, lon, bestguess=0):
    """Fetch height above the geoid (m) from Google Maps at the given coordinates (decimal degrees.)

    Falls back to third parameter bestguess (default 0) if the HTTP request fails.
    """
    gurl = 'https://maps.googleapis.com/maps/api/elevation/json'
    params = {'locations': '{},{}'.format(lat, lon)}
    if config.GKEY:
        params['key'] = config.GKEY
    try:
        resp = jsonfetch(gurl, **params)
        return resp['results'][0]['elevation']
    except (HTTPError, KeyError, IndexError) as e:
        info(e)
        return bestguess

def orthometric_ht(lat, lon):
    """Fetch ellipsoid height - geoid height from online UNAVCO calculator.

    Given lat, lon in decimal degrees, return orthometric height (m).
    Add the result to the geoid height to get ellipsoid height.
    """
    url = 'http://jules.unavco.org/Geoid/Geoid'
    dat = urlencode(dict(lat=lat, lon=lon, gpsheight=0)).encode('ascii')
    try:
        with urlopen(url, dat) as resp:
            enc = resp.info().get_content_charset('utf-8')
            text = resp.read().decode(enc)
    except HTTPError as e:
        info(e)
        return None
    mch = re.search(r'Orthometric height [^0-9]* EGM96 [^0-9]* (-?[0-9]+\.[0-9]+)', text)
    if not mch:
        return None
    return float(mch.group(1))

def get_ellipsoid_ht(lat, lon):
    """Use Google Maps API and UNAVCO calculator to determine ellipsoid height (m).

    lat, lon in decimal degrees.
    """
    orth = orthometric_ht(lat, lon) or 0
    # if we fail to get the difference, we just use the geoid height
    return round(google_ht(lat, lon) + orth, 3)
    # Unavco only returns three digits, keeping the other digits is silly

def earthnormal_xyz(x, y=None, z=None):
    """Unit normal vector to the WGS84 ellipsoid at the given ECEF coordinates."""
    if y is None and z is None:
        x, y, z = x
    wa = (WGS84.a / 2500)**2
    wb = (WGS84.a * (1 - WGS84.f) / 2500)**2
    nrml = np.array((x / wa, y / wa, z / wb))
    return nrml / sqrt(nrml @ nrml)

def earthnormal_llh(lat, lon=None, ht=None):
    """Unit normal vector to the WGS84 ellipsoid at the given coordinates

    Input is geodetic latitude and longitude in radians, and height above
    the ellipsoid in meters.
    """
    if lon is None and ht is None:
        lat, lon, ht = lat
    xyz1 = np.array(llh2xyz(lat, lon, ht))
    xyz2 = np.array(llh2xyz(lat, lon, ht + 1))
    return xyz2 - xyz1

def deg2rad(lat, lon=None, ht=None):
    """Convert latitude and longitude from degrees to radians.

    Height is returned unchanged, if given.
    """
    if lon is not None and ht is None:
        return lat * pi / 180, lon * pi / 180
    if lon is None and ht is None:
        if len(lat) == 2:
            return lat[0] * pi / 180, lat[1] * pi / 180
        if len(lat) == 3:
            return lat[0] * pi / 180, lat[1] * pi / 180, lat[2]
    return lat * pi / 180, lon * pi / 180, ht

def rad2deg(lat, lon=None, ht=None):
    """Convert latitude and longitude from radians to degrees.

    Height is returned unchanged, if given.
    """
    if lon is not None and ht is None:
        return lat * 180 / pi, lon * 180 / pi
    if lon is None and ht is None:
        if len(lat) == 2:
            return lat[0] * 180 / pi, lat[1] * 180 / pi
        if len(lat) == 3:
            return lat[0] * 180 / pi, lat[1] * 180 / pi, lat[2]
    return lat * 180 / pi, lon * 180 / pi, ht

def xyz2lat(x, y, z, tol=1e-10):
    p2 = x*x + y*y
    oe2z2 = (1 - WGS84.e2) * z*z
    k = 1/(1 - WGS84.e2)
    err = 1
    while err > tol:
        c = (p2 + oe2z2*k*k)**(3/2)/WGS84.a/WGS84.e2
        k1 = 1 + (p2 + oe2z2*k*k*k)/(c - p2)
        err = abs(k1 - k)
        k = k1
    return atan2(z*k, sqrt(p2))

def xyz2llh(x, y=None, z=None, tol=1e-10):
    """Given ECEF Cartesian coordinates in meters, return
    latitude & longitude in radians, and height above ellipsoid in meters.
    """
    if y is None and z is None:
        x, y, z = x
    lon = atan2(y, x)
    lat = xyz2lat(x, y, z, tol)
    p = sqrt(x*x + y*y)
    hgt = p / cos(lat) - WGS84.ellnormal(lat)
    return lat, lon, hgt

def llh2xyz(lat, lon=None, ht=None):
    """Given (geodetic) latitude & longitude in radians,
    and height above ellipsoid in meters, return ECEF Cartesian
    coordinates in meters.
    """
    if lon is None and ht is None:
        lat, lon, ht = lat
    enlat = WGS84.ellnormal(lat)
    x = (enlat + ht) * cos(lat) * cos(lon)
    y = (enlat + ht) * cos(lat) * sin(lon)
    z = ((1 - WGS84.e2) * enlat + ht) * sin(lat)
    return x, y, z

def _enutrans(lat, lon, vec):
    """Matrix to help turn ECEF coordinates to ENU coordinates at given reference point."""
#    [[-sin(lon),           cos(lon),          0       ],
#     [-sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat)],
#     [ cos(lat)*cos(lon),  cos(lat)*sin(lon), sin(lat)]])
    x = -sin(lon)*vec[0] + cos(lon)*vec[1]
    y = -sin(lat)*cos(lon)*vec[0] - sin(lat)*sin(lon)*vec[1] + cos(lat)*vec[2]
    z =  cos(lat)*cos(lon)*vec[0] + cos(lat)*sin(lon)*vec[1] + sin(lat)*vec[2]
    return x, y, z


def xyz2enu(base, pt):
    """Transform ECEF coordinates to local East-North-Up coordinates.

    Given ECEF Cartesian coordinates (m) for a base point and _pt_
    (as numpy arrays), return an array of three coordinates (m).
    """
    lat, lon, ht = xyz2llh(base)
    return _enutrans(lat, lon, pt - base)

def enu2azel(e, n=None, u=None):
    """Given East-North-Up coordinates, return azimuth and elevation.

    Azimuth is clockwise from north. Both are in radians.
    """
    if n is None and u is None:
        e, n, u = e
    az = atan2(e, n)
    if az < 0:
        az += 2*pi
    p = sqrt(e*e + n*n)
    el = atan2(u, p)
    return az, el






