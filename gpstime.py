# Created by Nick Matteo <kundor@kundor.org> June 9, 2009
"""Utility time functions useful for GPS work.

LeapSeconds is a class providing a dictionary of UTC datetimes and
corresponding leapsecond offsets from TAI.
LeapSeconds.update() is a class method to update the leap seconds information.
`leapseconds' is an instantiation of the leap second dictionary.

NB: Standard Python datetime objects are only precise to 1 microsecond.
"""
# These classes do NOT account for:
#  - Difference in GPS time vs. UTC(USNO) (Currently sync'd once per day.)
#  - Difference in UTC(USNO) vs. UTC (sync'd about once a month.)
#  - Difference in UTC vs. GLONASS
#  - Difference in TAI vs. Galileo
# I don't know about UTC vs. GLONASS, but the other errors remain on the order
# of a few nanoseconds.  The sum of all such errors should remain well below a
# microsecond (the limit of Python's datetime precision.)

import re
import os
from os import path
from urllib.request import urlopen
from urllib.error import URLError
from datetime import datetime, timedelta, timezone
import time
from warnings import warn
from collections.abc import Sequence
from numbers import Number
import numpy as np
from utility import info


URL1 = 'http://maia.usno.navy.mil/ser7/tai-utc.dat'
URL2 = 'http://hpiers.obspm.fr/iers/bul/bulc/UTC-TAI.history'

GPSepoch = np.datetime64('1980-01-06', 'us')
GPSepochdt = datetime(1980, 1, 6, 0, 0, 0, 0, timezone.utc)

def gpsweekgps(ndt):
    """GPS week number of a numpy datetime64 in GPS time."""
    return int((ndt - GPSepoch) / np.timedelta64(1, 'W'))

def gpsweekutc(ndt):
    """Compute GPS week number of a numpy datetime64 in UTC time."""
    return gpsweekgps(ndt + leapsecs(ndt))

def gpstotsecgps(ndt):
    """GPS seconds since epoch for a numpy datetime64 in GPS time."""
    return (ndt - GPSepoch) / np.timedelta64(1, 's')

def gpssowgps(ndt):
    """GPS second of week of a numpy datetime64 in GPS time."""
    return gpstotsecgps(ndt) % (60*60*24*7)

def gpssowutc(ndt):
    """GPS second of week of a numpy datetime64 in UTC time."""
    return gpssowgps(ndt + leapsecs(ndt))


def gpsweek(dt):
    """Given aware datetime, return GPS week number (adding leapseconds)."""
    dt += timedelta(seconds=gpsleapsecsutc(dt))
    return int((dt - GPSepochdt) / timedelta(weeks=1))

def _dow(dt):
    return dt.isoweekday() % 7

def gpsdow(dt):
    """Given aware datetime, return GPS day of week (adding leapseconds)."""
    return _dow(dt.astimezone(timezone.utc) + timedelta(seconds=gpsleapsecsutc(dt)))

def _sod(dt):
    return (dt.hour*60 + dt.minute)*60 + dt.second + dt.microsecond/1000000

def gpssod(dt):
    """Given aware datetime, return GPS second of day (adding leapseconds)."""
    return _sod(dt.astimezone(timezone.utc) + timedelta(seconds=gpsleapsecsutc(dt)))

def _sow(dt):
    return _dow(dt)*86400 + _sod(dt)

def gpssow(dt):
    """Given aware datetime, return GPS second of week (adding leapseconds)."""
    return _sow(dt.astimezone(timezone.utc) + timedelta(seconds=gpsleapsecsutc(dt)))

def isnaive(dt):
    """Return true if input is a naive datetime."""
    return isinstance(dt, datetime) and (dt.tzinfo is None or dt.tzinfo.utcoffset(dt) is None)

def dhours(hrs):
    """Convenience function: returns timedelta of given # of hours."""
    return timedelta(hours=hrs)

def getutctime(dt=None, tz=timezone.utc):
    """Convert time to a UTC-aware datetime object, or naive if tz=None.

    Accepts datetime, struct_time, tuple (Y,M,D,[H,M,S,μS]), tuple (GPS week,
    GPS second of week), numpy datetime64, or POSIX timestamp.
    Input is assumed to already be in desired timezone unless the timezone is included (aware
    datetime objects and struct_times).
    """
    if dt is None and tz:
        return datetime.now(tz)
    elif dt is None:
        return datetime.utcnow()
    elif isinstance(dt, np.datetime64):
        return dt.tolist().replace(tzinfo=tz)
    elif isinstance(dt, time.struct_time):
        nt = datetime(*dt[:6])
        if dt.tm_gmtoff is not None:
            nt -= timedelta(seconds=dt.tm_gmtoff)
        if tz:
            return nt.replace(tzinfo=timezone.utc).astimezone(tz)
        return nt
    elif isinstance(dt, Sequence) and len(dt) == 2: # list or tuple: GPS Week, GPS second of week
        gt = GPSepochdt + timedelta(weeks=dt[0], seconds=dt[1])
        if tz:
            return gt.astimezone(tz)
        return gt.replace(tzinfo=None)
    elif isinstance(dt, Sequence): # list or tuple: Year, Month, Day, H, M, S, μS
        return datetime(*dt, tzinfo=tz)
    elif isinstance(dt, Number):
        return datetime.fromtimestamp(dt, tz)
    elif isnaive(dt):
        return dt.replace(tzinfo=tz)
    elif isinstance(dt, datetime) and tz:
        return dt.astimezone(tz)
    elif isinstance(dt, datetime):
        return dt.astimezone(timezone.utc).replace(tzinfo=None)
    raise ValueError("Don't know how to interpret this as time")

class LeapSeconds(dict):
    """A dictionary of datetimes : leap second adjustment.

    Uses data file leapseco.dat, in same directory as the code.
    TAI differs from UTC by the adjustment at the latest datetime before the given epoch.
    NB: For dates before 1972, there is secular variation in the adjustment,
    which is NOT accounted for.
    """
    infofile = path.join(path.dirname(path.abspath(__file__)), 'leapseco.dat')

    def __init__(self):
        """Load and parse leap seconds data file."""
        dict.__init__(self)
        try:
            lfile = open(self.infofile)
        except IOError:
            warn('Leap seconds data file not found.  Attempting download.')
            if self.update():
                lfile = open(self.infofile)
            else:
                raise RuntimeError('Leap seconds data file not available.')
        for line in lfile:
            match = re.match('^([0-9:/-]+) : ([0-9.-]+)$', line)
            if match:
                dt = datetime(*(time.strptime(match.group(1), '%Y/%m/%d-%H:%M:%S')[0:6]))
                self[dt] = float(match.group(2))
        lfile.close()

    @classmethod
    def timetoupdate(cls):
        """Attempt to verify whether January 1 or July 1 has passed since last update.

        Otherwise, don't bother updating (new leap seconds only occur on those dates.)
        """
        now = datetime.utcnow()
        if not os.access(cls.infofile, os.R_OK):
            return True  # If file isn't there, try update
        try:
            fid = open(cls.infofile)
        except IOError:
            return True  # ditto
        try:
            updtime = datetime(*(time.strptime(fid.readline(),
                                               'Updated: %Y/%m/%d\n')[0:6]))
        except ValueError:
            warn('Leap second data file in invalid format.')
            return True
        fid.close()
        if updtime > now:
            warn('Leap second data file is from the future.')
            return False
        elif updtime.month <= 6:
            target = datetime(updtime.year, 7, 1)
        else:
            target = datetime(updtime.year + 1, 1, 1)
        if now <= target:
            return False
        return True

    @classmethod
    def update(cls):
        """Download and parse new leap second information from reliable web sources."""
        if not cls.timetoupdate():
            info('No potential leap second has occurred since last update.')
            return False
        if not os.access(path.dirname(cls.infofile), os.W_OK):
            raise IOError('Leap second data file cannot be written.')
        try:
            upd = urlopen(URL1)
            form = 1
        except URLError:
            upd = urlopen(URL2)
            form = 2
        mons = ['', b'JAN', b'FEB', b'MAR', b'APR', b'MAY', b'JUN', b'JUL',
                b'AUG', b'SEP', b'OCT', b'NOV', b'DEC']
        newfile = path.join(path.dirname(__file__), 'leapseco.new')
        lfile = open(newfile, 'w')
        lfile.write('Updated: ' + datetime.utcnow().strftime('%Y/%m/%d\n'))
        year = 1961
        for line in upd:
            if form == 1:
                year = int(line[0:5])
                month = mons.index(line[6:9].upper())
                day = int(line[10:13])
                adjust = float(line[36:48])
            elif form == 2:
                if len(line) < 36 or re.match(b' ?-* ?$| RELATIONSHIP| Limits', line):
                    continue
                if line[0:6].strip() != b'':
                    year = int(line[0:6])
                month = mons.index(line[7:10].upper())
                day = int(line[12:15].rstrip(b'. '))
                adjust = float(line[31:47].replace(b' ', b'').rstrip(b's\t\n'))
            lfile.write(datetime(year, month, day).strftime('%Y/%m/%d-%H:%M:%S')
                        + ' : ' + str(adjust) + '\n')
        upd.close()
        lfile.close()
        if path.exists(cls.infofile):
            os.remove(cls.infofile)
        os.rename(newfile, cls.infofile)
        return True

leapseconds = LeapSeconds()

def leapsecs(dt, cmp):
    """# of leapseconds at datetime dt.

    Whether dt exceeds the UTC time in the leapseconds dict is determined by
    the provided function cmp.
    """
    dt = dt.replace(tzinfo=None)
    if dt.year < 1958:
        raise ValueError('TAI vs UTC is unclear before 1958; unsupported.')
    try:
        return leapseconds[max([l for l in leapseconds if cmp(l, dt)])]
    except ValueError:
        return 0 # before 1961-Jan-01, TAI = UTC

def leapsecsutc(utc):
    """# of TAI-UTC leapseconds at UTC datetime."""
    return leapsecs(utc, lambda l, dt : l <= dt)

def gpsleapsecsutc(utc):
    """# of GPS-UTC leapseconds at UTC datetime."""
    return leapsecs(utc, lambda l, dt : l <= dt) - 19

def leapsecstai(tai):
    """# of TAI-UTC leapseconds at TAI datetime."""
    return leapsecs(tai, lambda l, dt : leapseconds[l] <= (dt - l).total_seconds())

def leapsecsnp(ndt):
    """GPS leap seconds at a numpy datetime64 in UTC time, as timedelta64.

    Doesn't work for dates before 1972.
    """
    ls = leapseconds[max(l for l in leapseconds if l <= ndt)] - 19
    return np.timedelta64(int(ls), 's')
