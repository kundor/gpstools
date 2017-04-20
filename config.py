"""Configuration variables.

These are always used as config.NAME (never with `from config import NAME`)
so that changing the values interactively will be respected.

If a file called site_config.py is present, any values defined there
will override the values in this file. site-config.py is ignored by version control,
so you may put settings there specific for your machine without messing up
the defaults. It is also where private API keys are defined, to avoid exposure
on github.
"""

import sys
import numpy as np

SITE = 'VB001'
"""Name of the VAPR site to process."""

LOGFILE = sys.stderr
"""Where informative messages are sent. (Set to a file for headless use.)"""

DEBUG = False
"""Whether to output debugging messages (generally not of interest.)"""

BINEX_FILES = '/data/VAPR/{site}/%Y/%j/data.bin'
"""Location of BINEX streams, in strftime format; {site} is replaced by config.SITE."""

SNR_RANGE = (30, 56)
"""Range of values to use for snr plots, to enable visual comparison between plots."""

VOLT_RANGE = (3.3, 4.5)
"""Range to use for voltage plots. (This will be extended to accommodate extreme values.)"""

TEMP_RANGE = (-5, 35)
"""Range to use for temperature plots. (This will be extended to accommodate extreme values.)"""

MINELEV = 10
"""Minimum elevation to use for plotting SNR or no. tracked satellites."""

PLOTDIR = '/usr/local/adm/config/apache/htdocs/i/vapr/{site}'
"""Where to save plot images for the web site; {site} is replaced by config.SITE."""

DAILYDIR = PLOTDIR + '/%Y/%j'
"""Location of archive web page created every UTC midnight, in strftime format; {site} is replaced by config.SITE."""

PLOT_IVAL = np.timedelta64(5, 'm')
"""How often to create plots."""

PLOT_SNR_HOURS = 4
"""How many hours of data to show in SNR plots."""

PLOT_HK_HOURS = 24
"""How many hours of data to show in housekeeping plots (temp/volt, number of satellites, mean snr.)"""

SP3DIR = '/scratch/sp3'
"""Where sp3 orbit position files are stored."""

ERR_EMAIL = 'vapr@kundor.org'
"""Where to email exception reports."""

FROM_EMAIL = 'vapr@localhost'
"""Address to use as the from address in emailed reports."""

SMTP_HOST = 'localhost'
"""SMTP server for sending emails."""

GKEY = None

PUSHER = None

try:
    from site_config import *
except ImportError:
    pass
