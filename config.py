import sys
import numpy as np

LOGFILE = sys.stderr
DEBUG = True

BINEX_FILES="/data/VAPR/VB001/%Y/%j/data.bin"
"""Location of BINEX streams, in strftime format."""

SNR_MIN = 20
SNR_MAX = 56 # range to use for snr plots, to enable visual comparison between plots

PLOTDIR='/usr/local/adm/config/apache/htdocs/i/vapr/VB001'
PLOT_IVAL = np.timedelta64(5, 'm')
PLOT_SNR_HOURS = 4
PLOT_HK_HOURS = 24

SP3DIR = '/scratch/sp3'