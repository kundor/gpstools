import sys
import numpy as np

LOGFILE = sys.stderr
DEBUG = False

BINEX_FILES="/data/VAPR/VB001/%Y/%j/data.bin"
"""Location of BINEX streams, in strftime format."""

SNR_RANGE = (20, 56) # range to use for snr plots, to enable visual comparison between plots
VOLT_RANGE = (3.3, 4.5) # range to use for voltage plots
TEMP_RANGE = (-5, 35) # ... and temperature plots

PLOTDIR='/usr/local/adm/config/apache/htdocs/i/vapr/VB001'
PLOT_IVAL = np.timedelta64(5, 'm')
PLOT_SNR_HOURS = 4
PLOT_HK_HOURS = 24

SP3DIR = '/scratch/sp3'