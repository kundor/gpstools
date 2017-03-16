import sys
import numpy as np

LOGILE = sys.stderr
DEBUG = True

SNR_MIN = 20
SNR_MAX = 56 # range to use for snr plots, to enable visual comparison between plots

PLOTDIR='/usr/local/adm/config/apache/htdocs/i/vapr/VB001'
PLOT_IVAL = np.timedelta64(5, 'm')
PLOT_HOURS = 4

SP3DIR = '/scratch/sp3'