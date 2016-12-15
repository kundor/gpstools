#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 14:12:59 2016

@author: nima9589
"""
from ascii_read import readsp3, vfilter

def gen_dps(lineiter):
    """Get datetime64, prn, SNR triples from each line in VAPR ASCII format.

    Lines should all be valid, so lineiter should either be a prefiltered
    file object, or run through vfilter(filename).
    """


def gen_dps89(filenames):
    """Get datetime64, prn, SNR triples from a list of files in snr89 format."""