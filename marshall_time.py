#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Scrape Marshall Field temp from http://www.ral.ucar.edu/projects/winter/sites/marshall/

Created on Tue May  9 10:50:06 2017

@author: nima9589
"""
from urllib.request import urlopen
import re
import datetime
import time

def get_timetemp():
    """Fetch temperature from Marshall Field site."""
    url = 'http://www.ral.ucar.edu/projects/winter/sites/marshall/'
    with urlopen(url) as resp:
        enc = resp.info().get_content_charset('utf-8')
        text = resp.read().decode(enc)
    mch = re.search(r'<td>Time<br /><span class="current_weather_val">([^<]*)</span> UTC</td>.*'
            r'<td>Temperature<br /><span class="current_weather_val">(-?[0-9]+\.[0-9]+)</span> C</td>', text, re.DOTALL)
    if not mch:
        raise RuntimeError('Could not find values')
#    dtim = datetime.datetime.strptime(mch.group(1), '%Y-%m-%d %H:%M:%S')
    return mch.groups

def scrapelog(file='marshall_temp.txt'):
    with open(file, 'wt') as fid:
        while True:
            vals = get_timetemp()
            print(' '.join(vals))
            fid.write(' '.join(vals))
            time.sleep(600)
