#!/usr/bin/python
#
#Usage: ascii2snr.py <nmea_filename> <station_name>
###############################################################################
import sys
import datetime

fileName = sys.argv[1]
stationName = sys.argv[2]
file = open(fileName, mode='r')
file_h_v = []
doy_v = []

#Parse input file
gps_time = 0
gps_day = 0
gps_month = 0
gps_year = 0
i = 0
for line in file:
    line = line.strip()
    i = i + 1
    #print i
    if 'HK' in line or len(line) < 7:
        continue
    rmc = line.split(',')
    try:
        gps_time_utc = rmc[0]
        if len(gps_time_utc) < 9: # bad line
            continue
        gps_time = (float(gps_time_utc[0:2])*3600
                    + float(gps_time_utc[2:4])*60
                    + float(gps_time_utc[4:10]))
    except ValueError as e:
        print(e)
        continue
    #Get Date
    gps_date = rmc[1]
    if gps_date.isdigit() and len(gps_date) == 6:
        try:
            gps_day = int(gps_date[0:2])
            gps_month = int(gps_date[2:4])
            gps_year = 2000 + int(gps_date[4:6])
            doy = (datetime.date(gps_year, gps_month, gps_day).timetuple().tm_yday)
            doy_s = "%03d" % doy
            #Check if DOY exists
            if int(doy_s) in doy_v:
                file_h = file_h_v[doy_v.index(int(doy_s))]
            else:
                doy_v.append(int(doy_s))
                snr_fname = stationName + doy_s + '0.'+ str(gps_year)[2:4] + '.snr89'
                print(snr_fname)
                file_h = open(snr_fname, mode='w')
                file_h_v.append(file_h)
            for j in range(2, len(rmc)-1, 2):
                prn_s = rmc[j]
                snr = rmc[j+1]
                if not prn_s.isdigit():
                    continue
                try:
                    float(snr)
                except ValueError:
                    continue
                if 0 < int(prn_s) < 33:
                    file_h.write(prn_s + ' 0 0 ' + \
                            str(gps_time) + ' 0 0 ' + str(snr) + '\n')
        except ValueError as e:
            print("Value Error:", e, "on line", i)

file.close()
for file_h in file_h_v:
    file_h.close()
