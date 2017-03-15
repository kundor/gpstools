#!/bin/bash

rm rx050740.17.snr88, rx060740.17.snr88, HK_074.17.txt
../parser.py data.bin
diff rx050740.17.snr88 rx050740.17.snr88.GOOD
diff rx060740.17.snr88 rx060740.17.snr88.GOOD
diff HK_074.17.txt HK_074.17.txt.GOOD
