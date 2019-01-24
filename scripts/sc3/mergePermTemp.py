#!/usr/bin/env python

import obspy
import sys

stPerm = obspy.read(sys.argv[1])
stTemp = obspy.read(sys.argv[2])

for tr in stTemp:
	stPerm += tr

stPerm.write(sys.argv[3], format='MSEED')
