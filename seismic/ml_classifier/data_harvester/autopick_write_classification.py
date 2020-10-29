#write a new ensemble file for the automated picks, including a column indicating
#if they were accepted as genuine by the CNN.
import sys
sys.path.append('/g/data1a/ha3/rlt118/hiperseis/')

import numpy as np
import re

loadDir='/g/data/ha3/rlt118/neural-datasets/autopicks/'
with open(loadDir+'ensemble.s.txt','r') as f:
    picks=f.readlines()
    pickCtr=len(picks)-1

with open(loadDir+'ensemble.s.verified.txt', 'r') as f:
    vpicks = f.readlines()
    vpicks = vpicks[1:]



vpdict = {}
for pick in vpicks:
    plist = pick.split()
    snr = float(plist[17])
    evID = plist[0]
    net = plist[6]
    st = plist[7]
    ch = plist[8]
    vpdict[evID+net+st+ch] = 1


outfile = loadDir + "ensemble.s.classified.txt"
with open(outfile, "w") as f:
    #header
    f.write("#eventID originTimestamp mag originLon originLat originDepthKm net sta cha pickTimestamp phase stationLon stationLat az baz distance ttResidual snr qualityMeasureCWT domFreq qualityMeasureSlope bandIndex nSigma CNNresult\n")
    for pick in picks:
        plist = pick.split()
        cwt = float(plist[18])
        snr = float(plist[17])

        evID = plist[0]
        net=plist[6]
        st=plist[7]
        ch=plist[8]
        if not (evID+net+st+ch in vpdict):
            #pick was not verified by CNN
            f.write(pick.rstrip()+" 0\n")
        else:
            #pick was verified by CNN
            f.write(pick.rstrip()+" 1\n")
