#script to import and process waveforms containing both P-wave and S-wave
#picks from analysts. Produces a set of files which could be used to train a
#machine-learning algorithm to recognise P and S waves from seismic traces.


#this one does teleseismic (>10 degrees) S-wave picks from the ISC catalogue

from mat4py import *
import numpy as np

from obspy.core import UTCDateTime

from getwave import getWave

import csv

#initialise IRIS client to query if the desired channel is not in our database.
#this client needs to be initialised here because the initialiser spawns multiple threads.
#this is forbidden on import in Python 2 so it cannot be initialised in the getwave module
from obspy.clients.fdsn.client import Client
irisclient=Client("IRIS")

#load ISC pick catalogue, stored in 4 separate .mat files

Sctr=0
wfctr=0
with open('/g/data/ha3/Passive/Events/BabakHejrani/ISC.csv') as ISCpicks:
    pickrdr=csv.reader(ISCpicks,delimiter=',')
    for pick in pickrdr:
        st=pick[0].strip()
        if st=='#ISC': #ignore event metadata
            continue
        ph=pick[7].strip()
        dist=float(pick[-1].strip())
        if ph=='S' and dist > 10:
            Sctr+=1
            ch=pick[1].strip()
            if '?' in ch:
                continue
            at=[float(num.strip()) for num in pick[8:-1]]
            at.append((at[-1]-int(at[-1]))*1e6)
            at=[int(num) for num in at]
            at=UTCDateTime(*at)
            starttime=at-20
            endtime=at+40
            succ=getWave(irisclient,st,ch,starttime,endtime)
            if succ:
                wfctr+=1
                        
            

    #for ind,phases in enumerate(ISCpicks['ph']):
    #    for sub,ph in enumerate(phases):
    #        if ph=='S' and ISCpicks['dist'][ind][sub][0]>10:
#                Stime=[int(num) for num in ISCpicks['at'][ind][sub]]
#                Stime=Stime+[int((ISCpicks['at'][ind][sub][-1]-Stime[-1])*1e6)]
#                Stime=UTCDateTime(*Stime)
#                #save a minute around the S pick
#                starttime=Stime-20
#                endtime=Stime+40
#                
#                #try to save a waveform                
#                succ=getWave(ISCpicks['st'][ind][sub],ISCpicks['ch'][ind][sub],starttime,endtime)
#                if succ:
#                    Sctr+=1
#print Sctr
