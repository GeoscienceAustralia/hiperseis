#script to import and process waveforms containing both P-wave and S-wave
#picks from analysts. Produces a set of files which could be used to train a
#machine-learning algorithm to recognise P and S waves from seismic traces.

#requires python >= 2.7

#this one does teleseismic (>10 degrees) S-wave picks from the ISC catalogue
import sys


from mat4py import *
import numpy as np

from obspy.core import UTCDateTime
from obspy.core.stream import Stream


from getwave import getWave

import csv

#initialise IRIS client to query if the desired channel is not in our database.
#this client needs to be initialised here because the initialiser spawns multiple threads.
#this is forbidden on import in Python 2 so it cannot be initialised in the getwave module
from obspy.clients.fdsn.client import Client
from obspy.geodetics.base import gps2dist_azimuth as distaz
irisclient=Client("IRIS")
#load ISC pick catalogue CSV

saveDir="/g/data/ha3/rlt118/neural-datasets/categoriser-teleseismic/smallset-2/"

Sctr=0
wfctr=0
with open('/g/data/ha3/Passive/Events/MergeCatalogs/ISC.csv') as ISCpicks, open(saveDir+'picklog.csv','w') as picklog:
    pickrdr=csv.reader(ISCpicks,delimiter=',')
    pickwtr=csv.writer(picklog,delimiter=',')
    event=""
    for pick in pickrdr:
        st=pick[0].strip()
        if st=='#ISC': #store event metadata then don't try to process this as a pick
            event=pick
            evlat=float(event[6].strip())
            evlong=float(event[7].strip())
            continue
        ph=pick[7].strip()
        dist=float(pick[-1].strip())
        #get all S picks
        if ph=='S' and dist > 10:
            print(ph)
            Sctr+=1
            ch=pick[1].strip()

            at=[float(num.strip()) for num in pick[8:-1]]
            at.append((at[-1]-int(at[-1]))*1e6)
            at=[int(num) for num in at]
            at=UTCDateTime(*at)
            starttime=at-20
            endtime=at+40
            succ=False
            if '?' in ch:
                #skip this
                continue
                #Many (perhaps most) of the picks have no assigned channel.
                #use the iris client's ability to parse lists of channels to get
                #an appropriate waveform.

                #compute back-azimuth
                stalong=float(pick[4].strip())
                stalat=float(pick[5].strip())

                baz=distaz(evlat,evlong,stalat,stalong)[2]

                #get appropriate N/E waveforms to rotate

                #TODO convert this try/except block to use the IRIS client
                try:
                    #search for N/1 channel with high-frequency response
                    inv=irisclient.get_stations(starttime=starttime,endtime=endtime,station=st,level="channel",channel="BHN,BH1,HHN,HH1,SHN,SH1")
                    #find appropriate channels in the inventory (ASDF returns list, not inventory type)
                    chan=inv.get_contents()['channels'][0].split('.')

                    stream=irisclient.get_waveforms(chan[0],chan[1],chan[2],chan[3],starttime,endtime)
                    wf1=stream[0]
                    if wf1.stats['channel'][-1] == '1':
                        wf2=irisclient.get_waveforms(chan[0],chan[1],chan[2],chan[3][0:-1]+'2',starttime,endtime)[0]
                        wfz=irisclient.get_waveforms(chan[0],chan[1],chan[2],chan[3][0:-1]+'Z',starttime,endtime)[0]

                        stream=Stream(traces=[wfz,wf1,wf2])
                        #get an inventory including the Z and E/2 channels
                        IRISinv=irisclient.get_stations(starttime=starttime,endtime=endtime,station=st,level="channel")
                        stream=stream.rotate(method="->ZNE",inventory=IRISinv)
                        
                    else:
                        wf2=irisclient.get_waveforms(chan[0],chan[1],chan[2],chan[3][0:-1]+'E',starttime,endtime)[0]
                        stream=Stream(traces=[wf1,wf2])
                    stream=stream.rotate(method="NE->RT",back_azimuth=baz)
                    wf=None
                    for trywf in stream:
                        if trywf.stats['channel'][-1]=='T':
                            wf=trywf
                            break
                    if wf and len(wf)>100:
                        wf.resample(100)
                        wf.detrend()
                        wf.normalize()
                        wf.write(saveDir+str(wfctr)+'_S.pkl',format="PICKLE")
                        np.save(saveDir+str(wfctr)+'_S.npy',wf.data)
                        succ=True
                        print "Got a rotated horizontal waveform"
                        
                        
                except Exception as e:
                    #handle data missing
                    print >> sys.stderr, e
                    stream=None
            else:
                wf=getWave(irisclient,st,ch,starttime,endtime,saveDir=saveDir)
                if wf and len(wf)>100:
                    wf.write(saveDir+str(wfctr)+'_S.pkl',format="PICKLE")
                    np.save(saveDir+str(wfctr)+'_S.npy',wf.data)
                    succ=True
            if succ:
                wfctr+=1
                pickwtr.writerow(pick)
                        
            

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
