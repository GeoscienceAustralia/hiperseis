#script to import and process waveforms containing both P-wave and S-wave
#picks from analysts. Produces a set of files which could be used to train a
#machine-learning algorithm to recognise P and S waves from seismic traces.
from obspy.core import *
from mat4py import *
import numpy as np
import sys
#ASDF database library
sys.path.append('/g/data1a/ha3/rakib/seismic/pst/passive-seismic')
from ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet

#initialise waveform dataset
fds = FederatedASDFDataSet('/g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt', variant='db', use_json_db=True,
                           logger=None)
import os

#load GA's pick database
GApicks=loadmat("GA.mat")['GA']['picks']

#define a method to take pick parameters and return the corresponding trace
def genTS(net,st,ch,loc,starttime,endtime):
    
    #get the trace
    try:
        waveforms=fds.get_waveforms(net, st, loc, ch, 
                      starttime, endtime,
                      automerge=True, trace_count_threshold=10)
        if len(waveforms)==1: #discard empty streams or those containing multiple traces
            ret=waveforms[0]
        else:
            ret=None
    except:
        ret=None
    return ret


saveDir='/g/data/ha3/rlt118/neural-datasets/categoriser-local/'

#look for individual channels with S picks
Sctr=0
pSctr=0

for event,phases in enumerate(GApicks['ph']):
    for index,phase in enumerate(phases):
        #if found an S pick
        if phase=='P':
            Stime=[int(num) for num in GApicks['at'][event][index]]
            Stime=Stime+[int((GApicks['at'][event][index][-1]-Stime[-1])*1e6)]
            Stime=UTCDateTime(*Stime)
            #save a minute around the S pick
            starttime=Stime-20
            endtime=Stime+40

            #get the waveform
            Sne=GApicks['ne'][event][index].strip()
            Sst=GApicks['st'][event][index].strip()
            Schloc=GApicks['ch'][event][index].split('_')
            Sch=Schloc[0].strip()
            if len(Schloc)>1:
                Sloc=Schloc[1].strip()
            else:
                Sloc=''
            wf=genTS(Sne,Sst,Sch,Sloc,starttime,endtime)
            pSctr+=1
            #do some checks to see if the returned waveform is not empty
            if wf is not None and len(wf)>100:
                #resample to 100 Hz, detrend and normalise
                wf.resample(100)
                wf.detrend()
                wf.normalize()
                Sctr+=1
                os.system('clear')
                print len(wf)
                wf.write(saveDir+str(Sctr)+'_P.pkl',format="PICKLE")
                np.save(saveDir+str(Sctr)+'_P.npy',wf.data)

            print str(Sctr)+'/'+str(pSctr)
                    

print str(Sctr)+'/'+str(pSctr),"waveforms harvested"
