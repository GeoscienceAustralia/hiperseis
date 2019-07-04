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

#define a method to take P and S pick parameters and return the corresponding trace
def genTS(net,st,ch,loc,starttime,endtime):
    
    #get the trace
    print "querying database..."
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
    print "done"
    return ret


#look for individual channels with both P and S picks
Sctr=0
pSctr=0
Pctr=0
pPctr=0

for event,phases in enumerate(GApicks['ph']):
    for index,phase in enumerate(phases):
        #TODO harvest some not-S picks from similar regional distances, which should be far enough away from
        #the S picks at these distances that we can harvest a trace which doesn't include the S wave as well
        
        #focus on regional events for now. Many of the S picks are local, only one or two degrees away
        
        #if found a P pick (and the dataset is not swamped with P picks)
        if phase=='P' and Pctr < Sctr:

            if GApicks['dist'][event][index][0]>5 and GApicks['dist'][event][index][0]<30:
                #convert the time in GA.mat to a format that is compatible with
                #the UTCDateTime class of obspy, which requires partial seconds to be expressed as
                #integer-valued microseconds
                Stime=[int(num) for num in GApicks['at'][event][index]]
                Stime=Stime+[int((GApicks['at'][event][index][-1]-Stime[-1])*1e6)]
                Stime=UTCDateTime(*Stime)
                starttime=Stime-20
                endtime=Stime+40

                #get the waveform. Stripping the whitespace is necessary to actually match anything
                #when querying the ASDF database
                Sne=GApicks['ne'][event][index].strip()
                Sst=GApicks['st'][event][index].strip()
                Schloc=GApicks['ch'][event][index].split('_')
                Sch=Schloc[0].strip()
                if len(Schloc)>1:
                    Sloc=Schloc[1].strip()
                else:
                    Sloc=''
                wf=genTS(Sne,Sst,Sch,Sloc,starttime,endtime)
                pPctr+=1
                #do some checks to see if the returned waveform is not empty and actually contains the P arrival
                if wf is not None and len(wf)>100:
                    #resample to 100 Hz, detrend and normalise
                    wf.resample(100)
                    wf.detrend()
                    wf.normalize()
                    Pctr+=1
                    wf.write('data2/'+str(Pctr)+'_P.pkl',format="PICKLE")
                    np.save('data2/'+str(Pctr)+'_P.npy',wf.data)
                    
        elif phase=='S':        
            if GApicks['dist'][event][index][0]>5 and GApicks['dist'][event][index][0]<30:
                #convert the time in GA.mat to a format that is compatible with
                #the UTCDateTime class of obspy, which requires partial seconds to be expressed as
                #integer-valued microseconds
                Stime=[int(num) for num in GApicks['at'][event][index]]
                Stime=Stime+[int((GApicks['at'][event][index][-1]-Stime[-1])*1e6)]
                Stime=UTCDateTime(*Stime)
                starttime=Stime-20
                endtime=Stime+40

                #get the waveform. Stripping the whitespace is necessary to actually match anything
                #when querying the ASDF database
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
                #do some checks to see if the returned waveform is not empty and actually contains the S arrival
                if wf is not None and len(wf)>100:
                    #resample to 100 Hz, detrend and normalise
                    wf.resample(100)
                    wf.detrend()
                    wf.normalize()
                    Sctr+=1
                    wf.write('data2/'+str(Sctr)+'_S.pkl',format="PICKLE")
                    np.save('data2/'+str(Sctr)+'_S.npy',wf.data)
                
                

print str(Sctr)+'/'+str(pSctr),"waveforms harvested"
