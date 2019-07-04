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


#load GA's pick database
GApicks=loadmat("GA.mat")['GA']['picks']

#define a method to take P and S pick parameters and return the corresponding trace
def genTS(net,st,ch,loc,ptime,stime):
    
    #add some fuzz to the start and end times of the traces to avoid the neural net training itself to just pick a constant time
    starttime=ptime-(stime-ptime)/2+np.random.uniform(-(stime-ptime)/4,(stime-ptime)/4)
    endtime=stime+(stime-ptime)-np.random.uniform(-(stime-ptime)/2,(stime-ptime)/2)

    waveforms=fds.get_waveforms(net, st, loc, ch, 
                      starttime, endtime,
                      automerge=True, trace_count_threshold=10)
    if len(waveforms)==1: #discard empty streams or those containing multiple traces
        ret=waveforms[0]
    else:
        ret=None
    return ret
        
#remove the orientation code from channel data and return resulting location and generalised channel
def stripChannel(channel):
    chloc=channel.split('_')
    chloc[0]=chloc[0].rstrip()[0:-1]
    return '_'.join(chloc)



#look for individual channels with both P and S picks
simulctr=0

for event,phases in enumerate(GApicks['ph']):
    for index,phase in enumerate(phases):
        #if found an S pick
        if phase=='S':
            #search for a matching P pick
            
            nel,stl,chl=GApicks['ne'][event],GApicks['st'][event],GApicks['ch'][event],
            neS,stS,chS=nel[index],stl[index],chl[index]
            #look for any component of the same channel (N/E/Z)
            chS_strip=stripChannel(chS)
            for trialPind in range(len(nel)):
                if phases[trialPind]=='P' and nel[trialPind]==neS and stl[trialPind]==stS and chS_strip==stripChannel(chl[trialPind]):
                    
                    simulctr +=1
                    
                    break

print simulctr
            
"""

#process the dictionary of multi-phase picks to extract all the data required to query the waveform database
print('Processing waveforms for all multi-phase picks...')
wfctr=0
timediffs=[]
for ch in allpickdict:
    #take the string we created before and convert it back into a list of network, station, channel parameters
    idArray=ch.split()
    net=idArray[0]
    st=idArray[1]
    chloc=idArray[2].split('_')#the pick database has the channel and location as one string e.g. 'BHZ_00', or just 'BHZ' if the location is '--'.
    chan=chloc[0]
    if len(chloc)>1:
        loc=chloc[1]
    else:
        loc='--'
    for pick in allpickdict[ch]:
        if pick[0][0]=='P':
            ptimeint=[int(i) for i in pick[1][0]]
            pmsec=int((pick[1][0][-1]-ptimeint[-1])*10**6)#microsecond component of the time, gets truncated when converting to int
            ptime=UTCDateTime(*(ptimeint+[pmsec]))
            stimeint=[int(i) for i in pick[1][1]]
            smsec=int((pick[1][1][-1]-stimeint[-1])*10**6)
            stime=UTCDateTime(*(stimeint+[smsec]))
        else:
            ptimeint=[int(i) for i in pick[1][1]]
            pmsec=int((pick[1][1][-1]-ptimeint[-1])*10**6)
            ptime=UTCDateTime(*(ptimeint+[pmsec]))
            stimeint=[int(i) for i in pick[1][0]]
            smsec=int((pick[1][0][-1]-stimeint[-1])*10**6)
            stime=UTCDateTime(*(stimeint+[smsec]))
        #get a waveform to train against
        wf=genTS(net,st,chan,loc,ptime,stime)
        if wf is not None and len(wf)>100 and ptime>wf.times("utcdatetime")[0]:#discard bad or short waveforms extracted from the database
            wfctr+=1
            #resample the waveforms to 1024 points, detrend and normalise. The extra 0.01 ensures that the resulting trace
            #does in fact have 1024 points
            wf.resample(1024.01/wf.stats.npts*wf.stats.sampling_rate)
            wf.detrend()
            wf.normalize()
	        #generate the pick distributions and normalise them
            pdist=np.exp(-np.power(wf.times()-(ptime-wf.times("utcdatetime")[0]),2)/(0.02))#sigma is 0.1
            sdist=np.exp(-np.power(wf.times()-(stime-wf.times("utcdatetime")[0]),2)/(0.02))
            pdist=pdist/np.sum(pdist)
            sdist=sdist/np.sum(sdist)
            #anything that isn't P or S should be noise
            noisedist=np.ones(1024)-(pdist+sdist)
            dists=np.array([pdist,sdist,noisedist]).transpose()
            #plot the resulting final waveform to a file
            outfname='_'.join((net,st,chan,loc,ptime.ctime(),stime.ctime()))
            wf.plot(outfile='wftestplot/'+outfname+'.png')
            #save training data and metadata regarding the actual seismic station etc.
            #define a couple of extra attributes to store pick times in the trace object
            wf.ptime=ptime
            wf.stime=stime
            if wfctr<=750:
                #pickle the whole trace data. This should now include
                #the pick times wf.plot() at least shows channel and station data too so I guess
                #this is there somewhere, if not directly accessible.
                wf.write('trainset/'+str(wfctr)+'.pkl',format="PICKLE")
                np.save('trainset/'+str(wfctr)+'_trc.npy',wf.data)
                np.save('trainset/'+str(wfctr)+'_msk.npy',dists)
            else:
                #pickle the whole trace data. This should now include
                #the pick times wf.plot() at least shows channel and station data too so I guess
                #this is there somewhere, if not directly accessible.
                wf.write('valset/'+str(wfctr)+'.pkl',format="PICKLE")
                np.save('valset/'+str(wfctr)+'_trc.npy',wf.data)
                np.save('valset/'+str(wfctr)+'_msk.npy',dists)
print('successfully extracted '+str(wfctr)+' waveforms from the database.')

"""
