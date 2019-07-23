#loads the autopicks from Marcus and returns an iterator for autopick data
from obspy.core import UTCDateTime
from obspy.core.stream import Stream

import re

import sys

import numpy as np

#initialiser arguments:
#dataset - FederatedASDFDataSet
#irisclient - FDSN client for station metadata
#[fname] - file where picks are stored

class pickLoader:
    def __init__(self,dataset,irisclient,fname='/g/data/ha3/mwh547/Picks/20190320/ensemble.s.txt'):
        self.f=open(fname,'r')
        self.fds=dataset
        self.ic=irisclient

    def __iter__(self):
        return self    

    def next(self):
        pick=self.f.readline().split()
        if not pick:
            self.f.close()
            raise StopIteration
        data=_getData(pick,self.fds,self.ic)
        if data is not None:
            ret=[pick,data]
        else:
            ret=self.next()
        return ret

class pickLoaderRand:
    def __init__(self,dataset,irisclient,fname='/g/data/ha3/mwh547/Picks/20190320/ensemble.s.txt'):
        #read the whole file
        with open(fname,'r') as f:
            self.picks = [line.split() for line in f]
        self.fds=dataset
        self.ic=irisclient
    
    def getPick(self,ind): #get the pick on line number ind
        pick=self.picks[ind]
        data=_getData(pick,self.fds,self.ic)
        return [pick,data]
        

    def __len__(self): #length of the autopick db
        return len(self.picks)


#global method to get a waveform from a pick
def _getData(pick,fds,ic):
    if pick[0][0]=='#': #this line is commented
        return None
    net=pick[6]
    st=pick[7]
    ch=pick[8]
    time=UTCDateTime(float(pick[9]))

    starttime=time-20
    endtime=time+40



    if ch=='00T':

        baz=float(pick[14])
        #get all the horizontal short-period channels at loc '00' or '' (assuming the 00 in 00T refers to
        #loc)
        try:

            inv=fds.get_stations(starttime,endtime,network=net,station=st)
            #find appropriate channels in the inventory (ASDF returns list, not inventory type)
            for chan in inv:
                if re.search('^[HBS].[N1]$',chan[3]):
                    break
            if not re.search('^[HBS].[N1]$',chan[3]):
                raise Exception('no appropriate horizontal channels found')

            stream=fds.get_waveforms(net,st,chan[2],chan[3],starttime,endtime)
            wf1=stream[0]
            if wf1.stats['channel'][-1] == '1':
                wf2=fds.get_waveforms(net,st,chan[2],chan[3][0:-1]+'2',starttime,endtime)[0]
                wfz=fds.get_waveforms(net,st,chan[2],chan[3][0:-1]+'Z',starttime,endtime)[0]

                stream=Stream(traces=[wfz,wf1,wf2])
                #this is a necessary hack - ASDF inventories do not have orientation data for
                #the 1 and 2 channels so it is impossible to rotate to ZNE accurately. Therefore
                #we get station data off IRIS. All the OA installations have ZNE data provided,
                #so only the permanent stations need to rotate 12->NE, and being permanent they should
                #have metadata available through IRIS.
                IRISinv=ic.get_stations(network=net,station=st,channel=ch,level="channel")
                stream=stream.rotate(method="->ZNE",inventory=IRISinv)
                
            else:
                wf2=fds.get_waveforms(net,st,chan[2],chan[3][0:-1]+'E',starttime,endtime)[0]
                stream=Stream(traces=[wf1,wf2])
        except Exception as e:
            #handle data missing
            print >> sys.stderr, e
            stream=None
            
        wf=None
        if stream:
            try:
                stream=stream.rotate(method='NE->RT',back_azimuth=baz)
                for trywf in stream:
                    if trywf.stats['channel'][-1]=='T':
                        wf=trywf
                        break
            except Exception as e:
                print >>sys.stderr, e
                wf=None

    else:
        try:
            #just grab the channel that the pick was on
            inv=fds.get_stations(starttime,endtime,network=net,station=st,channel=ch)
            loc=inv[0][2]
            wf=fds.get_waveforms(net,st,loc,ch,starttime,endtime)[0]
        except Exception as e:
            print >> sys.stderr, e
            wf=None
    if wf and len(wf)>100:
        wf.resample(100)
        wf.detrend()
        wf.normalize()
        return wf.data
    else:
        print >> sys.stderr, 'Waveform not found. Skipping pick...'
        return None
