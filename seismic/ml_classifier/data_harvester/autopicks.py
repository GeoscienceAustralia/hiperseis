#loads the autopicks from Marcus and returns an iterator for autopick data
from obspy.core import UTCDateTime
from obspy.core.stream import Stream

import re

import sys

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
        if pick[0][0]=='#':
            #this is a commented line
            return self.next()
        
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

                inv=self.fds.get_stations(starttime,endtime,network=net,station=st)
                #find appropriate channels in the inventory (ASDF returns list, not inventory type)
                for chan in inv:
                    if re.search('^[HBS].[N1]$',chan[3]):
                        break
                if not re.search('^[HBS].[N1]$',chan[3]):
                    raise Exception('no appropriate horizontal channels found')

                stream=self.fds.get_waveforms(net,st,chan[2],chan[3],starttime,endtime)
                wf1=stream[0]
                if wf1.stats['channel'][-1] == '1':
                    wf2=self.fds.get_waveforms(net,st,chan[2],chan[3][0:-1]+'2',starttime,endtime)[0]
                    wfz=self.fds.get_waveforms(net,st,chan[2],chan[3][0:-1]+'Z',starttime,endtime)[0]

                    stream=Stream(traces=[wfz,wf1,wf2])
                    #this is a necessary hack - ASDF inventories do not have orientation data for
                    #the 1 and 2 channels so it is impossible to rotate to ZNE accurately. Therefore
                    #we get station data off IRIS. All the OA installations have ZNE data provided,
                    #so only the permanent stations need to rotate 12->NE, and being permanent they should
                    #have metadata available through IRIS.
                    IRISinv=self.ic.get_stations(network=net,station=st,channel=ch,level="channel")
                    stream=stream.rotate(method="->ZNE",inventory=IRISinv)
                    
                else:
                    wf2=self.fds.get_waveforms(net,st,chan[2],chan[3][0:-1]+'E',starttime,endtime)[0]
                    stream=Stream(traces=[wf1,wf2])
            except Exception as e:
                #handle data missing
                print >> sys.stderr, e
                stream=None
                
            wf=None
            if stream:
                stream=stream.rotate(method='NE->RT',back_azimuth=baz)
                for trywf in stream:
                    if trywf.stats['channel'][-1]=='T':
                        wf=trywf
                        break

        else:
            try:
                #just grab the channel that the pick was on
                inv=self.fds.get_stations(starttime,endtime,network=net,station=st,channel=ch)
                loc=inv[0][2]
                wf=self.fds.get_waveforms(net,st,loc,ch,starttime,endtime)[0]
            except Exception as e:
                print >> sys.stderr, e
                wf=None
        if wf and len(wf)>100:
            wf.resample(100)
            wf.detrend()
            wf.normalize()
            ret=[pick,wf.data]
        else:
            print >> sys.stderr, 'Waveform not found. Skipping pick...'
            ret=self.next()
        
        return ret



