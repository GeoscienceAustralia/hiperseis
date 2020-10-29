#script to collect waveforms corresponding to the automated S-wave picks and save them with an index to a directory
import sys
sys.path.append('/g/data1a/ha3/rlt118/hiperseis/')
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
from autopicks_approx import pickLoader
from obspy.core import UTCDateTime, Stream

fds=FederatedASDFDataSet('/g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt', logger=None)


import numpy as np
import re

#global method to get a waveform from a pick
def getStream(pick,fds,mintime):
    if pick[0][0]=='#': #this line is commented
        return None
    net=pick[6]
    st=pick[7]
    ch=pick[8]
    time=UTCDateTime(float(pick[9]))
    if mintime:
        if time < mintime:
            return None

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
            wf2=fds.get_waveforms(net,st,chan[2],chan[3][0:-1]+'E',starttime,endtime)[0]
            #this is a hack and wrong. Hopefully 1 and 2 are only off from N and E by a few degrees max
            wf1.id=wf1.id[:-1]+'N'
            wf2.id=wf2.id[:-1]+'E'
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
        return wf
    else:
        print >> sys.stderr, 'Waveform not found. Skipping pick...'
        return None



loadDir='/g/data/ha3/rlt118/neural-datasets/autopicks/'
saveDir='/g/data/ha3/rlt118/lowslope-autopicks/'
with open(loadDir+'ensemble.s.txt','r') as f:
    picks=f.readlines()
    pickCtr=len(picks)-1

with open(loadDir+'ensemble.s.verified.txt', 'r') as f:
    vpicks = f.readlines()
    vpicks = vpicks[1:]

vpdict = {}
for pick in vpicks:
    plist = pick.split()
    slope = float(plist[20])
    if slope < 1:
        evID = plist[0]
        net = plist[6]
        st = plist[7]
        ch = plist[8]
        vpdict[evID+net+st+ch] = 1
    

outfile = saveDir + "ensemble.s.unverified.txt"
with open(outfile, "w") as f:
    for pick in picks:
        plist = pick.split()
        cwt = float(plist[18])
        slope = float(plist[20])
        if slope < 1:
            evID = plist[0]
            net=plist[6]
            st=plist[7]
            ch=plist[8]
            if evID+net+st+ch in vpdict:
                s = getStream(plist, fds, None)
                if s:
                    #save mseed
                    cleanevID = re.sub('[^A-Za-z0-9]+', '', evID)
                    s.write(saveDir+cleanevID+net+st+ch+".mseed", format="MSEED")
                    #write pick to out file
                    f.write(pick+'\n')



