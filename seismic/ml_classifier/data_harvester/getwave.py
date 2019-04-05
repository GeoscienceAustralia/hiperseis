import numpy as np


#ASDF database library
#from ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet

#initialise waveform dataset
#fds = FederatedASDFDataSet('/g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt', variant='db', use_json_db=True,
#                           logger=None)





wfctr=0

#define a method to take pick parameters and return the corresponding trace
def _genTS(fdsnclient,st,ch,starttime,endtime,loc=None,net=None):
    
    #get the trace

    try:
        if not net:
            statinfo=fdsnclient.get_stations(starttime=starttime,endtime=endtime,station=st,channel=ch,loc=loc,level="channel").get_contents()['channels'][0].split('.')
            net=statinfo[0]
        if not loc:
            loc=statinfo[2]
        waveforms=fdsnclient.get_waveforms(net,st,loc,ch,starttime,endtime)
        if len(waveforms)>0:
            print "Got an IRIS waveform!"
            ret=waveforms[0]
        else:
            ret=None
    except:
        ret=None
    

    if not ret:
        print "No waveform found."

    return ret

#external-facing function. This thing does all the required processing with obspy routines (detrend, normalise,
#resample) and returns the preprocessed trace object.

#'ISC' keyword argument - is the channel string in 'ISC format' (no location)? If not, split the channel string
#up using an underscore. Splitting is required for, e.g., Babak's GA pick database.

def getWave(webclient,st,chloc,starttime,endtime,saveDir="/g/data/ha3/rlt118/neural-datasets/categoriser-teleseismic/",phase='S',ISC=True,network=None):
    global wfctr
    print starttime
    st=st.strip()
    if '?' in chloc:
        return False

    if ISC:
        ch=chloc.strip()
        wf=_genTS(webclient,st,ch,starttime,endtime)
    else:
        chloc=chloc.split('_')
        ch=chloc[0].strip()
        loc=chloc[1].strip() if len(chloc)>1 else ''
        wf=_genTS(webclient,st,ch,starttime,endtime,loc=loc,net=network)
    #do some checks to see if the returned waveform is not empty
    if wf is not None and len(wf)>100:
        #resample to 100 Hz, detrend and normalise
        wf.resample(100)
        wf.detrend()
        wf.normalize()
        wfctr+=1
        wf.write(saveDir+str(wfctr)+'_'+phase+'.pkl',format="PICKLE")
        np.save(saveDir+str(wfctr)+'_'+phase+'.npy',wf.data)
        return True
    else:
        return False
    
