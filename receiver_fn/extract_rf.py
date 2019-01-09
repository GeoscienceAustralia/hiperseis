import numpy as np
from scipy.signal import hilbert
import rf
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
def phase_weights(stream):

    tphase=[]
    for tr in stream: 
            analytic = hilbert(tr.data) 
            angle = np.angle(analytic) 
            iPhase = np.exp(1j * angle) 
            tphase.append(iPhase) 
    tphase=np.array(tphase)
    tphase=np.abs(np.mean(tphase,axis=0))
    return tphase/np.max(tphase)

def count_groups(stream):
    group=[]
    for trace in stream:
        group.append(trace.stats.rf_group)
    group=np.array(group)
    return group

#-------------Main---------------------------------

if __name__=='__main__':

    ''' @package extract_rf
    This code contains different approaches to extract RFs from H5 file in stacked form.
    Output is preapared for trans-dimensional inversion in ASCII format

    Currently there are two methods of stacking
    1. rf stacked by similarity 
    2. all rf stacked 
     
    '''

    print "Reading the input file..."
    # Input file
    stream=rf.read_rf('DATA/7X-ZRT-R-cleaned.h5','H5')
    print "Reading is done..."

    station_list=[]
    group_list=[]
    # here we collect station names

    for i in xrange(len(stream)):
        station_list.append(stream[i].stats.station.encode('utf-8'))
        group_list.append(stream[i].stats.rf_group)

    group_list=np.array(group_list)
    station_list=np.array(station_list)

    # we need to find the largest number of groups for each uniqe station
    gidx=np.argsort(-group_list)
    group_list=group_list[gidx]
    station_list=station_list[gidx]

    # unique will return first occurence of the station sorted in descending order of group number
    station_list,idx=np.unique(station_list,return_index=True)
    group_list=group_list[idx]

    print "Gathered ",len(station_list)," stations"
    for i in xrange(station_list.shape[0]):
        print station_list[i],group_list[i] 

    estat=''
    while station_list[estat==station_list].shape[0]==0:
          estat=raw_input("Station to extract: ")
    
    net=stream[0].stats.network.encode('utf-8')

    filter_type='bandpass'
    freqmin=0.03
    freqmax=0.5

    station=stream.select(station=estat,component='R').moveout()

#we use a zero-phase-shift band-pass filter using 2 corners. This is done in two runs forward and backward, so we end up with 4 corners de facto.
    station=station.filter(filter_type, freqmin=freqmin, freqmax=freqmax,corners=2,zerophase=True).interpolate(10).trim2(-5.,20.)
# somehow trim2 doesn't trim down to 20 sec after arrival. Should be investigated further
#   print station[0].stats.delta,station[0].stats.npts
    

    if len(station)>1:


       # first we get stacks - normal and phase weighted
       phase_w=phase_weights(station)
       copy_st=station.copy()
       stacked=station.copy().stack()
       ph_weighted=copy_st.stack()
       ph_weighted=ph_weighted[0].data*phase_w
       time=stacked[0].stats.delta*np.array(list(xrange(stacked[0].stats.npts)))


       # then we take the same for each similarity groups

       groups=count_groups(station)
       max_grp=np.max(groups)
       print "Max grp ",max_grp
       # however first we define general plotting scheme and plot previous results
       fig = plt.figure(figsize=(11.69,8.27))
       columns=2
       rows=np.int(np.ceil(float(max_grp)/float(columns)))+1
       grid=gridspec.GridSpec(columns,rows,wspace=0.2,hspace=0.2)
       ax=plt.subplot(grid[0])
       ax.plot(time,stacked[0].data)
       ax.set_title('Stacked')
       ax=plt.subplot(grid[1])
       ax.plot(time,ph_weighted)
       ax.set_title('Phase weighted stack')
     
       frame=2
       for i in xrange(max_grp):

             grp_stream=rf.RFStream()

             for trace in station:
                 if trace.stats.rf_group==i:
                    grp_stream.append(trace)
             
             grp_stacked=grp_stream.copy().stack()
             phase_w=phase_weights(grp_stream)
             grp_stacked_wght=grp_stacked.copy()[0].data*phase_w
             grp_time=grp_stacked[0].stats.delta*np.array(list(xrange(grp_stacked[0].stats.npts)))
             ax=plt.subplot(grid[i+frame])
             ax.plot(grp_time,grp_stacked_wght)
             ax.set_title('Group '+str(i))

 

       plt.show()

       
       with open(net+'-'+estat+'-rf.txt','w') as text_file:
            for i in xrange(time.shape[0]):
                  text_file.write(str(time[i]-5)+'   '+str(ph_weighted[i])+'\n')

       text_file.close()
       
       
