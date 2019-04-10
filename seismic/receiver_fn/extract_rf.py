import numpy as np
from scipy.signal import hilbert
import rf
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
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
     
    Note the parameters of gaussian pulse and its width where

Value of "a" | Frequency (hz) at which G(f) = 0.1 |  Approximate Pulse Width (s)

10                      4.8                                0.50
5                       2.4                                0.75
2.5                     1.2                                1.00
1.25                    0.6                                1.50
1.0                     0.5                                1.67 (5/3)
0.625                   0.3                                2.10
0.5                     0.24                               2.36
0.4                     0.2                                2.64
0.2                     0.1                                3.73

    '''

    print "Reading the input file..."
    # Input file
#   stream=rf.read_rf('DATA/7X-ZRT-R-1cleaned.h5','H5')
    stream=rf.read_rf('DATA/7X-ZRT-R-ma12-cleaned.h5','H5')
    print "Reading is done..."

    net=stream[0].stats.network.encode('utf-8')
    # output directory
    out_dir=net+"-INV/"

    # inversion programs use 1Hz pulse width, therefore higher corner should be not lower than that
    filter_type='bandpass'
    freqmin=0.1
    freqmax=1.0

    # Trimming window
    tstart=-5.
    tend=40.



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
    sstat=[]

#   while station_list[estat==station_list].shape[0]==0:
#         estat=raw_input("Station to extract: ")
    estat=raw_input("Station to extract [All]: ")
    if station_list[estat==station_list].shape[0]==0: 
       sstat=station_list
       plot=False
    else:
       sstat.append(estat)
       plot=True


    
    for estat in sstat:

         station=stream.select(station=estat,component='R').moveout()

#        we use a zero-phase-shift band-pass filter using 2 corners. This is done in two runs forward and backward, so we end up with 4 corners de facto.
#        print station[0].stats.delta,station[0].stats.npts
         
            

         if len(station)>1:

            for trace in station:
                # preserve original amplitude to rescale later to preserve proportions relative to source
                if trace.stats.amax > 0:
                   amp_max=trace.stats.amax
                   print "*"
                else:
                   amp_max=np.max(trace.data)
                trace.taper(0.01)
                # 6.25 is the frequency hardwired into the inversion program
                trace=trace.filter(filter_type, freqmin=freqmin, freqmax=freqmax,corners=2,zerophase=True).interpolate(6.25)
#               trace=trace.interpolate(6.25)
                trace.data=trace.data*(amp_max/np.amax(trace.data))


            # first we get stacks - normal and phase weighted
            copy_st=station.copy()
            stacked=station.copy().stack()
            stacked.trim2(tstart,tend,'onset')
            time_s=stacked[0].stats.delta*np.array(list(xrange(stacked[0].stats.npts)))+tstart

            amp_max=np.max(stacked[0].data)

            phase_w=phase_weights(station)
            ph_weighted=copy_st.stack()
            ph_weighted[0].data=ph_weighted[0].data*phase_w
            # Note - weighting changes the real amplitude and it must be rescaled back to origin
            ph_weighted.trim2(tstart,tend,'onset')
            time_p=ph_weighted[0].stats.delta*np.array(list(xrange(ph_weighted[0].stats.npts)))+tstart
            zero=ph_weighted[0].data[time_p<0.]
            idx=np.max(np.where(zero<=0.)[0])
            ph_weighted[0].data[:idx+1]=0.
#           ph_weighted.filter(filter_type, freqmin=freqmin, freqmax=freqmax,corners=1,zerophase=True)
            ph_weighted[0].data=ph_weighted[0].data*(amp_max/np.max(ph_weighted[0].data))


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
            ax.plot(time_s,stacked[0].data)
            ax.set_title(estat+' Stacked')
            ax=plt.subplot(grid[1])
            ax.plot(time_p,ph_weighted[0].data)
            ax.set_title('Phase weighted stack')
          
            frame=2
            for i in xrange(max_grp):

                  grp_stream=rf.RFStream()

                  for trace in station:
                      if trace.stats.rf_group==i:
                         grp_stream.append(trace)
                  print "Group: ",i," number of records: ", len(grp_stream)
                  grp_stacked=grp_stream.copy().stack()
                  grp_stck_max=np.max(np.abs(grp_stacked.copy()[0].data))
#                 grp_stck_max=amp_max
                  phase_w=phase_weights(grp_stream)
                  grp_stacked_wght=grp_stacked.copy()[0].data*phase_w
                  grp_stacked_wght=grp_stacked_wght*(grp_stck_max/np.max(np.abs(grp_stacked_wght)))  
                  
                  grp_time=grp_stacked[0].stats.delta*np.array(list(xrange(grp_stacked[0].stats.npts)))+tstart
                  ax=plt.subplot(grid[i+frame])
                  ax.plot(grp_time,grp_stacked_wght)
                  ax.set_title('Group '+str(i))

 


            if not os.path.exists(out_dir):
                  os.makedirs(out_dir) 
                  os.makedirs(out_dir+'PDF')

            if plot:
                  plt.show()
            else:
                  fig.savefig(out_dir+'PDF/'+net+'-'+estat+'-rf2-ph_weighted.pdf',format='PDF')
                  plt.close('all')

            
            with open(out_dir+net+'-'+estat+'-rf2-ph_weighted.dat','w') as text_file:
                  for i in xrange(time_p.shape[0]):
                        text_file.write(str(time_p[i])+'   '+str(ph_weighted[0].data[i])+'\n')

            text_file.close()
       
