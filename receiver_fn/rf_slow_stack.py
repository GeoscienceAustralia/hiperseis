#!/usr/bin/env python
import numpy as np
import matplotlib
# comment out the line below if you want to visualise the figure
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib.pyplot as plt
from obspy.core.utcdatetime import UTCDateTime as UTC

from joblib import Parallel, delayed
import functools
import matplotlib.tri as tri
import numpy.ma as ma
import os

# Here are the libraries to deal with RFSTREAM, it uses obspy classes for event and station
#from rf import RFStream
import rf
from obspy.core.event.event import Event
from obspy.core.inventory.station import Station
from obspy.core.util import AttribDict

from rf.profile import profile
from tqdm import tqdm
from rf.imaging import plot_profile_map
from rf import get_profile_boxes, iter_event_data, IterMultipleComponents
from obspy.signal.array_analysis import array_processing

from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize

from obspy.imaging.cm import obspy_sequential
from matplotlib import cm
import matplotlib.dates as mdates

import matplotlib.gridspec as gridspec

KM_PER_DEG = 111.1949


# Definition of the simple 1D Earth model, rememebr each interface will give one Ps conversion
# you can add shallow layer to study Ps conversions in sedimentary basins
z=np.array( [0,   35,  35,  165])
vp=np.array([5.8 ,6.5 ,8.04,8.17])
vs=np.array([3.5,3.9,4.48,4.51])
simple_model=rf.simple_model.SimpleModel(z,vp,vs)

# For sake of simplicity we will generate separate model for sedimentary basin 
zb=np.array( [0,   5,   5,   300])
vp=np.array([4.5 ,4.5,6.5 , 8.17])
vs=np.array([2.5, 2.5,3.5,  4.51])
basin_model=rf.simple_model.SimpleModel(zb,vp,vs)



# its not advised to use standard models because they create large number of conversions at each layer
#simple_model=rf.simple_model.load_model(fname='iasp91')


#-------------Main---------------------------------

if __name__=='__main__':


    ''' This program composes vespagrams to identify RF converted phases and their multiples
        please refer to Tian et al. GRL 2005 VOL. 32, L08301, doi:10.1029/2004GL021885 for good examples
        
        input  - H5 file with receiver functions
        output - PDF files to print

        Dependencies - rf and obspy packages beside other standard python packages
        The image is composed using triangulation. It gives good results but block median or mean must be implemented at some stage to reduce size of PDF.
    '''

    stream=rf.read_rf('DATA/7X-rf_zrt_cleaned.h5','H5')
    rf_type='LQT-Q '
    filter_type='bandpass'
    freqmin=0.03
    freqmax=0.5
#we use a zero-phase-shift band-pass filter using 2 corners. This is done in two runs forward and backward, so we end up with 4 corners de facto.
# Lets assume we have LQT orientation
    selected_stream=stream.select(component='Q').filter(filter_type, freqmin=freqmin, freqmax=freqmax,corners=2,zerophase=True).interpolate(10)


# if none lets try ZRT
    if len(selected_stream)<=0:
         selected_stream=stream.select(component='R').filter(filter_type, freqmin=freqmin, freqmax=freqmax,corners=2,zerophase=True).interpolate(10)

         rf_type='ZRT-R '

    if len(selected_stream)<=0:
         print "Tried Q and R components but neither faund, quitting..."
         exit(-1)
       

    station_list=[]

    # here we collect station names but maybe ID is more appropriate in case of having the same station names in different deployments
    
    for i in xrange(len(selected_stream)):
         station_list.append(selected_stream[i].stats.station.encode('utf-8'))
         net=selected_stream[i].stats.network.encode('utf-8')

    pdffile=net+'-'+rf_type.strip()+'-rf-vespagrams.pdf'

    exists=os.path.isfile(pdffile)

    if exists:
         # we add the time stamp to identify the file that can be read in linux command line as date -d @number (date -d @1542926631)
         pdffile=net+'-'+rf_type.strip()+'-'+str(int(UTC.now()._get_timestamp()))+'-rf-vespagrams.pdf'
       
    pdf=PdfPages(pdffile)
    pdf.attach_note(rf_type+filter_type+' '+str(freqmin)+'-'+str(freqmax)+' Hz') 
    d = pdf.infodict()
    d['Title'] = rf_type+'RF vespagrams of '+net+' network'
    d['Keywords'] =  rf_type+filter_type+' '+str(freqmin)+'-'+str(freqmax)+' Hz'

    station_list=np.unique(np.array(station_list))
    print "Gathered ",len(station_list)," stations"

#   Define layout of the page outer_grid
    columns=3
    rows=2
    frame=0
    figure=1
# Main loop here over all stations ------------------------------------------

    for i in xrange(station_list.shape[0]):
         if frame==0:
            printed=False 
            fig = plt.figure(figsize=(11.69,8.27))
            outer_grid=gridspec.GridSpec(columns,rows,wspace=0.2,hspace=0.2)
             
         print "Station ",station_list[i],i+1," of ",station_list.shape[0]
         traces=selected_stream.select(station=station_list[i])
         print 'Contains: ',len(traces),' events'

         # we choose short RF to simplify and speed up the processing
         # from -5 to 20 seconds and sloweness range from 5 to 9 s/deg
         # its enough to see multiples and possible LAB conversion at ~19 sec (~160km)

         traces=traces.trim2(-5,20,'onset')
         moved=[] 
         slow=[]

         for tr in traces:
                tr.normalize()

# This 'if' block is designed to check correct data placement on vespagram to trace the logic (debugging purposes)
#               if tr.stats.slowness > 6. and tr.stats.slowness < 7.:
#                  print 'altered'
#                  data=tr.data.copy()
#                  print data.shape,tr.stats.delta
#                  500 below corresponds to 0 with sampling rate of 100Hz
#                  data[500:800]=1.
#                  moved.append(data)
#               else:
                moved.append(tr.data.copy()/np.max(np.abs(tr.data)))
                slow.append(tr.stats.slowness)



         print "Sloweness min and max: ",np.min(slow),np.max(slow)
         slow.append(np.min(slow)-0.1)
         moved.append(np.zeros(traces[0].data.shape))
         slow.append(np.max(slow)+0.1)
         moved.append(np.zeros(traces[0].data.shape))


         slow=np.array(slow) 
         idx=np.argsort(slow)
         moved=np.nan_to_num(np.array(moved))
#        moved=np.array(moved)
         slow=slow[idx]
         z=moved[idx,:]
#        print 'minmax',np.min(z),np.max(z)
         x=np.array(list(range(moved.shape[1])))*traces[0].stats.delta-5.
         x=np.repeat([x],moved.shape[0],axis=0)
         y=np.ones((moved.shape[0],moved.shape[1]))

         phase_Ps=[]
         phase_Pms=[]
         phase_PpPmS=[]
         phase_PpSmS=[]
         phase_slow=[]
         # basin part
         phase_Pbs=[]
         phase_PpPbs=[]

         for j in xrange(slow.shape[0]):
                y[j,:]=y[j,:]*slow[j]
                phase_Ps.append(simple_model.calculate_delay_times(slow[j],phase='PS'))
                phase_Pms.append(simple_model.calculate_delay_times(slow[j],phase='PmS'))
                phase_PpPmS.append(simple_model.calculate_delay_times(slow[j],phase='PpPmS'))
                phase_PpSmS.append(simple_model.calculate_delay_times(slow[j],phase='PpSmS'))
                phase_slow.append(np.ones(phase_Ps[-1].shape[0])*slow[j])

                # basin, we will use reflection at the top layer only
                if len(zb)>0:
                   phase_Pbs.append(basin_model.calculate_delay_times(slow[j],phase='PS'))
                   phase_PpPbs.append(basin_model.calculate_delay_times(slow[j],phase='PpPmS'))
                   

         xi = np.linspace(-5,20,200)
         yi = np.linspace(0,9,400)
         
#        Gridding the data using triangulation. standard gridding doesn't work well here
         triang = tri.Triangulation(x.flatten(), y.flatten())
         interpolator = tri.LinearTriInterpolator(triang, z.flatten())
         xi, yi = np.meshgrid(xi, yi)
         zi = interpolator(xi, yi)
         
#        Define two plots as inner_grid to place them inside one cell of outer_grid
         inner_grid=gridspec.GridSpecFromSubplotSpec(2,1,subplot_spec=outer_grid[frame],wspace=0.5,hspace=0.)
         ax1=plt.Subplot(fig,inner_grid[0])
         ax2=plt.Subplot(fig,inner_grid[1],sharex=ax1)
         lim=np.max(np.abs(zi[zi<0])*0.5)
         levels = np.linspace(-lim,lim,15)
#        print "Levels ",-lim,lim
         cmap= plt.cm.jet
         cs=ax1.contourf(xi,yi,zi,levels=levels,extend='both',cmap=cmap)
         cs.cmap.set_under('k')
         cs.cmap.set_over('k')
         ax1.set_ylim(5,9)
         ax1.set_xlim(-5,20)
         ax1.plot(phase_Ps,slow,color='black')      # direct conversion, positive amplitude
         ax1.plot(phase_PpPmS,slow,color='crimson') # multiples,         positive amplitude
         ax1.plot(phase_PpSmS,slow,color='purple')  # multiples,         negative amplitude
         
         ax1.annotate('Pms', xy=(phase_Ps[-1][0],9.1),xycoords='data',ha='center', va='bottom', \
                rotation=0.,annotation_clip=False,fontsize=7)
         ax1.annotate('Ps LAB', xy=(phase_Ps[-1][-1],9.1),xycoords='data',ha='center', va='bottom', \
                rotation=0.,annotation_clip=False,fontsize=7)

         if len(phase_Pbs)>0:
                ax1.annotate('Pbs', xy=(phase_Pbs[-1][0],9.1),xycoords='data',ha='center', va='bottom', \
                rotation=0.,annotation_clip=False,fontsize=7)
                ax1.plot(phase_Pbs,slow,color='black')
                ax1.plot(phase_PpPbs,slow,color='crimson')
            


         ax1.spines['bottom'].set_visible(False)
         ax1.tick_params(labelbottom='off')
         ax1.spines['bottom'].set_visible(False)
         ax1.yaxis.tick_right()
         ax1.yaxis.set_label_position("right")
         xlabels = ax1.get_xticklabels()
         ylabels = ax1.get_yticklabels()
         for label in xlabels:
             label.set_rotation(90)
             label.set_fontsize(7)
         for label in ylabels:
             label.set_rotation(90)
             label.set_fontsize(7)

         ax1.annotate(station_list[i], xy=(-0.08,0), ha='left', va='center', \
             xycoords='axes fraction', textcoords='offset points',rotation=90.)

         start, end = ax1.get_ylim()
         ax1.yaxis.set_ticks(np.arange(start+1, end+1, 1))

         cs=ax2.contourf(xi,-1.*yi,zi,levels=levels,extend='both',cmap=cmap)
         cs.cmap.set_under('k')
         cs.cmap.set_over('k')
         ax2.spines['top'].set_visible(False)
         ax2.set_ylim(-9,-5)
         ax2.set_xlim(-5,20)
         ax2.yaxis.tick_right()
         ax2.yaxis.set_label_position("right")
         ylabels = ax2.get_yticklabels()
         ax2.plot(phase_Ps,-slow,color='black')
         ax2.plot(phase_PpPmS,-slow,color='crimson')
         ax2.plot(phase_PpSmS,-slow,color='purple')
        

         ax2.annotate('+PpPms', xy=(phase_PpPmS[-1][0],-9.1),xycoords='data',ha='center', va='top', \
                rotation=0., annotation_clip=False,fontsize=7,color='crimson')
         ax2.annotate('-PpSms', xy=(phase_PpSmS[-1][0],-9.1),xycoords='data',ha='center', va='top', \
                rotation=0., annotation_clip=False,fontsize=7,color='purple')

         if len(phase_PpPbs)>0:
                ax2.annotate('+PpPbs', xy=(phase_PpPbs[-1][0],-9.1),xycoords='data',ha='center', va='top', \
                rotation=0., annotation_clip=False,fontsize=7,color='crimson')
                
                ax2.plot(phase_PpPbs,-slow,color='crimson')
                ax2.plot(phase_Pbs,-slow,color='black')



         for label in ylabels:
             label.set_rotation(90)
             label.set_fontsize(7)

         if frame > 3:
             xlabels = ax2.get_xticklabels()
             for label in xlabels:
                 label.set_rotation(90)
                 label.set_fontsize(7)
             ax2.set_xlabel('Time (sec.)')
         else:
             ax2.set_xticklabels([])

         if frame%2!=0:
            ax2.annotate('Sloweness s/deg', xy=(1.2,1), ha='left', va='center', \
                xycoords='axes fraction', textcoords='offset points',rotation=90.)
                


         start, end = ax2.get_ylim()
         ax2.yaxis.set_ticks(np.arange(start, end, 1))
         traces.moveout()
         x=np.array(list(range(traces[0].data.shape[0])))*traces[0].stats.delta-5.
         y=traces.stack()
# Some amplitude scaling to have nice plot
         y=y[0].data/1.5-5.
         ax2.plot(x,y,clip_on=False,linewidth=3,color='white')
         ax2.plot(x,y,clip_on=False,linewidth=1)

         fig.add_subplot(ax1)
         fig.add_subplot(ax2)
         
         frame=frame+1
         print 'frame',frame
         if frame >= rows*columns:
            cb_ax=fig.add_axes([0.25,0.95,0.5 ,0.02])
            labels=fig.colorbar(cs,cax=cb_ax,ticks=[-1,0,1],orientation='horizontal',extend='neither',extendfrac=0.00001,extendrect=True,drawedges=False)
            labels.ax.set_xticklabels(['-','0','+',''])
            pdf.savefig()
            figure=figure+1
            frame=0
            printed=True
            plt.close()
#           plt.show()

if not printed:

         cb_ax=fig.add_axes([0.25,0.95,0.5 ,0.02])
         labels=fig.colorbar(cs,cax=cb_ax,ticks=[-1,0,1],orientation='horizontal',extend='neither',extendfrac=0.00001,extendrect=True,drawedges=False)
         labels.ax.set_xticklabels(['-','0','+',''])
         pdf.savefig()
         plt.close()

pdf.close()
print "No worries, mate..."
