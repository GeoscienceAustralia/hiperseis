import xcorqc
import argparse
import datetime
import glob, os
from os.path import join, exists
import json

# =========================== User Input Required =========================== #

#Path to the data
data_path = '/g/data/ha3/Passive/'

#IRIS Virtual Ntework name
virt_net = '_ANU'

# FDSN network identifier (2 Characters)
FDSNnetwork = '7G(2013-2015)'

# Component to analyse
comp = 'BHZ'

network_path = join(data_path, virt_net, FDSNnetwork) 

years=('2014','2015','2016')

refstn_code = 'INKA'

pltoutdir = '/g/data/ha3/Passive/xcor4/'

arraydatadir = '/g/data/ha3/Passive/_ANU/7G(2013-2015)/ASDF/'
arrayjsonfile = '7G(2013-2015)_raw_dataDB.json'
arrayasdffile = '7G(2013-2015).h5'
refdatadir='/g/data/ha3/Passive/Ref/'
refjsonfile='INKA.json' #FIXME put the ref station here
refasdffile = 'INKA.h5'

# =========================================================================== #
# Load data

array_sdb = SeisDB(arraydatadir+arrayjsonfile)
array_asdf = pyasdf.ASDFDataSet(arraydatadir+arrayasdffile)

# Create the ref station data NOTE: this only needs to be done once
ref_sdb = SeisDB(refjsonfile + refjsonfile)
ref_asdf = pyasdf.ASDFDataSet(refdatadir+refasdffile)

# =========================================================================== #

st = Stream()
refst = Stream()

#ASHBY: Can you please write the ASDF waveform loading code from the above files
#       and funnel the streams into arrayst and refst respectively?
refst.merge(method=1,fill_value=0)
refst.print_gaps()
print "Merging deployment streams"
st.merge(method=1,fill_value=0)
st.print_gaps()
print "Running xcor"
#plotname = stn_+'_'+('%04d%02d'%(year,month))+"-"+('%03d%03d'%(doys[0],doys[len(doys)-1]))
[xcl,xcxl,complist] = IntervalStackXCorr(refst,st)
saveXCorr(xcl,xcxl,xcorrout,figname) # Saves numpy array
saveXCorrPlot(xcl,xcxl,plotout,figname) # Saves plot as png
