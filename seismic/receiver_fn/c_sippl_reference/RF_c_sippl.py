#!/usr/bin/env python
# coding=utf-8
# module RF
"""
Collection of functions for RF processing, data preparation and use of different routines.

This code used for this publication by Christian Sippl:
  https://www.sciencedirect.com/science/article/pii/S0040195116300245
"""

# pylint: skip-file

import os
import sys
import subprocess

if sys.version_info[0] == 2:
  import cPickle as pkl
else:
  import pickle as pkl

import obspy
import iris
import basic_c_sippl
from util_c_sippl import KM_PER_DEG
from obspy.geodetics import gps2dist_azimuth
from pylab import *
from numpy import *
from glob import glob
# from obspy.taup import taup  # broken, doesn't exist
from obspy.core import read, Stream, UTCDateTime, Trace
# This is based on an *ancient* version of obspy (0.10.2), not compatible with obspy 1.x.x.
# from obspy.sac import sacio
from obspy.taup import TauPyModel
#from toeplitz import sto_sl

def stack(inlist,outfile):
  """
  stack seismogram traces (before test if lengths are identical)
  in: Streams (either one-component or more), given in inlist

  TODO: response equalization before??
  """
  # list assumed to be e.g. glob output (paths)
  e = obspy.core.Stream() #pre-define empty component Streams
  n = obspy.core.Stream()
  z = obspy.core.Stream()
  for j in inlist:
    trace_new = read(j)
    if trace_new[0].stats.channel[-1] == 'E':
      e += trace_new
    elif trace_new[0].stats.channel[-1] == 'N':
      n += trace_new
    elif trace_new[0].stats.channel[-1] == 'Z':
      z += trace_new

  # check if all Streams have content..
  if len(e) != 0:
    #stack E
    e.detrend(type='demean')
    e_stack = e[0].data
    len_e = len(e[0].data)
    for trace in e[1:]:
      if len(trace.data) == len_e:
        e_stack += trace.data
      else:
        print('Trace lengths do not fit - trace not added!')

  if len(n) != 0:
    #stack N
    n.detrend(type='demean')
    n_stack = n[0].data
    len_n = len(n[0].data)
    for trace in n[1:]:
      if len(trace.data) == len_n:
        n_stack += trace.data
      else:
        print('Trace lengths do not fit - trace not added!')

  if len(z) != 0:
    #stack Z
    z.detrend(type='demean')
    z_stack = z[0].data
    len_z = len(z[0].data)
    for trace in z[1:]:
      if len(trace) == len_z:
        z_stack += trace.data
      else:
        print('Trace lengths do not fit - trace not added!')

  return e_stack, n_stack, z_stack         

def prepare(eventfile,station,datapath,sampfreq,filt=False,rotate='2D',info_file='/home/sippl/info_file',local=False):
  """
  prepare data for RF deconvolution: returns dictionary containing event information and R,T,Z traces, downsampled to sampfreq

  """

  #read in eventinfo file, get the data
  evfile = open(eventfile,'r')
  data = evfile.readlines()
  evfile.close()

  sta = basic_c_sippl.Station(station,info_file)
  stat_lat = sta.lat
  stat_lon = sta.lon
  stat_elev = sta.elev

  #how to handle station coordinates?? -- infofile (get template)
  model = TauPyModel(model='ak135')

  bigdict = {} 
  readflag = True
  for j in data:
    eid, yr, mn, dy, time, lat, lon, dep, typ, ttime, tof = j.strip('\n').split(None)
    
    bigdict[int(eid)] = {}
    bigdict[int(eid)]['station'] = station
    bigdict[int(eid)]['event_lat'] = float(lat)
    bigdict[int(eid)]['event_lon'] = float(lon)
    bigdict[int(eid)]['event_depth'] = float(dep)
    bigdict[int(eid)]['station_lat'] = stat_lat
    bigdict[int(eid)]['station_lon'] = stat_lon
    bigdict[int(eid)]['station_elev'] = stat_elev
    bigdict[int(eid)]['Traces'] = {}
    bigdict[int(eid)]['orig_time'] = UTCDateTime(int(yr),int(mn),int(dy),int(time.split(':')[0]),int(time.split(':')[1]),float(time.split(':')[2]))

    #get data files
    print(eid)

    try:
      ev_folder = glob(datapath+'/'+eid+'__*')[0]
      str_orig = read(ev_folder+'/*H?*',format='MSEED')
      e = str_orig.select(component='E')
      n = str_orig.select(component='N')
      z = str_orig.select(component='Z')
    except IndexError:
      continue

    try:
      assert len(e[0]) == len(n[0]) == len(z[0])
    except AssertionError:
      print("Different trace lengths: not processed")
      del bigdict[int(eid)]
      continue

    #check that none of the components is flat
    if max(e[0].data) - min(e[0].data) < 5:
      print("Flat component: not processed")
      del bigdict[int(eid)]
      eid = str(int(eid) - 1)
      continue
    if max(n[0].data) - min(n[0].data) < 5:
      print("Flat component: not processed")
      del bigdict[int(eid)]
      eid = str(int(eid) - 1)
      continue

    if max(z[0].data) - min(z[0].data) < 5:
      print("Flat component: not processed")
      del bigdict[int(eid)]
      eid = str(int(eid) - 1)
      continue          

   
#    try:
#      assert e[0].data.all() != n[0].data.all()
#      assert e[0].data.all() != z[0].data.all()
#      assert n[0].data.all() != z[0].data.all()   
#    except AssertionError:
#      print("Two or more components are identical: not processed")
#      continue

    distm,az,baz = gps2dist_azimuth(float(lat),float(lon),stat_lat,stat_lon)

    dist = distm/(1e3 * KM_PER_DEG)

    bigdict[int(eid)]['Distance'] = round(dist,2)
    bigdict[int(eid)]['Azimuth'] = round(az,2)
    bigdict[int(eid)]['Backazimuth'] = round(baz,2)

    #get ray parameter

    if local:
      raytr = iris.ttime_dict(bigdict[int(eid)]['event_lat'], bigdict[int(eid)]['event_lon'],
                              bigdict[int(eid)]['event_depth'], bigdict[int(eid)]['station_lat'],
                              bigdict[int(eid)]['station_lon'], phase='p')
    else:
      raytr = iris.ttime_dict(bigdict[int(eid)]['event_lat'], bigdict[int(eid)]['event_lon'],
                              bigdict[int(eid)]['event_depth'], bigdict[int(eid)]['station_lat'],
                              bigdict[int(eid)]['station_lon'], phase='P')
    rayp = round(raytr['rayp'],5)
    print('RAYP:')
    print(rayp)
    bigdict[int(eid)]['ray parameter'] = rayp

    #calculate Moho piercing point
    if not local:
      arr = model.get_pierce_points(float(dep),round(dist,2),phase_list='P')
    else:
      arr = model.get_pierce_points(float(dep),round(dist,2),phase_list='p')

    pierce_dist = (arr[0].pierce[-1][2] - arr[0].pierce[-3][2])*180./pi*KM_PER_DEG
    print('PIERCE DIST:')
    print(pierce_dist)
    bigdict[int(eid)]['Pierce_distance_P'] = round(pierce_dist,2)

    #compute rotation
    if rotate == '2D':
      str_rot = str_orig.rotate('NE->RT',back_azimuth=baz)
      r = str_rot.select(component='R')[0]
      t = str_rot.select(component='T')[0]      

    elif rotate == '3D':
      arr2 = model.get_travel_times(bigdict[int(eid)]['event_depth'],bigdict[int(eid)]['Distance'],phase_list=['P'])
      inc = arr2[0].incident_angle
      str_rot = str_orig.rotate('ZNE->LQT',back_azimuth=baz,inclination=inc)
      r = str_rot.select(component='Q')[0]   
      t = str_rot.select(component='T')[0]
      z = str_rot.select(component='L')[0]   

    #downsampling
    factor = (z[0].stats.sampling_rate / sampfreq)
    if not factor%1 == 0:
      print('Error: non-integer decimation factor!!')
      del bigdict[int(eid)]
      continue
    else:
      z.decimate(5)
      r.decimate(5)
      t.decimate(5)
      z.decimate(5)
      r.decimate(5)
      t.decimate(5)
      z[0].detrend(type='demean')
      r.detrend(type='demean')
      t.detrend(type='demean')

      if filt:
        z.filter('bandpass',freqmin=0.1,freqmax=5.)
        r.filter('bandpass',freqmin=0.1,freqmax=5.)
        t.filter('bandpass',freqmin=0.1,freqmax=5.)
#
      try:
        assert z[0] != zeros(len(z[0]))
        assert r != zeros(len(r))
        assert t != zeros(len(t))
      except AssertionError:
        print("At least one traces is only zeros: not processed")
        continue

      z[0].detrend(type='linear')
      r.detrend(type='linear')
      t.detrend(type='linear')

    if rotate == '2D':
      #try:
      #  z[0].slice(z[0].stats.starttime+55.,z[0].stats.endtime-60.)
      #  z[0].taper(max_percentage=0.5,type='hann')
      #  r.slice(r.stats.starttime+55.,r.stats.endtime-60.)
      #  r.taper(max_percentage=0.5,type='hann')
      #  t.slice(t.stats.starttime+55.,t.stats.endtime-60.)
      #  t.taper(max_percentage=0.5,type='hann')
      #except ValueError: #data missing
      #  del bigdict[int(eid)]
      #  continue 
      bigdict[int(eid)]['Traces']['Z'] = z[0]
      bigdict[int(eid)]['Traces']['R'] = r
      bigdict[int(eid)]['Traces']['T'] = t
    elif rotate == '3D':
      bigdict[int(eid)]['Traces']['L'] = z[0]
      bigdict[int(eid)]['Traces']['L'].stats.channel = 'BHL'
      bigdict[int(eid)]['Traces']['Q'] = r
      bigdict[int(eid)]['Traces']['Q'].stats.channel = 'BHQ'
      bigdict[int(eid)]['Traces']['T'] = t
      bigdict[int(eid)]['Traces']['T'].stats.channel = 'BHT'

  return bigdict


def deconv_3D(indict,source_component='L',method='time'):
  """
  deconvolution method(s) for 3D deconvolution; modified from Tom Richter's code (https://github.com/trichter/rf)
  """
  winsrc = (-10,30,5)
  winrsp = (-20,80)
  winrf = (-20,80)

  #somehow get that into processing stream...
  #make stream of traces
  for j in indict.keys():
    #create stream of the 3 components
    st = Stream()
    try:
      st += indict[j]['Traces']['L'].copy()
      st += indict[j]['Traces']['Q'].copy()
      st += indict[j]['Traces']['T'].copy() 
    except KeyError:
      del indict[j]
      continue

    samp = st[0].stats.sampling_rate
    onset = st[0].stats.starttime + 60.
 
    src = st.select(component=source_component)[0]
    src_index = st.traces.index(src)
    src = src.copy()
    src.trim(onset + winsrc[0], onset + winsrc[1])
    src.taper(max_percentage=None, max_length=winsrc[2])
    src = src.data
    rsp = [st[(i + src_index) % 3].data for i in range(3)]
    
    if method == 'time':
      time_rf = winrf[1] - winrf[0]
      shift = int(round(samp * ((winrsp[1] - winrsp[0] - winsrc[1] + winsrc[0] - time_rf) / 2 + winrsp[0] - winsrc[0] - winrf[0])))
      length = int(round(time_rf * samp)) + 1
      rf_resp = deconvt(rsp, src, shift, length=length)
      tshift = -winrf[0]
      for tr in st:
        tr.stats.tshift = tshift  
    else:
      rf_resp = deconvf(rsp, src, samp, **kwargs)
    for i in range(3):
      st[(i + src_index) % 3].data = rf_resp[i].real
    st_q = st.select(component='Q')[0]
    st_t = st.select(component='T')[0]
    q = st.traces.index(st_q)
    t = st.traces.index(st_t)
    #cut to same parameters as with 2D??
    st[q].filter('bandpass',freqmin=0.1,freqmax=1.0,corners=4)
    st[t].filter('bandpass',freqmin=0.1,freqmax=1.0,corners=4)
    indict[j]['Traces']['QRF'] = st[q].slice((onset - 5.),(onset + 80.))
    indict[j]['Traces']['TRF'] = st[t].slice((onset - 5.),(onset + 80.))

  i = 1
  indict_new = {}
  for j in sorted(indict.keys()):
    indict_new[i] = indict[j]
    i += 1
  
  return indict_new

def CCP_master(startpoint, endpoint, width, spacing, depth, v_background='ak135', info_file='/home/sippl/info_file',
               rfdicpath='/home/sippl/sandbox/RF/ALFREX/dicts_a2.5'):
  """
  fully automatic plotting of CCP stacks, selection of profile parameters (starting point, azimuth, length, width) in inout, stations selected based on width
  point in (lat,lon), decimal degrees
  az in degrees from north
  length, width, spacing, depth in km
  
  plotting parameters?
  """

  #check which stations are in --> compute profile range plus width range
  #get angle and distance
  ybig = max(startpoint[0],endpoint[0])
  ysmall = min(startpoint[0],endpoint[0])
  xbig = max(startpoint[1],endpoint[1])
  xsmall = min(startpoint[1],endpoint[1])
  dx = (xbig - xsmall) * KM_PER_DEG * cos((ybig+ysmall)/2. * pi / 180.)
  dy = (ybig - ysmall) * KM_PER_DEG
  try:
    az = (arctan(dy/(float(dx)))*180.)/pi  # Why not using atan2 here?
  except ZeroDivisionError:
    az = 0.
  length = sqrt(dx**2 + dy**2)

  #find quadrant (additive term to angle...)
  if startpoint[0] >= endpoint[0] and startpoint[1] < endpoint[1]:
    add = 90.
  elif startpoint[0] < endpoint[0] and startpoint[1] <= endpoint[1]:
    add = 0.
  elif startpoint[0] > endpoint[0] and startpoint[1] >= endpoint[1]:
    add = 180.
  elif startpoint[0] <= endpoint[0] and startpoint[1] > endpoint[1]:
    add = 270.

  az += add

  #set velocity model (other options can be added) 
  if v_background == 'ak135':
    b_lay = [0.,20.,35.,77.5,120.,165.,210.,260.,310.,360.,410.,460.,510.,560.,610.,660.]
    vp = [5.8,6.5,8.04,8.045,8.05,8.175,8.3,8.4825,8.665,8.8475,9.03,9.36,9.528,9.696,9.864,10.032,10.2]
    vs = [3.46,3.85,4.48,4.49,4.5,4.509,4.518,4.609,4.696,4.783,4.87,5.08,5.186,5.292,5.398,5.504,5.61]

    vmod = (b_lay,vp,vs)

  #first set up a matrix
  matrx,depstep,lenstep = setup_ccp_profile(length,spacing,depth)

  matrx_entry = matrx.copy()
  #go through info_file, get lat,lon, check here...
  info = open(info_file,'r')
  info_dat = info.readlines()
  info.close()

  #create map plot file
  from mpl_toolkits import basemap
  m = basemap.Basemap(projection='merc',urcrnrlat=ybig+1.,urcrnrlon=xbig+1.,llcrnrlon=xsmall-1.,llcrnrlat=ysmall-1.,resolution='h')
  m.drawcoastlines()
  x1,y1 = m(startpoint[1],startpoint[0])
  x2,y2 = m(endpoint[1],endpoint[0])

  m.plot([x1,x2],[y1,y2],'r--') 

  stationlist = []

  for line in info_dat:
    stat,dum,dum,dum,dum,dum,dum,lat,lon,dum,dum,dum,dum,dum,dum,dum,dum = line.split(None)
    sta_lat = float(lat)
    sta_lon = float(lon)

    #print(stat,sta_lat,sta_lon)

    angle_norm = az%90

    #calculate position on profile (length)
    sta_offset_n = ((startpoint[0] - sta_lat) * KM_PER_DEG) / sin(angle_norm*pi/180.)
    sta_offset_e = (-1) * ((startpoint[1] - sta_lon) * KM_PER_DEG * cos(startpoint[0] * pi /180.)) / sin((90.-angle_norm)*pi/180.)
    
    """
    if abs(az-90.) > 30. and abs(az-270.) > 30.:
      sta_offset = (-1) * sta_offset_n / cos(az*pi/180.)
    else: #calculate the other way...
      sta_offset = (-1) * sta_offset_e / cos((90. - az)*pi/180.)

    """
   
    xstart = 0
    ystart = 0
    xend = (endpoint[1] - startpoint[1]) * KM_PER_DEG * cos((endpoint[1] + startpoint[1])/2. * pi / 180.)
    yend = (endpoint[0] - startpoint[0]) * KM_PER_DEG
    xstat = (sta_lon - startpoint[1]) * KM_PER_DEG * cos((sta_lon + startpoint[1])/2. * pi / 180.)
    ystat = (sta_lat - startpoint[0]) * KM_PER_DEG

    dist = ((yend - ystart) * xstat - (xend - xstart) * ystat + xend*ystart - yend*xstart) / sqrt((yend - ystart)**2 + (xend - xstart)**2)

    if abs(dist) <= width:

      if abs(az-90.) > 30. and abs(az-270.) > 30.:
        n_correction = dist / tan(angle_norm*pi/180.)
        sta_offset = sta_offset_n + n_correction

      else:
        e_correction = dist / tan((90 - angle_norm)*pi/180.)
        sta_offset = sta_offset_e - e_correction

      if sta_offset > 0 and sta_offset < length:
        print("Station "+stat+" included!!")
        print(sta_offset, dist)
        #add station to CCP stack
        x,y = m(sta_lon,sta_lat)
        m.plot(x,y,'ro')
        #stationlist.append(stat)
        
        try:
          indict = pkl.load(open(glob(rfdicpath+'/'+stat+'_[Rr][Ff]*')[0],'rb'))
          for i in indict.keys():
            azi = indict[i]['Backazimuth']
            model = TauPyModel(model='ak135')
            arr_p = model.get_travel_times(indict[i]['event_depth'],indict[i]['Distance'],phase_list=['P'])
            inc_p = arr_p[0].incident_angle
            matrx,matrx_entry = add_ccp_trace(indict[i]['Traces']['RF'],inc_p,matrx,matrx_entry,vmod,depstep,lenstep,sta_offset,azi)
  
        except IndexError: #no dictionary for this station...
          pass
        
      else:
        x,y = m(sta_lon,sta_lat)
        m.plot(x,y,'k^')

    else:
      x,y = m(sta_lon,sta_lat)
      m.plot(x,y,'k^')
    # end if abs(dist) <= width

  #normalize, then plot
  matrx_norm = (matrx / matrx_entry).transpose()
  plot_ccp(matrx_norm,length,depth,spacing,ofile='test.pdf')


def setup_ccp_profile(length,spacing,maxdep):
  """
  construct the grid for a CCP stacking profile
  """
  #calculate number of cells in x and y direction 
  n_y = int(round((maxdep/spacing),0))
  n_x = int(round(length/spacing,0)) 

  #get center values
  depstep = arange(0 + round(spacing/2.,1), maxdep,spacing)
  lenstep = arange(0 + round(spacing/2.,1), length,spacing)

  #create matrix
  mtrx = zeros([n_x,n_y])

  return mtrx,depstep,lenstep

def add_ccp_trace(trace,inc_p,matrx,matrx_entry,vmod,depstep,lenstep,sta_offset,az):
  """
  project amplitudes from all RFs in indict onto the profile...2D rot: 
  """
  #start at zero: inc_p given, inc_s needs to be calculated  
  h = 0
  c = 0
  d = 0
  tpz = 0
  tsz = 0
  rpz = 0
  rsz = 0

  for j in range(1,len(depstep)):
    c = d
    if j == 1:
      h = h_tot = 1
    else:
      h = depstep[j] - depstep[j-1]
      h_tot += h
    #check in velocity model
    for f in range(len(vmod[0])):
      if vmod[0][f] < depstep[j]:
        d = f      
        
    #derive P incidence from previous P incidence, then current S from current P
    inc_p = arcsin((sin(inc_p *pi/180.) * vmod[1][d]) / vmod[1][c]) * 180 / pi  
    inc_s = arcsin((sin(inc_p *pi/180.) * vmod[2][d]) / vmod[1][d]) * 180 / pi
      
    #horizontal distances (attention: still needs azimuth normalization)
    rpz += h*tan(inc_p/180.*pi) 
    rsz += h*tan(inc_s/180.*pi)
   
    rd = rpz - rsz
 
    tpz += h/cos(inc_p/180.*pi) * (1/vmod[1][d])
    tsz += h/cos(inc_s/180.*pi) * (1/vmod[2][d]) 
    
    td = sin(inc_p/180.*pi)/vmod[1][d] * rd
    #check if velocity jump, if yes get new angles
    tps = tsz + td - tpz

    amp = get_amplitude(trace,tps)
    
    #project, put into correct bin in matrix
    xsz = rsz * cos(az*pi/180.) #relative to station
    indx_x,indx_y = matrx_lookup(xsz,sta_offset,h_tot,depstep,lenstep)
    
    matrx[indx_x,indx_y] += amp
    matrx_entry[indx_x,indx_y] += 1

  return matrx,matrx_entry


def matrx_lookup(xsz,sta_offset,h,depstep,lenstep):
  """
  return index values for amplitude contrbution in profile matrix
  """
  distance_offset = sta_offset - xsz #because zero is in the north

  diff_x = 999.
  diff_y = 999. 
  indx_x = 0
  indx_y = 0

  #find horizontal position
  for j in range(len(lenstep)):
    if abs(lenstep[j] - distance_offset) < diff_x:
      diff_x = abs(lenstep[j] - distance_offset)
      indx_x = j

  for k in range(len(depstep)):
    if abs(depstep[k] - h) < diff_y:
      diff_y = abs(depstep[k] - h)
      indx_y = k

  return indx_x, indx_y


def get_amplitude(trace,time,rf_offset=5.):
  """
  retrieve amplitude value
  """
  indx = (time + rf_offset) * trace.stats.sampling_rate
  amp = trace.data[indx]
  return amp 


def plot_ccp(matrx,length,maxdep,spacing,ofile='bla.pdf'):
  """
  plot results of CCP stacking procedure
  """
  tickstep_x = 50
  tickstep_y = 25

  figure()
  imshow(matrx,aspect='equal',vmin=-0.15,vmax=0.15)

  ylim([int(maxdep/spacing),0]) 
 
  xlabel('distance [km]')
  ylabel('depth [km]') 

  xticks(range(0,int(length/spacing),int(tickstep_x/spacing)),range(0,int(length),tickstep_x))
  yticks(range(0,int(maxdep/spacing),int(tickstep_y/spacing)),range(0,int(maxdep),tickstep_y))

  savefig(ofile)

  #amplitude weighting
  #right scaling
  

def rf_viewer(indict,time_win=[0,35]):
  """
  shows plot of RF, interactive choice whether the trace should be kept or discarded
  """
  from builtins import input
  dict_new = {}
  count = 1
  for j in sorted(indict.keys()):
    clf()
    xvals = (arange(len(indict[j]['Traces']['RF'])))/indict[j]['Traces']['RF'].stats.sampling_rate
    plot(xvals,indict[j]['Traces']['RF'],'k-')
    xlim(time_win)
    print("Input needed: choose (r) to reject the trace or (k) to keep it!")
    key = input()

    if key == 'r':
      continue
    elif key == 'k':
      dict_new[count] = indict[j].copy()
      count += 1
  return dict_new

def deconvt(rsp_list, src, shift, spiking=1., length=None, normalize=True):
    """
    Time domain deconvolution.
    Deconvolve src from arrays in rsp_list.
    Calculate Toeplitz auto-correlation matrix of source, invert it, add noise
    and multiply it with cross-correlation vector of response and source.
    In one formula: ::
        RF = (STS + spiking*I)^-1 * STR
        N... length
            ( S0   S-1  S-2 ... S-N+1 )
            ( S1   S0   S-1 ... S-N+2 )
        S = ( S2   ...                )
            ( ...                     )
            ( SN-1 ...          S0    )
        R = (R0 R1 ... RN-1)^T
        RF = (RF0 RF1 ... RFN-1)^T
        S... source matrix (shape N*N)
        R... response vector (length N)
        RF... receiver function (deconvolution) vector (length N)
        STS = S^T*S = symetric Toeplitz autocorrelation matrix
        STR = S^T*R = cross-correlation vector
        I... Identity
    :param rsp_list: either a list of arrays containing the response functions
        or a single array
    :param src: array of source function
    :param shift: shift the source by that amount of samples to the left side
        to get onset in RF at the desired time (negative -> shift source to the
        right side)\n
        shift = (middle of rsp window - middle of src window) +
        (0 - middle rf window)
    :param spiking: random noise added to autocorrelation (eg. 1.0, 0.1)
    :param length: number of data points in results
    :param normalize: normalize all results so that the maximum of the first
        result is 1
    :return: (list of) array(s) with deconvolution(s)
    """
    if length is None:
        length = len(src)
    flag = False
    RF_list = []
    STS = _acorrt(src, length)
    STS = STS / STS[0]
    STS[0] += spiking
    if not isinstance(rsp_list, (list, tuple)):
        flag = True
        rsp_list = [rsp_list]
    for rsp in rsp_list:
        STR = _xcorrt(rsp, src, length // 2, shift)
        if len(STR) > len(STS):
            STR = np.delete(STR, -1)
        RF = _toeplitz_real_sym(STS, STR)
        RF_list.append(RF)
    if normalize:
        norm = 1 / np.max(np.abs(RF_list[0]))
        for RF in RF_list:
            RF *= norm
    if flag:
        return RF
    else:
        return RF_list

def deconvf(rsp_list, src, sampling_rate, water=0.05, gauss=2., tshift=10.,
            pad=0, length=None, normalize=True, normalize_to_src=False,
            return_dict=False):
    """
    Frequency-domain deconvolution using waterlevel method.
    Deconvolve src from arrays in rsp_list.
    :param rsp_list: either a list of arrays containing the response functions
        or a single array
    :param src: array of source function
    :param sampling_rate: sampling rate of the data
    :param water: waterlevel to stabilize the deconvolution
    :param gauss: Gauss parameter of Low-pass filter
    :param tshift: delay time 0s will be at time tshift afterwards
    :param pad: multiply number of samples used for fft by 2**pad
    :param length: number of data points in results, optional
    :param normalize: if results are normalized
    :param normalize_to_src: True ->  normalized so that the maximum of a
        deconvolution of the source with itself is 1\n
        False -> normalized so that the maximum of the deconvolution of the
        first response array in rsp_list is 1
    :param return_dict: return additionally a lot of different parameters in a
        dict for debugging purposes
    :return: (list of) array(s) with deconvolution(s)
    """
    if length is None:
        length = len(src)
    N = length
    nfft = nextpow2(N) * 2 ** pad
    freq = np.fft.fftfreq(nfft, d=1. / sampling_rate)
    gauss = np.exp(np.maximum(-(0.5 * 2 * pi * freq / gauss) ** 2, -700.) -
                   1j * tshift * 2 * pi * freq)

    spec_src = fft(src, nfft)
    spec_src_conj = np.conjugate(spec_src)
    spec_src_water = np.abs(spec_src * spec_src_conj)
    spec_src_water = np.maximum(spec_src_water, max(spec_src_water) * water)

    if normalize_to_src:
        spec_src = gauss * spec_src * spec_src_conj / spec_src_water
        rf_src = ifft(spec_src, nfft)[:N]
        norm = 1 / max(rf_src)
        rf_src = norm * rf_src

    flag = False
    if not isinstance(rsp_list, (list, tuple)):
        flag = True
        rsp_list = [rsp_list]
    rf_list = [ifft(gauss * fft(rsp, nfft) * spec_src_conj / spec_src_water,
                    nfft)[:N] for rsp in rsp_list]
    if normalize:
        if not normalize_to_src:
            norm = 1. / max(rf_list[0])
        for rf in rf_list:
            rf *= norm
    if return_dict:
        if not normalize_to_src:
            spec_src = gauss * spec_src * spec_src_conj / spec_src_water
            rf_src = ifft(spec_src, nfft)[:N]
            norm = 1 / max(rf_src)
            rf_src = norm * rf_src
        ret_dict = {'rf_src': rf_src, 'rf_src_conj': spec_src_conj,
                    'spec_src_water': spec_src_water, 'freq': freq,
                    'gauss': gauss, 'norm': norm, 'N': N, 'nfft': nfft}
        return rf_list, ret_dict
    elif flag:
        return rf
    else:
        return rf_list

def deconv(bigdict,nit,psh=5,minerr=0.01,a=2.5,mode='R',wdir='/home/sippl/programs/Rftn.Codes.Export/IterDecon'):
  """
  wrap of C. Ammon's iterdecon code
  I/O in obspy's data format...
  Gaussian filter contained
  mode: choose whether R or T receiver functions are to be calculated
  """

  os.chdir(wdir)

  for j in bigdict.keys():

    #produce SAC files from obspy traces
    try:
      bigdict[j]['Traces'][mode].write('numerator','SAC')
      bigdict[j]['Traces']['Z'].write('denominator','SAC')
    except KeyError:
      continue
    #setup input temp file

    xx = open('xxx','w')
    xx.write('numerator\ndenominator\n'+str(nit)+'\n'+str(psh)+'\n'+str(minerr)+'\n'+str(a)+'\n1\n0')
    xx.close()

    #call code with input params
    process = subprocess.Popen("iterdecon < xxx", shell=True,stdout=subprocess.PIPE)
    out = process.communicate()[0]
    fit =  out.split('\n')[-3].split(None)[-4].strip('%')
    
    bigdict[j]['fit'] = fit

    os.system("rm xxx")

    #read in output
    rf = read('decon.out')
    observed = read('observed')
    predicted = read('predicted')

    bigdict[j]['Traces']['RF'] = rf[0]
    bigdict[j]['Traces']['Predicted'] = predicted[0]
    bigdict[j]['Traces']['Observed'] = observed[0]

  return bigdict
  #evaluate: overlay plot obs + pre, Xcorr


#plotting routines!!

def get_az_plot(indict,minfit=85.,azbin=45): #stacked RF every azbin degrees (if available)...quality criteria applied before
  """
  fill out negative wiggles!!
  """
  nbins = 360 / azbin

  typlen = len(indict[1]['Traces']['RF'])
  binar = []
  poslist = []

  #fill with empty Traces
  for k in range(nbins):
    binar.append(Trace())
    
  for trace in binar:
    trace.data = zeros(typlen) 
 
  #missing: if statements for limiting/quality control 
  for j in indict.keys():
    pos = int(indict[j]['Backazimuth'] / azbin)
    if pos == nbins:
      pos -= 1
    if not isnan(indict[j]['Traces']['RF'].data[0]):
      if float(indict[j]['fit']) > minfit:
        binar[pos].data += indict[j]['Traces']['RF'].data 
        poslist.append(pos)
  #amplitude normalization 
  figure(figsize=(12,12))
  sampr = indict[1]['Traces']['RF'].stats.sampling_rate
  for q in range(len(binar)):
    binar[q].data /= (max(max(binar[q].data),abs(min(binar[q].data))))/(2*azbin)
    x = binar[q].data + azbin*q + (0.5*azbin)
    plot(x)
    numb = poslist.count(q)
    text(45*sampr,(azbin*q + (0.5*azbin)),numb)
  xlim([0,50*sampr])
  ylabel('Backazimuth')
  xlabel('time [s]')
  xticks([0,5*sampr,10*sampr,15*sampr,20*sampr,25*sampr,30*sampr,35*sampr,40*sampr,45*sampr,50*sampr],
         ['0','5','10','15','20','25','30','35','40','45','50'])

  show()

  savefig('azplot.pdf')


def get_dist_plot(indict,minfit=85.,distbin=5): #same thing for distance bins
  """
  distance in deg!! 
  """
  #minimum should be 30 degrees, maximum 95 
  nbins = 70 / distbin
  typlen = len(indict[1]['Traces']['RF'])

  binar = []
  poslist = []

  #fill with empty Traces
  for k in range(nbins):
    binar.append(Trace())

  for trace in binar:
    trace.data = zeros(typlen)


  #missing: if statements for limiting/quality control 
  for j in indict.keys():
    pos = int((float(indict[j]['Distance']) - 30.) / float(distbin))
    if not isnan(indict[j]['Traces']['RF'].data[0]):
      if float(indict[j]['fit']) > minfit:
        binar[pos].data += indict[j]['Traces']['RF'].data
        poslist.append(pos)
  #amplitude normalization 
  figure(figsize=(12,12))
  sampr = indict[1]['Traces']['RF'].stats.sampling_rate
  for q in range(len(binar)):
    binar[q].data /= (max(max(binar[q].data),abs(min(binar[q].data))))/(2*distbin)
    x = binar[q].data + distbin*q + (0.5*distbin) + 30
    plot(x)
    numb = poslist.count(q)
    text(45*sampr,(distbin*q + (0.5*distbin))+30,numb)
  xlim([0,50*sampr])
  ylabel('Distance (deg)')
  xlabel('time [s]')
  xticks([0,5*sampr,10*sampr,15*sampr,20*sampr,25*sampr,30*sampr,35*sampr,40*sampr,45*sampr,50*sampr],
         ['-5','0','5','10','15','20','25','30','35','40','45'])

  show()

  savefig('distplot.pdf')


def az_select(indict,azrange=[0.,360.]):
  """
  delete non-acceptable events (new numbering!!)
  """
  newdict = {}
  counter = 1
  for j in indict.keys():
    if indict[j].has_key('Backazimuth'):
      if float(indict[j]['Backazimuth']) > azrange[0] and float(indict[j]['Backazimuth']) < azrange[1]:
        newdict[counter] = indict[j]
        counter += 1

  return newdict

def dist_select(indict,distrange=[30.,95.]):
  """
  for distance
  """
  newdict = {}
  counter = 1
  for j in indict.keys():
    if indict[j].has_key('Distance'):
      if float(indict[j]['Distance']) > distrange[0] and float(indict[j]['Distance']) < distrange[1]:
        newdict[counter] = indict[j]
        counter += 1

  return newdict

def amp_select(indict):
  """
  throw out RFs with too big amplitudes (as they are ususally crap)
  """
  newdict = {}
  counter = 1
  for j in indict.keys():
    if indict[j]['Traces']['RF'].data.max() <= 2. and indict[j]['Traces']['RF'].data.min() >= -2.:
      newdict[counter] = indict[j]
      counter += 1
  return newdict

def fit_select(indict,minfit=85.):
  """
  same thing for fit parameter
  """
  newdict = {}
  counter = 1
  for j in indict.keys():
    if indict[j].has_key('fit'):
      if float(indict[j]['fit']) > minfit:
        newdict[counter] = indict[j]
        counter += 1

  return newdict


def cc_matrix(indict,Chi=0.8,Tau=0.35,write=False,outfolder='.',rotate='2D'):
  """
  CC matrix method for RF selection (see Tkalcic 2011)
  acceptable ranges of az, dist and fit as input parameters
  """
  #set up matrix
  Chimat = zeros([len(indict),len(indict)])
  if rotate == '2D':
    trmark = 'RF'
  elif rotate == '3D':
    trmark = 'QRF'
  for j in range(len(indict)):
    Chimat[j][j] = 1 #no need for any calculation here...
    if j == len(indict) - 1:
      break
    for k in range(j+1,len(indict)):
      #compute CC coeffs
      Chimat[j][k] = obspy.signal.cross_correlation.xcorr(indict[j+1]['Traces'][trmark].data,indict[k+1]['Traces'][trmark].data,0)[1]
      Chimat[k][j] = Chimat[j][k]

  #run through matrix, select good traces
  gooddict = {}
  num = 1
  for l in range(len(Chimat)):
    #sum all occurrences of Chmiat(l) > Chi
    count = 0
    for m in range(len(Chimat)):
      if Chimat[l][m] >= Chi:
        count += 1
    if count/float(len(Chimat)) >= Tau:
      #good trace
      gooddict[num] = indict[l+1]
      num += 1

  if write: #write out SAC files for all good traces; plus create file
    filename = gooddict[1]['station']+'__allRFs.list'
    outf = open(outfolder+'/'+filename,'w')
    for i in gooddict.keys():
      outname = '%04i%02i%02i__%02i%02i%02i__%4s__RF.SAC' % (gooddict[i]['orig_time'].year, gooddict[i]['orig_time'].month,
        gooddict[i]['orig_time'].day, gooddict[i]['orig_time'].hour, gooddict[i]['orig_time'].minute, 
        gooddict[i]['orig_time'].second,gooddict[i]['station'])
      #outname = gooddict[i]['name']+'__sel'
      gooddict[i]['Traces'][trmark].write(outfolder+'/'+outname,'SAC')
       
      sacio().SetHvalueInFile(outfolder + '/' + outname, 'dist', gooddict[i]['Distance'])
      sacio().SetHvalueInFile(outfolder + '/' + outname, 'az', gooddict[i]['Azimuth'])
      sacio().SetHvalueInFile(outfolder + '/' + outname, 'baz', gooddict[i]['Backazimuth'])
      sacio().SetHvalueInFile(outfolder + '/' + outname, 'evla', gooddict[i]['event_lat'])
      sacio().SetHvalueInFile(outfolder + '/' + outname, 'evlo', gooddict[i]['event_lon'])
      sacio().SetHvalueInFile(outfolder + '/' + outname, 'evdp', gooddict[i]['event_depth'])

      raytr = iris.ttime_dict(gooddict[i]['event_lat'], gooddict[i]['event_lon'], gooddict[i]['event_depth'],
                              gooddict[i]['station_lat'], gooddict[i]['station_lon'], phase='P')
      outf.write(outname+'  '+str(round(raytr['rayp'],5))+'  0.00  '+str(round(gooddict[i]['Backazimuth'],3))+'  '+str(round(raytr['dist'],3))+'\n')

      #set header entries
    outf.close()
  return gooddict

def plot_ccmat(befdict,outdict):
  """
  associated plotting method
  """
  
  #calculate average
  sum_tr = outdict[1]['Traces']['RF']
  for k in range(2,len(outdict)+1):
    sum_tr.data += outdict[k]['Traces']['RF'].data
  sum_tr.data /= float(len(outdict))


  figure(figsize=(14,7))

  for j in befdict.keys():
    plot(befdict[j]['Traces']['RF'],'k',linewidth=0.5)
  for q in outdict.keys():
    plot(outdict[q]['Traces']['RF'],'b',linewidth=1)
  plot(sum_tr,'r',linewidth=2)
  sampr = befdict[1]['Traces']['RF'].stats.sampling_rate
  xlim([0,35*sampr])
  xticks([0,5*sampr,10*sampr,15*sampr,20*sampr,25*sampr,30*sampr,35*sampr],['-5','0','5','10','15','20','25','30'])
  xlabel('time [sec]')


  savefig('CC_plot.pdf')
  return sum_tr


def create_modelfile(deplist,vplist,vslist,outpath='./model.mod'):  #OK, works...
  """
  deplist: expects lower layer boundary here (i.e. not starting with zero)!
  create multi-layer model file for (e.g.) input in respknt
  columns: number, vp, vs, rho, thickness, qp, qs, strk, dip (last 4: default 0), poisson
  needed as input for creation of synthetic RF
  """
  outf = open(outpath,'w')
  #length check
  if not len(deplist) == len(vplist) == len(vslist):
    print("Entered lists have different lengths! No model file created.")
    return
  headerline = '%3i Generic                         \n' % (len(deplist))
  outf.write(headerline)
  dep_init = 0
  for j in range(len(deplist)):
    if j == (len(deplist) - 1):
      dep_init = deplist[j]
    poisson = (((vplist[j]/vslist[j])**2) - 2)/(((vplist[j]/vslist[j])**2) - 1)*(1/2.)  
    rho = 0.77+0.32*(vplist[j]) #taken from icmod.f, C.Ammon
    line = '%3i %8.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n' % ((j+1),vplist[j],vslist[j],rho,deplist[j]-dep_init,0,0,0,0,poisson)
    outf.write(line)
    dep_init = deplist[j]
  outf.close() 


def output_bodin(RF,psh=5,cutoff=30.,outfile='RF_obs.dat',norm_fac=1.):
  """
  create ASCII file of a single receiver function, ready for input in T. Bodin's code
  """
  #shift so that 0 sec is at P arrival
  # cut off after 30 sec
  #format
  outf = open(outfile,'w')

  RF.slice(RF.stats.starttime,(RF.stats.starttime+cutoff+psh)) 
     
  for j in range(len(RF.data)-1):
    if j <= ((cutoff+psh)*RF.stats.sampling_rate):
      outf.write(' %9.6f    %.7g\n' % ((j/RF.stats.sampling_rate)-psh,(RF.data[j])/norm_fac))

  outf.close()


def bodin_eval(outdir,rf='RF_obs.dat',outpath='.',mode='RF'):
  """
  plotting utilities?
  Voronoi cells, Chain, single models,...
  Models, misfits,...  
  mode: either "RF" or "Joint"

  """
  rfin = open(rf,'r')
  rf_dat = rfin.readlines()
  rfin.close()

  if mode == 'RF':
    rf_synth = open('data_best.out','r')
  elif mode == 'Joint':
    rf_synth = open('data_bestrg.out','r')
  synth_dat = rf_synth.readlines()
  rf_synth.close()

  xar_obs = []
  xar_synth = []
  yar_obs = []
  yar_synth = []

  for j in range(len(rf_dat)):
    x,y = rf_dat[j].strip('\n').split(None)
    xar_obs.append(float(x))
    yar_obs.append(float(y))
    a,b = synth_dat[j].strip('\n').split(None)
    xar_synth.append(float(a))
    yar_synth.append(float(b))


  figure()
  plot(xar_obs,yar_obs,'b',label='observed')
  plot(xar_synth,yar_synth,'r',label='synthetic')
  legend()

  savefig('Curvefit.pdf')

  if mode == 'Joint':
    swdin = open(outdir+'/SWD_obs.dat','r')
    swd_dat = swdin.readlines()
    swdin.close()

    swdsy = open(outdir+'/data_bestdg.out','r')
    swd_synth = swdsy.readlines()
    swdsy.close()

    xar = []
    swd = []
    swsy = []
    for j in range(len(swd_dat)):
      xar.append(float(swd_dat[j].strip('\n').split(None)[0]))    
      swd.append(float(swd_dat[j].strip('\n').split(None)[1]))
      swsy.append(float(swd_synth[j].strip('\n')))

    clf()
    semilogx(xar,swd,'ko')
    semilogx(xar,swd,'k--')
    semilogx(xar,swsy,'ro')
    semilogx(xar,swsy,'r--')

    savefig('SWD_fit.pdf')


  #XXX major differences from here!
  #read in stuff

  if mode == 'RF':
    posterior = open(outdir+'/Posterior.out','r')
  elif mode == 'Joint':
    posterior = open(outdir+'/posteriorg.out','r')

  post_dat = posterior.readlines()
  prof, disd, d_max = post_dat[0].strip('\n').split(None)
  beta_min,beta_max,disv,width = post_dat[1].strip('\n').split(None)
  post_rest = post_dat[2:]

  if mode == 'RF':
    conv_misfit = open(outdir+'/Convergence_misfit.out','r')
    misf_dat = conv_misfit.readlines()
    c_misfit = zeros([2,len(misf_dat)])
    for y in range(len(misf_dat)):
      c_misfit[0][y], c_misfit[1][y] = misf_dat[y].strip('\n').split(None)

    conv_layers = open(outdir+'/Convergence_nb_layers.out','r')
    convlay_dat = conv_layers.readlines()
    c_layers = zeros([2,len(convlay_dat)])
    for a in range(len(convlay_dat)):
      c_layers[0][a], c_layers[1][a] = convlay_dat[a].strip('\n').split(None)


    conv_sigma = open(outdir+'/Convergence_sigma.out','r')
    csigma_dat = conv_sigma.readlines()
    csig = zeros([2,len(csigma_dat)])
    for b in range(len(csigma_dat)):
      csig[0][b], csig[1][b] = csigma_dat[b].strip('\n').split(None)


    conv_misfit.close()
    conv_layers.close()
    conv_sigma.close()

  if mode == 'RF':
    layers = open(outdir+'/NB_layers.out','r')
  elif mode == 'Joint':
    layers = open(outdir+'/Evidenceg.out','r')
  layers_dat = layers.readlines()

  if mode == 'RF':
    sigma = open(outdir+'/Sigma.out','r')
  elif mode == 'Joint':
    sigma = open(outdir+'/ML_Arg.out','r')
  sigmadat = sigma.readlines()
  sigmar = zeros([2,len(sigmadat)])
  for w in range(len(sigmadat)):
    sigmar[0][w], sigmar[1][w] = sigmadat[w].strip('\n').split(None)


  if mode == 'RF':
    average = open(outdir+'/Average.out','r')
  elif mode == 'Joint':
    average = open(outdir+'/Averageg.out','r')
  ave_data = average.readlines()
  
  ave = zeros([2,len(ave_data)])
  for q in range(len(ave_data)):
    ave[0][q], ave[1][q] = ave_data[q].strip('\n').split(None) 


  if mode == 'RF':
    cp = open(outdir+'/Change_points.out','r')
  elif mode == 'Joint':
    cp = open(outdir+'/CPg.out','r')

  cp_data = cp.readlines()
 
  cpar = zeros([2,len(cp_data)]) 
  for t in range(len(cp_data)):
    cpar[0][t], cpar[1][t] = cp_data[t].strip('\n').split(None)


  sigma.close()
  average.close()
  posterior.close()
  cp.close()
  layers.close()

  #Main plot

  count = 0
  P = zeros([int(disd),int(disv)])
  maxx = zeros([int(disd),1])
  for i in range(int(disd)):
    for j in range(int(disv)):
      P[i][j] = float((post_rest)[count].strip('\n')) 
      count += 1
    maxx[i] = P[i].argmax()


  figure(figsize=(14,9))
  x = [float(beta_min),float(beta_max)]
  y = [0,float(prof)]
 
  maxx = (maxx - 0.5) * ((float(beta_max) - float(beta_min))/float(disv)) + float(beta_min)

  v = 31
  n = (v*float(disd)/float(prof)) + 0.5

  
  subplot(131)
  imshow(P, cmap='jet', extent = [x[0],x[1],y[1],y[0]],aspect='auto')
  xlabel('Vs [km/s]')
  ylabel('Depth [km]')  
  xlim([float(beta_min),float(beta_max)])
  ylim([0,float(prof)])
  gca().invert_yaxis()

  subplot(132)

  plot(ave[1],ave[0],'r')
  xlim([float(beta_min),float(beta_max)]) 
  plot([float(beta_min),(float(beta_max)-2*float(width))],[0,float(d_max)],'k')
  plot([float(beta_min)+2*float(width),float(beta_max)],[0,float(d_max)],'k')
  plot(maxx,ave[0],'k')
  gca().invert_yaxis() 

  subplot(133)
  plot(cpar[1]/max(cpar[1]),cpar[0],'k')
  fill_betweenx(cpar[0],cpar[1]/max(cpar[1]),0,color='darkgrey')
  xlim([0,1])
  xlabel('P(transition)')
  gca().invert_yaxis()

  savefig(outpath+'/main.pdf')

  # Histograms

  figure()

  subplot(211)
  bar(left=array(range(len(layers_dat)))-0.5,height=(array(layers_dat,dtype='int'))/sum(array(layers_dat,dtype='float')),width=1)
  xlim([1,len(layers_dat)])
  xlabel('# layers')
  ylabel('P(# layers)')

  subplot(212)
  bar(left=sigmar[0],height=(array(sigmar[1],dtype='float')/sum(array(sigmar[1],dtype='float'))),width=(sigmar[0][1] - sigmar[0][0]))  
  xlim([sigmar[0][0],sigmar[0][-1]])
  xlabel('$\sigma$ for RF')
  ylabel('p($\sigma$)')

  savefig('Histograms.pdf')

  # More plots...

  figure(figsize=(13,9))
  
  subplot(311)

  plot(c_misfit[0][1:])
  plot(c_misfit[1][1:])
  #plot([c_misfit[0][0],0],[c_misfit[0][0],c_misfit[1][0]],'r')
  ylabel('Misfit') #unit = ??


  subplot(312)

  plot(c_layers[0][1:])
  plot(c_layers[1][1:])
  ylabel('# layers')


  subplot(313)
  plot(csig[0][1:]) 
  plot(csig[1][1:]) 
  xlabel('Iteration number')
  ylabel('$\sigma$')
 
  savefig('Convergence_plots.pdf')

def get_pierce(bigdict):
  """
  add piercing distances to dictionary
  """
  model = TauPyModel(model='ak135')
  for i in bigdict.keys():
    dep = bigdict[i]['event_depth']
    dist = bigdict[i]['Distance'] 
    
    arr = model.get_pierce_points(float(dep),round(dist,2),phase_list='P')
    pierce_dist = (arr[0].pierce[-1][2] - arr[0].pierce[-3][2])*180./pi*KM_PER_DEG
    bigdict[i]['Pierce_distance'] = round(pierce_dist,2)

  return bigdict


def RF_plot(indict,sort='baz'):
  """
  routine to plot extracted RFs sorted by backazimuth (sort='baz') or ray parameter (sort='rayp')
  """
  
  figure(figsize=(5,12))
  if sort == 'baz':
    for j in indict.keys():
      y1 = indict[j]['Backazimuth']
      y2 = ((indict[j]['Traces']['RF'].data)/indict[j]['Traces']['RF'].data.max())/0.03 + indict[j]['Backazimuth']
      fill_between(range(len(indict[j]['Traces']['RF'].data)),y1,y2,where=y2>y1,color='red')
      fill_between(range(len(indict[j]['Traces']['RF'].data)),y1,y2,where=y2<y1,color='blue')
      plot(((indict[j]['Traces']['RF'].data)/indict[j]['Traces']['RF'].data.max())/0.03 + indict[j]['Backazimuth'],'k')
    xlim([0,300])
    ylim([-10,380])
    ylabel('Backazimuth [degrees]')
    xlabel('time [s]')
    xticks([50,150,250],[0,10,20]) 
    savefig('RF_baz.pdf') 

  elif sort == 'rayp':
    for j in indict.keys():
      raytr = iris.ttime_dict(indict[j]['event_lat'], indict[j]['event_lon'], indict[j]['event_depth'],
                              indict[j]['station_lat'], indict[j]['station_lon'], phase='P')
      rayp = round(raytr['rayp'],5) 
      y1 = rayp 
      y2 = ((indict[j]['Traces']['RF'].data)/indict[j]['Traces']['RF'].data.max())/2. + rayp
      fill_between(range(len(indict[j]['Traces']['RF'].data)),y1,y2,where=y2>y1,color='red')
      fill_between(range(len(indict[j]['Traces']['RF'].data)),y1,y2,where=y2<y1,color='blue')
      plot(((indict[j]['Traces']['RF'].data)/indict[j]['Traces']['RF'].data.max())/2. + rayp,'k')
      
    xlim([0,300])
    ylim([4,9.5])
    ylabel('ray parameter')
    xlabel('time [s]')
    xticks([50,150,250],[0,10,20]) 
    savefig('RF_rayp.pdf')

def RF_matrix_plot(indict,sort='baz'):
  """
  plot RFs as matrix (amplitude normalization?), sorted (only sorted) by baz or ray parameter
  """
 
  from matplotlib import gridspec
 
  ny = len(indict[1]['Traces']['RF'])
  nx = len(indict)

  matrx = zeros([ny,nx])

  if sort == 'baz':
    #collect all baz values
    bazlst = []
    for j in indict.keys():
      bazlst.append(indict[j]['Backazimuth'])
    bazlst.sort()
    count = 0
    for i in bazlst:
      #find which trace this corresponds to   
      for k in indict.keys():
        if indict[k]['Backazimuth'] == i:
          matrx[:,count] = indict[k]['Traces']['RF'].data
          count += 1
          break
    figure(figsize=(6,8))
    gs = gridspec.GridSpec(2,1,height_ratios=[1,4])
    subplot(gs[0])
    plot(bazlst)
    xticks([])
    xlim([0,len(bazlst)])
    ylim([0,360])
    yticks([0,100,200,300])
    ylabel('baz')

    subplot(gs[1])
    imshow(matrx[0:375,:],vmax=0.3,vmin=-0.3,aspect='auto')

    xlabel('Receiver function number')
    ylabel('time [s]')
    yticks([0,125,250,375],['-5','0','5','10'])
    savefig('baz__.pdf')

  elif sort == 'rayp':
    #collect all rayp values
    rayplst = []
    for j in indict.keys():
      rayplst.append(indict[j]['ray parameter'])
    rayplst.sort()
    count = 0
    for i in rayplst:
      #find which trace this corresponds to   
      for k in indict.keys():
        if indict[k]['ray parameter'] == i:
          matrx[:,count] = indict[k]['Traces']['RF'].data
          count += 1
          break
    figure(figsize=(6,8))
    gs = gridspec.GridSpec(2,1,height_ratios=[1,4])
    subplot(gs[0])
    plot(rayplst)
    xticks([])
    xlim([0,len(rayplst)])
    ylim([4,9.5])
    yticks([5,7,9])
    ylabel('rayp')

    subplot(gs[1])
    imshow(matrx[0:350,:],vmax=0.3,vmin=-0.3,aspect='auto')

    xlabel('Receiver function number')
    ylabel('time [s]')
    yticks([0,50,100,150,200,250,300,350],['-5','0','5','10','15','20','25','30'])
    savefig('rayp__.pdf')

def autocorr_teleseis(indict,stat,outdir,use_SNR=False):
  """
  try autocorrelating the teleseisms from the RF dictionaries
  """
  #selection criteria?? SNR or some such...
  #stack them afterwards
  #cut some more...

  allcorr = zeros(1001)
  count = 0

  for ky in indict.keys():
    try:
      ztrace = indict[ky]['Traces']['Z']
    except KeyError: #Z traces does not exist
      continue
    #check data length
    if len(ztrace.data) == 1801:

      #noise window max amp
      range_noise = ztrace.data[300:350].max() - ztrace.data[300:350].min()

      #signal window max amp
      range_signal = ztrace.data[600:650].max() - ztrace.data[600:650].min()

      SNR = range_signal / range_noise
      #discriminate
      if use_SNR:
        print(SNR)
        if SNR > 5:
          print("used")
          corr = obspy.signal.cross_correlation.xcorr(ztrace.data, ztrace.data, 50*int(ztrace.stats.sampling_rate),full_xcorr=True)[2]
          allcorr += corr
          count += 1
      else:
        corr = obspy.signal.cross_correlation.xcorr(ztrace.data, ztrace.data, 50*int(ztrace.stats.sampling_rate),full_xcorr=True)[2]
        allcorr += corr
        count += 1
    else:
      print(len(ztrace.data))
      continue
    
  #normalize 
  allcorr /= float(count)
  #save as dict into folder
  outdict = {}
  outdict[stat+'-'+stat] = {}
  outdict[stat+'-'+stat][1] = allcorr
  pkl.dump(outdict,open(outdir+'/'+stat+'-'+stat,'wb'))  


def matrix_all_rayp(diclist,spacing=0.25,freq=5.):
  """
  plot all RFs in given dictionaries in binned stack matrix plot against ray parameter
  freq - lowest sample frequency that everything can be downsampled to
  """
  #collect all values of ray parameter, get min and max, then 
  #create shrunk dict, only ray parameter and RF
  dic_new = {}
  count = 1

  for dic in diclist:
    dicc = pkl.load(open(dic,'rb'))
    for k in dicc.keys():
      dic_new[count] = {}
      dic_new[count]['rayp'] = dicc[k]['ray parameter']
      dic_new[count]['RF'] = dicc[k]['Traces']['RF']
      count += 1

  rayplst = []
  for j in dic_new.keys():
    rayplst.append(dic_new[j]['rayp'])
    dec_fac = int(dic_new[j]['RF'].stats.sampling_rate/freq)
    if not dec_fac == 1:
      dic_new[j]['RF'].decimate(dec_fac)

  mn = array(rayplst).min()
  mx = array(rayplst).max() 

  bins = arange(mn,mx,spacing)
  #then create matrix; downsample everything to 5 Hz
  mtrx = zeros([len(bins),freq*110])

  for rf in dic_new.keys():
    rayp = dic_new[rf]['rayp']
    dist = 999.
    for bn in bins:
      if abs(rayp-bn) < dist:
        dist = abs(rayp-bn)
        bb = bn
      pos = list(bins).index(bb)  
      mtrx[pos,:] += dic_new[rf]['RF'].data[0:110*freq].transpose()

  for col in range(len(bins)):
    mtrx[col,:] /= mtrx[col,:].max()

  #works, but plotting parameters are mostly hardcoded (not nice)
  imshow(mtrx.transpose(),aspect='auto',vmax=0.035,vmin=-0.035)
  yticks([25,75,125,175,225,275,325,375,425,475,525],[0,10,20,30,40,50,60,70,80,90,100])
  ylabel('time [s]')
  xlabel('Slowness')
  xticks([0,10,20,30,40],[bins[0],bins[10],bins[20],bins[30],bins[40]])
  xlim([0,42])
  return mtrx

def _add_zeros(a, num, side='both'):
    """Add num zeros at side of array a"""
    return np.hstack([np.zeros(num)] * (side in ('both', 'left')) + [a] +
                     [np.zeros(num)] * (side in ('both', 'right')))


def _acorrt(a, num):
    """
    Not normalized auto-correlation of signal a.
    Sample 0 corresponds to zero lag time. Auto-correlation will consist of
    num samples. Correlation is performed in time domain by scipy.
    :param a: Data
    :param num: Number of returned data points
    :return: autocorrelation
    """
    return correlate(_add_zeros(a, num - 1, 'right'), a, 'valid')


def _xcorrt(a, b, num, zero_sample=0):
    """
    Not normalized cross-correlation of signals a and b.
    :param a,b: data
    :param num: The cross-correlation will consist of 2*num+1 samples.\n
        The sample with 0 lag time will be in the middle.
    :param zero_sample: Signals a and b are aligned around the middle of their
        signals.\n
        If zero_sample != 0 a will be shifted additionally to the left.
    :return: cross-correlation
    """
    if zero_sample != 0:
        a = _add_zeros(a, 2 * abs(zero_sample),
                       'left' if zero_sample > 0 else 'right')
    dif = len(a) - len(b) - 2 * num
    if dif > 0:
        b = _add_zeros(b, dif // 2)
    else:
        a = _add_zeros(a, -dif // 2)
    return correlate(a, b, 'valid')


#def _toeplitz_real_sym(a, b):
#    """
#    Solve linear system Ax=b for real symmetric Toeplitz matrix A.
#    :param a: first row of Toeplitz matrix A
#    :param b: vector b
#    :return: x=A^-1*b
#    """
#    return sto_sl(np.hstack((a, a[1:])), b, job=0)


def plot_events_earth(xlist,ylist,l_bound,u_bound,center=[-32.,123.]):
  """
  plotting routine into which all available events (e.g. from ISC catalog) are plotted on a globe 
  """

  from mpl_toolkits import basemap

  #m = basemap.Basemap(projection='ortho',lon_0=center[1],lat_0=center[0],resolution='i')
  m = basemap.Basemap(projection='robin',lon_0=center[1],resolution='i')
  m.drawcoastlines()

  m.fillcontinents(color='coral',lake_color='aqua')
  m.drawmapboundary(fill_color='aqua')
 
 
  for i in range(len(xlist)):
    #check whether list entries are in "allowed" range
    from obspy.core.util.geodetics import gps2DistAzimuth
    distm = gps2DistAzimuth(ylist[i],xlist[i],center[0],center[1])[0] 
    print(distm / (1000. * KM_PER_DEG))
    if distm / (1000. * KM_PER_DEG) > u_bound or distm / (1000. * KM_PER_DEG) < l_bound:
      print("discarded")
      continue
    print("In")
    x,y = m(xlist[i],ylist[i])
    plot(x,y,'ro',markersize=6)

  #get circles
  X1,Y1 = equi(m,center[1],center[0],l_bound*KM_PER_DEG)
  X1_new,Y1_new = m(X1,Y1)

  X2,Y2 = equi(m,center[1],center[0],u_bound*KM_PER_DEG)
  X2_new,Y2_new = m(X2,Y2)

  plot(X1_new,Y1_new,'k--')
  plot(X2_new,Y2_new,'k--')

  cenx,ceny = m(center[1],center[0])
  plot(cenx,ceny,'b*',markersize=15)

  m.drawparallels([-60,-30,0,30,60],labels=[1,0,0,0])
  m.drawmeridians([-120,-60,0,60,120,180],labels=[0,0,0,1])

  savefig('Blup.pdf')


def equi(m, centerlon, centerlat, radius, *args, **kwargs):
  """
  taken from tutorial by Thomas Lecocq
  """
  glon1 = centerlon
  glat1 = centerlat
  X = []
  Y = []
  for azimuth in range(0, 360):
    glon2, glat2, baz = shoot(glon1, glat1, azimuth, radius)
    X.append(glon2)
    Y.append(glat2)
  X.append(X[0])
  Y.append(Y[0])    
 
  return X,Y

def shoot(lon, lat, azimuth, maxdist=None):
  """Shooter Function
  Original javascript on http://williams.best.vwh.net/gccalc.htm
  Translated to python by Thomas Lecocq
  """
  glat1 = lat * np.pi / 180.
  glon1 = lon * np.pi / 180.
  s = maxdist / 1.852
  faz = azimuth * np.pi / 180.

  EPS= 0.00000000005
  if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
    alert("Only N-S courses are meaningful, starting at a pole!")

  a=6378.13/1.852
  f=1/298.257223563
  r = 1 - f
  tu = r * np.tan(glat1)
  sf = np.sin(faz)
  cf = np.cos(faz)
  if (cf==0):
    b=0.
  else:
    b=2. * np.arctan2 (tu, cf)

  cu = 1. / np.sqrt(1 + tu * tu)
  su = tu * cu
  sa = cu * sf
  c2a = 1 - sa * sa
  x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
  x = (x - 2.) / x
  c = 1. - x
  c = (x * x / 4. + 1.) / c
  d = (0.375 * x * x - 1.) * x
  tu = s / (r * a * c)
  y = tu
  c = y + 1
  while (np.abs (y - c) > EPS):

    sy = np.sin(y)
    cy = np.cos(y)
    cz = np.cos(b + y)
    e = 2. * cz * cz - 1.
    c = y
    x = e * cy
    y = e + e - 1.
    y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
          d / 4. - cz) * sy * d + tu

  b = cu * cy * cf - su * sy
  c = r * np.sqrt(sa * sa + b * b)
  d = su * cy + cu * sy * cf
  glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
  c = cu * cy - su * sy * cf
  x = np.arctan2(sy * sf, c)
  c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
  d = ((e * cy * c + cz) * sy * c + y) * sa
  glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi    

  baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)

  glon2 *= 180./np.pi
  glat2 *= 180./np.pi
  baz *= 180./np.pi

  return (glon2, glat2, baz)

def get_Ps_quartilerange(indir,outfile='medquarts_ALFREX'):
  """
  operates on list of dictionaries that contain Ps picks, for every station, the median and the
  inner quartile range of Ps times is computed and stored (for later plotting on a map etc)
  """

  inlist = glob(indir+'/*__new')
  otf = open(outfile,'w')

  stats = []
  Ps_med = []
  Ps_quarts = []

  for entry in inlist:
    stat = entry.split('/')[-1].split('_')[0]
    stats.append(stat)
    dat = pkl.load(open(entry,'rb'))
    Ps = []
    for j in dat.keys():
      if dat[j].has_key('Ps_time'):
        Ps.append(dat[j]['Ps_time'])
  
    Ps_ar = array(Ps)
    Ps_med.append(median(Ps))
    quart_upper = percentile(Ps,75.)
    quart_lower = percentile(Ps,25.)
    Ps_quarts.append((quart_lower,quart_upper))

  for l in range(len(stats)):
    otf.write(stats[l]+' '+str(round(Ps_med[l],2))+' '+str(round(Ps_quarts[l][0],2))+' '+str(round(Ps_quarts[l][1],2))+'\n')

  return stats, Ps_med, Ps_quarts
  








