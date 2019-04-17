#module util
#coding=utf-8

from pylab import *

import obspy
from obspy.core import UTCDateTime,Trace,Stream,read
from obspy.geodetics.base import gps2dist_azimuth

from copy import deepcopy as cp

import numpy
from numpy import arctan2

import subprocess,glob,os
from subprocess import Popen

import basic

def onclick(event):
  """  
  for capturing mouse clicks on a figure 
  """
  global ix, iy
  ix, iy = event.xdata, event.ydata
  print 'x = %d, y = %d'%(
      ix, iy)
  global coords
  coords = [ix, iy]

  return coords


def resample(self, newfsamp):
    from scipy.signal import resample
    n = round(len(self.data)*newfsamp/self.fsamp)
    outdata=resample(self.data, n, window='blackman')
    self.fsamp = newfsamp
    self.data = outdata

#obsolete
def get_stat_coords(stat,info_file):
  """
  get station coordinates from info_file
  stat: station string
  """
#  info_file = '/RAID/tipage/data/meta_data/info_file'
  cmd = 'awk \'/^'+stat+' / {print $8,$9,$10}\' '+info_file
  process = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
  lat_stat,lon_stat,stat_elev = process.communicate()[0].strip('\n').split(' ')
  lat_stat = float(lat_stat)
  lon_stat = float(lon_stat)
  elev_stat = float(stat_elev)
  return lat_stat,lon_stat,elev_stat

def read_cnv4CC(cnv_file):
  """
  modified (lighter on memory) version for usage by CC routines
  """
  infile = open(cnv_file,'r')
  data = infile.readlines()
  infile.close()
  ev_num=1
  cnv = basic.Catalog()
  for i in range(len(data)):
    if len(data[i]) > 35:
      if data[i][35] == 'e' or data[i][35] == 'E' or data[i][35] == 'w' or data[i][35] == 'W':
	#ev_dict = {}
	yr = '%02d' % (int(data[i][0:2]))
	year = int('20'+str(yr))
	month = int(data[i][2:4])
	day = int(data[i][4:6])
	hour = int(data[i][7:9])
	mn = int(data[i][9:11])
	sec = float(data[i][12:17])

	evlat = float(data[i][18:25])
	evlathem = data[i][25]
	if evlathem == 'S' or evlathem == 's':
	  evlat *= (-1)
	evlon = float(data[i][27:35])
	evlonhem = data[i][35]
	if evlonhem == 'W' or evlonhem == 'w':
	  evlon *= (-1)
	evdep = float(data[i][36:46])
	evmag = float(data[i][46:52])
	evgap = int(data[i][52:58])
	evres = float(data[i][58:].strip('\n'))
	orig_time = UTCDateTime(year,month,day,hour,mn,sec)
	print orig_time
	ev = basic.Event(orig_time,evlat,evlon,evdep)
        ev.Picks = {}
	ev.info['magnitude'] = evmag
	ev.info['gap'] = evgap
	ev.info['RMS_residual'] = evres
	stalist = []
	for a in range((i+1),(i+1000)):
	  if data[a] == ' \n' or data[a] == '\n': #first one for cnv files created by ifort-compiled Velest...
	    break
	  if a > (len(data) + 1):
	    break
	  linestring = data[a].strip('\n')
	  le = len(linestring)
	  str_len = 4+8
	  stat1 = linestring[0:str_len]
	  stalist.append(stat1)
	  if le > str_len:
	    stat2 = linestring[str_len:(2*str_len)]
	    stalist.append(stat2)      
	  else:
	    break  
	  if le > (2*str_len):
	    stat3 = linestring[(2*str_len):(3*str_len)]
	    stalist.append(stat3)        
	  else:
	    break
	  if le > (3*str_len):
	    stat4 = linestring[(3*str_len):(4*str_len)]
	    stalist.append(stat4)        
	  else:
	    break
	  if le > (4*str_len):
	    stat5 = linestring[(4*str_len):(5*str_len)]
	    stalist.append(stat5)        
	  else:
	    break
	  if le > (5*str_len):
	    stat6 = linestring[(5*str_len):(6*str_len)]
	    stalist.append(stat6)        
	  else:
	    break
	for stat in stalist:
	  stat_name = stat[0:4].strip('_')
	  stat_type = stat[4:5]
	  stat_time = float(stat[6:str_len])  
	  Flag = True
          try:
            ev.Picks[stat_name]
	  except KeyError:
            ev.Picks[stat_name] = {}
	    Flag = False
	  if stat_type == 'P' or stat_type == 'p':
	    if Flag:
              ev.Picks[stat_name]['P-Pick'] = stat_time
	    else:
              ev.Picks[stat_name] = {}
              ev.Picks[stat_name]['P-pick'] = stat_time
	  elif stat_type == 's' or stat_type == 'S':
	    if Flag:
              ev.Picks[stat_name]['S-pick'] = stat_time
	    else:
              ev.Picks[stat_name] = {}
              ev.Picks[stat_name]['S-pick'] = stat_time

	#ev.info['cumsum'] = get_cumsum(ev)
	#ev.info['cumsumS'] = get_cumsumS(ev)
        cnv.events[ev_num] = cp(ev)
	ev_num+=1
  cnv.info['Nev'] = ev_num
  cnv.info['source'] = cnv_file
  return cnv

def read_pha(infile):
  """
  reading a ph2dt .pha-file into Catalog format
  """
  cat = basic.Catalog()
  cat.events = {}

  infl = open(infile,'r')
  dat = infl.readlines()
  infl.close()

  for j in range(len(dat)):
    if dat[j][0] == '#': #new event
      print dat[j]
      dum,yr,mo,dy,hr,mn,sc,lat,lon,dep,mag,dum,dum,dum,ID = dat[j].split(None)
      orig = UTCDateTime(int(yr),int(mo),int(dy),int(hr),int(mn),float(sc))
      ev = basic.Event(orig,float(lat),float(lon),float(dep))
      ev.info['eventID'] = ID
      ev.info['magnitude'] = float(mag)
      ev.Picks = {}
      for k in range(j+1,j+1000):
        if k == len(dat):
          break
        if dat[k][0] == '#':
          break
        stat,tt,wgt,typ = dat[k].split(None)
        if not ev.Picks.has_key(stat):
          ev.Picks[stat] = {}
        if typ == 'P':
          ev.Picks[stat]['P-pick'] = float(tt)
          if float(wgt) == 1.:
            ev.Picks[stat]['P-weight'] = 0
          elif float(wgt) == 0.75:
            ev.Picks[stat]['P-weight'] = 1
          elif float(wgt) == 0.5:
            ev.Picks[stat]['P-weight'] = 2
          elif float(wgt) == 0.25:
            ev.Picks[stat]['P-weight'] = 3
        elif typ == 'S':
          ev.Picks[stat]['S-pick'] = float(tt)
          if float(wgt) == 1.0:
            ev.Picks[stat]['S-weight'] = 0
          elif float(wgt) == 0.75:
            ev.Picks[stat]['S-weight'] = 1
          elif float(wgt) == 0.5:
            ev.Picks[stat]['S-weight'] = 2
          elif float(wgt) == 0.25:
            ev.Picks[stat]['S-weight'] = 3

      cat.events[orig.isoformat()] = cp(ev)

  return cat


def write_pha(cat,outfile='blup.pha'):
  """
  produce .pha file for ph2dt from Catalog()
  """
  out = open(outfile,'w')
  num = 1
  for j in sorted(cat.events.keys()):
    orig = cat.events[j].info['origin_time']
    #event line
    out.write('# %4i %02i %02i %02i %02i %05.2f %7.4f %8.4f %6.2f %3.1f 0.0 0.0 0.0 %6s\n' % (orig.year,orig.month,orig.day,orig.hour,orig.minute,(orig.second+(orig.microsecond/1e6)),cat.events[j].info['event_lat'],cat.events[j].info['event_lon'],cat.events[j].info['event_dep'],cat.events[j].info['magnitude'],cat.events[j].info['eventID']))
    for stat in cat.events[j].Picks.keys():
      out.write('%-4s %4.2f %5.2f P\n' % (stat,cat.events[j].Picks[stat]['P-pick'],1-(cat.events[j].Picks[stat]['P-weight']/4.)))
      if cat.events[j].Picks[stat].has_key('S-pick'):
        out.write('%-4s %4.2f %5.2f S\n' % (stat,cat.events[j].Picks[stat]['S-pick'],1-(cat.events[j].Picks[stat]['S-weight']/4.)))  
    num += 1    

  out.close()


def read_cnv(cnv_file,date='date',pol=True,statlen=4,writeID=False):
  """
  function to read the various flavours of CNV files into a dictionary
  cnv_file: string; path to the CNV file to be read
  date: either 'doy' or 'date'; specifies if date is given as julian day (doy) or month and day (date)
  pol: bool; if True pick polarities are expected
  statlen: any number between 2 and 6 (3 or 4 are most common): length of station strings in CNV files
  writeID: if chosen as True, EventID (number) is saved as an additional item in cat.events[x].info

  returns a dictionary of events
  structure: {event_number:{origin time, lat, lon, ..., Picks:{Station1:{P time, P weight, S time, Sweight,...}}}}


  """
  infile = open(cnv_file,'r')
  data = infile.readlines()
  infile.close()
  ev_num=1
  cnv = basic.Catalog()
  for i in range(len(data)):
    if len(data[i]) > 35:
      if data[i][35] == 'e' or data[i][35] == 'E' or data[i][35] == 'w' or data[i][35] == 'W':
	#ev_dict = {}
	yr = '%02d' % (int(data[i][0:2]))
	year = int('20'+str(yr))
	if date == 'doy':
	  doy = int(data[i][3:6])
	elif date == 'date':
	  month = int(data[i][2:4])
	  day = int(data[i][4:6])
	hour = int(data[i][7:9])
	mn = int(data[i][9:11])
	sec = float(data[i][12:17])

	#time tests (seemingly bug in Velest when origin time moves over full minute/hour etc.)
	if sec >= 60.:
	  sec -= 60.
	  mn += 1.
	if mn >= 60:
	  mn -= 60
	  hour += 1
	if hour >= 24:
	  hour -= 24
	  if date == 'doy':
	    doy += 1
	  else:
	    day += 1

	evlat = float(data[i][18:25])
	evlathem = data[i][25]
	if evlathem == 'S' or evlathem == 's':
	  evlat *= (-1)
	evlon = float(data[i][27:35])
	evlonhem = data[i][35]
	if evlonhem == 'W' or evlonhem == 'w':
	  evlon *= (-1)
	evdep = float(data[i][36:46])
	evmag = float(data[i][46:52])
	evgap = int(data[i][52:58])
	evres = float(data[i][58:].strip('\n'))
	if date == 'doy':
	  orig_time = UTCDateTime(str(year)+'-'+str(doy)+'T'+str(hour)+':'+str(mn)+':'+str(sec))
	elif date == 'date':
	  orig_time = UTCDateTime(year,month,day,hour,mn,sec)
	print orig_time
	ev = basic.Event(orig_time,evlat,evlon,evdep)
        ev.Picks = {}
	ev.info['magnitude'] = evmag
	ev.info['gap'] = evgap
	ev.info['RMS_residual'] = evres
	stalist = []
	for a in range((i+1),(i+1000)):
	  if data[a] == ' \n' or data[a] == '\n': #first one for cnv files created by ifort-compiled Velest...
	    break
	  if a > (len(data) + 1):
	    break
	  linestring = data[a].strip('\n')
	  le = len(linestring)
	  if pol:
	    str_len = statlen+9
	  if not pol:
	    str_len = statlen+8
	  stat1 = linestring[0:str_len]
	  stalist.append(stat1)
	  if le > str_len:
	    stat2 = linestring[str_len:(2*str_len)]
	    stalist.append(stat2)      
	  else:
	    break  
	  if le > (2*str_len):
	    stat3 = linestring[(2*str_len):(3*str_len)]
	    stalist.append(stat3)        
	  else:
	    break
	  if le > (3*str_len):
	    stat4 = linestring[(3*str_len):(4*str_len)]
	    stalist.append(stat4)        
	  else:
	    break
	  if le > (4*str_len):
	    stat5 = linestring[(4*str_len):(5*str_len)]
	    stalist.append(stat5)        
	  else:
	    break
	  if le > (5*str_len):
	    stat6 = linestring[(5*str_len):(6*str_len)]
	    stalist.append(stat6)        
	  else:
	    break
	for stat in stalist:
	  stat_name = stat[0:statlen].strip('_')
	  stat_type = stat[statlen:(statlen+1)]
	  if pol:
	    stat_pol = stat[(statlen+1)]
	    stat_weight = int(stat[(statlen+2)])
	    stat_time = float(stat[(statlen+3):str_len])
	  if not pol:
	    stat_weight = int(stat[(statlen+1)])
	    stat_time = float(stat[(statlen+2):str_len])  
	  Flag = True
          try:
            ev.Picks[stat_name]
	  except KeyError:
            ev.Picks[stat_name] = {}
	    Flag = False
	  if stat_type == 'P' or stat_type == 'p':
	    if Flag:
              ev.Picks[stat_name]['P-Pick'] = stat_time
	      if pol:
                ev.Picks[stat_name]['P-pol'] = stat_pol
              ev.Picks[stat_name]['P-weight'] = int(stat_weight)
	    else:
              ev.Picks[stat_name] = {}
	      if pol:
                ev.Picks[stat_name]['P-Pol'] = stat_pol
              ev.Picks[stat_name]['P-pick'] = stat_time
              ev.Picks[stat_name]['P-weight'] = int(stat_weight)
	  elif stat_type == 's' or stat_type == 'S':
	    if Flag:
              ev.Picks[stat_name]['S-pick'] = stat_time
	      if pol:
                ev.Picks[stat_name]['S-pol'] = stat_pol
              ev.Picks[stat_name]['S-weight'] = int(stat_weight)
	    else:
              ev.Picks[stat_name] = {}
	      if pol:
                ev.Picks[stat_name]['S-pol'] = stat_pol
              ev.Picks[stat_name]['S-pick'] = stat_time
              ev.Picks[stat_name]['S-weight'] = int(stat_weight)


	ev.info['cumsum'] = get_cumsum(ev)
	ev.info['cumsumS'] = get_cumsumS(ev)
        if writeID:
          ev.info['eventID'] = ev_num
        cnv.events[orig_time.isoformat()] = cp(ev)
        #JUST FOR NOW!!!
        #cnv.events[ev_num] = cp(ev)
	ev_num+=1
  cnv.info['Nev'] = ev_num
  cnv.info['source'] = cnv_file
  return cnv

def read_ISCsimple(infile):
  """
  reading in a simplified version of ISC (with out without mb)
  """
  fl = open(infile,'r')
  dat = fl.readlines()
  fl.close()

  cat = basic.Catalog()
  cat.events = {}
  for line in dat:
    try:
      dy,time,lat,lon,dep,mag = line.split(None)
    except:
      dy,time,lat,lon,dep = line.split(None)
      mag = 0
    year,month,day = dy.split('/') 
    hour,minute,second = time.split(':')
    orig_time = UTCDateTime(int(year),int(month),int(day),int(hour),int(minute),float(second))

    ev = basic.Event(orig_time,float(lat),float(lon),float(dep))

    ev.Picks = {}
    ev.info['magnitude'] = float(mag)

    cat.events[orig_time.isoformat()] = cp(ev)
  cat.info['Nev'] = len(dat)
  cat.info['source'] = infile
  return cat

def read_CSNcat(infile):
  """
  routine to read in catalog format used by CSN (Centro Sismologico Nacional de Chile)
  """
  fl = open(infile,'r')
  dat = fl.readlines()
  fl.close()

  cat = basic.Catalog()
  cat.events = {}
  for line in dat:
    year = line[1:5]
    month = line[6:8]
    day = line[8:10]
    hour = line[11:13]
    if hour == '  ':
      hour = 0
    minute = line[13:15]
    second = line[16:20]

    if second == '    ':
      second = 0.

    try:
      if second == '60.0':
        if minute == '59':
          orig_time = UTCDateTime(int(year),int(month),int(day),int(hour)+1,0,0.)
        else:
          orig_time = UTCDateTime(int(year),int(month),int(day),int(hour),int(minute)+1,0.)
      else:
        orig_time = UTCDateTime(int(year),int(month),int(day),int(hour),int(minute),float(second))
    except ValueError:
      orig_time = UTCDateTime(int(year),int(month),int(day),int(hour),int(minute),int(second))
    lat = line[23:30]
    lon = line[30:38]
    dep = line[38:43]

    try:
      nsta = int(line[49:51])
    except ValueError:
      nsta = -99
    try:
      mag = float(line[56:59])
    except ValueError:
      mag = 2.5
    if lat == '       ' or lon == '        ' or dep == '     ':
      continue
    if lat == '' or lon == '' or dep == '':
      continue
    try:
      ev = basic.Event(orig_time,float(lat),float(lon),float(dep))
    except ValueError: #incomplete coordinates
      continue
    ev.Picks = {}
    ev.info['nobs'] = nsta
    ev.info['magnitude'] = mag

    cat.events[orig_time.isoformat()] = cp(ev)
  cat.info['Nev'] = len(dat)
  cat.info['source'] = infile

  return cat  



def get_cumsum(ev):
  """
  determine cumsum for one event
  """
  cumsum = 0
  for k in ev.Picks.keys():
    if ev.Picks[k].has_key('P-pick'):
      if int(ev.Picks[k]['P-weight']) == 0:
        cumsum += 4
      elif int(ev.Picks[k]['P-weight']) == 1:
        cumsum += 3
      elif int(ev.Picks[k]['P-weight']) == 2:
        cumsum += 2
      elif int(ev.Picks[k]['P-weight']) == 3:
        cumsum += 1
    if ev.Picks[k].has_key('S-pick'):
      if int(ev.Picks[k]['S-weight']) == 0:
        cumsum += 8
      elif int(ev.Picks[k]['S-weight']) == 1:
        cumsum += 6
      elif int(ev.Picks[k]['S-weight']) == 2:
        cumsum += 4
      elif int(ev.Picks[k]['S-weight']) == 3:
        cumsum += 2

  return cumsum

def get_cumsumS(ev):
  """
  determine cumsum for S arrivals only, for one event
  """
  cumsumS = 0
  for k in ev.Picks.keys():
    if ev.Picks[k].has_key('S-pick'):
      if int(ev.Picks[k]['S-weight']) == 0:
        cumsumS += 8
      elif int(ev.Picks[k]['S-weight']) == 1:
        cumsumS += 6
      elif int(ev.Picks[k]['S-weight']) == 2:
        cumsumS += 4
      elif int(ev.Picks[k]['S-weight']) == 3:
        cumsumS += 2

  return cumsumS

#needs adaption!!
def weedout_evdict(evdict,p_wt=[0,1,2,3,4],s_wt=[0,1,2,3,4]):
  """
  routine for sorting out picks from evdict (output of e.g. read_cnv()) based on their assigned weights
  p_wt: list of allowed P weights (i.e. all picks with weights not in this list will be deleted)
  s_wt: same thing for S
  """
  
  for j in range(1,len(evdict)+1):
    for k in evdict[j]['Picks'].keys():
      if evdict[j]['Picks'][k].has_key('P-weight'):
	if not evdict[j]['Picks'][k]['P-weight'] in p_wt:
	  del evdict[j]['Picks'][k]#delete whole station...
	else:
	  if evdict[j]['Picks'][k].has_key('S-weight'):
	    if not evdict[j]['Picks'][k]['S-weight'] in s_wt:
	      del evdict[j]['Picks'][k]['S-weight'] #only delete all S-pick-related stuff...
	      del evdict[j]['Picks'][k]['S-pick']

  return evdict

def read_catalogue(catalogue,magflag=True,cumsumflag=True):
  """
  read hypo71.py-derived catalogues into dictionary structure
  when reading catalogue output of hypo71.py, mode 'pre', set cumsumflag to False
  """

  fil = open(catalogue,'r')
  data = fil.readlines()
  fil.close()

  num = 1
  cat = basic.Catalog()
  cat.events = {}
  for i in range(len(data)):
    if data[i][0:1] == '1':
      if magflag and cumsumflag:
	julorig,year,mon,day,hr,mn,sec,lat,lon,dep,nobs,gap,rms,cumsum,mag = data[i].strip('\n').split(None)
      elif not magflag and cumsumflag:
	julorig,year,mon,day,hr,mn,sec,lat,lon,dep,nobs,gap,rms,cumsum = data[i].strip('\n').split(None)
      elif magflag and not cumsumflag: #should not exist in normal procedure...
	julorig,year,mon,day,hr,mn,sec,lat,lon,dep,nobs,gap,rms,mag = data[i].strip('\n').split(None)
      else:
	julorig,year,mon,day,hr,mn,sec,lat,lon,dep,nobs,gap,rms = data[i].strip('\n').split(None)
      orig_time = UTCDateTime(float(julorig))
      ev = basic.Event(orig_time,float(lat),float(lon),float(dep))
      ev.Picks = {}
      if magflag:
	ev.info['magnitude'] = float(mag)
      ev.info['gap'] = int(gap)
      ev.info['RMS_residual'] = float(rms)
      ev.info['nobs'] = int(nobs) #the following two are incorporated for eventual filtering purposes
      for j in range((i+1),(i+1000)):
        if data[j][0:1] == '1':
	  num += 1
	  break
	if data[j] == '\n':
	  break
	print j
	try:
	  stat,julP,Pres,Prel,Pw,dum,Ppol,Pweight,julS,Sres,Srel,Sw,dum,dum,Sweight = data[j].strip('\n').split(None)
          print Prel, Srel
	  sflag = True
	except:
	  stat,julP,Pres,Prel,Pw,dum,Ppol,Pweight = data[j].strip('\n').split(None) #no S phase
	  sflag = False
        ev.Picks[stat] = {}
        ev.Picks[stat]['P-pick'] = UTCDateTime(float(julP)) - UTCDateTime(float(julorig))
	ev.Picks[stat]['P-pol'] = Ppol
	ev.Picks[stat]['P-weight'] = Pweight
	ev.Picks[stat]['P-res'] = Pres
	ev.Picks[stat]['Pwgt'] = Pw
	if sflag:
	  ev.Picks[stat]['S-pick'] = Srel
	  ev.Picks[stat]['S-weight'] = Sweight
	  ev.Picks[stat]['S-pol'] = ''
	  ev.Picks[stat]['S-res'] = Sres
	  ev.Picks[stat]['Swgt'] = Sw
   
      ev.info['cumsum'] = get_cumsum(ev)
      ev.info['cumsumS'] = get_cumsumS(ev) 
      cat.events[orig_time.isoformat()] = cp(ev)

  cat.info['Nev'] = num
  cat.info['source'] = catalogue

  return cat


def get_nearest_sample(time,trace):
  """
  routine for finding sample nearest to a given time in a record
  """
  if time < trace.stats.starttime:
    print 'Error: given time lies before file onset'
    return trace.stats.starttime
  elif time > trace.stats.endtime:
    print 'Error: given time lies after file end'
    return trace.stats.endtime
  else:
    inc = round(((time - trace.stats.starttime)/trace.stats.delta),0)
    samptime = trace.stats.starttime + inc*trace.stats.delta

    return samptime


def write_cnv(cat,outfile,pol=True,statlen=4,date='date',velest=False,cumsum_c=-100,cumsumS_c=-100,nobs_c=-100,gap_c=360,rms_c=99.):
  """
  function to write cnv file from dictionary format
  """
  ofile = open(outfile,'w')
  for i in sorted(cat.events.keys()):#sum over events
    if cat.events[i].info['cumsum'] < cumsum_c:
      continue
    if cat.events[i].info['cumsumS'] < cumsumS_c:
      continue
    if len(cat.events[i].Picks.keys()) < nobs_c:
      continue
    if cat.events[i].info['gap'] > gap_c:#
      continue
    if cat.events[i].info['RMS_residual'] > rms_c:
      continue
    if cat.events[i].info['event_lat'] >= 0.:
      lat = cat.events[i].info['event_lat']
      id_NS = 'N'
    else:
      lat = cat.events[i].info['event_lat'] * (-1)
      id_NS = 'S'
    #cat.events[i].info['gap'] = 123
    #cat.events[i].info['RMS_residual'] = 0.99
    if cat.events[i].info['event_lon'] >= 0.:
      lon = cat.events[i].info['event_lon']
      id_EW = 'E'
      if velest:
	id_EW = 'W'
    else:
      lon = cat.events[i].info['event_lon'] * (-1)
      id_EW = 'W'
    if date == 'date':
      eventline = '%2d%02d%02d %2d%2d %5.2f %07.4f%1s %08.4f%1s  %6.2f   %4.2f    %3i     %5.2f\n' % (int((str(cat.events[i].info['origin_time'].year))[2:]),cat.events[i].info['origin_time'].month,cat.events[i].info['origin_time'].day,cat.events[i].info['origin_time'].hour,cat.events[i].info['origin_time'].minute,(cat.events[i].info['origin_time'].second+(cat.events[i].info['origin_time'].microsecond/1e6)),lat,id_NS,lon,id_EW,cat.events[i].info['event_dep'],cat.events[i].info['magnitude'],cat.events[i].info['gap'],cat.events[i].info['RMS_residual'])
    elif date == 'doy':
      eventline = '%2d %3d %2d%2d %5.2f %07.4f%1s %08.4f%1s  %6.2f   %4.2f    %3i     %5.2f\n' % (int((str(cat.events[i].info['origin_time'].year))[2:]),int(cat.events[i].info['origin_time'].julday),cat.events[i].info['origin_time'].hour,cat.events[i].info['origin_time'].minute,(cat.events[i].info['origin_time'].second+(cat.events[i].info['origin_time'].microsecond/1e6)),lat,id_NS,lon,id_EW,cat.events[i].info['event_dep'],cat.events[i].info['magnitude'],cat.events[i].info['gap'],cat.events[i].info['RMS_residual'])
    ofile.write(eventline)
    phase_strings = []
    for k in cat.events[i].Picks.keys():
      if not cat.events[i].Picks[k].has_key('P-pick'):
	continue
      if statlen == 4:
	if len(k) == 2:
	  ky = k+'__'
	elif len(k) == 3:
	  ky = k+'_'
	elif len(k) == 4:
	  ky = k
	else:
	  ky = k[0:4]
      elif statlen == 3:
	if len(k) == 2:
	  ky = k+'_'
	elif len(k) == 3:
	  ky = k
	else:
	  ky = k[0:3]
      try:
	assert cat.events[i].Picks[k]['P-pick'] != '*****'
      except AssertionError:
	continue
      if statlen == 4:
	if not pol:
          print cat.events[i].Picks[k]['P-weight'],cat.events[i].Picks[k]['P-pick']
	  string1 = '%4sp%1i%6.2f' % (ky,int(cat.events[i].Picks[k]['P-weight']),cat.events[i].Picks[k]['P-pick'])
	else:
	  string1 = '%4sp%1s%1i%6.2f' % (ky,cat.events[i].Picks[k]['P-pol'],cat.events[i].Picks[k]['P-weight'],cat.events[i].Picks[k]['P-pick'])
      elif statlen == 3:
	if not pol:
	  string1 = '%3sp%1i%6.2f' % (ky,cat.events[i].Picks[k]['P-weight'],cat.events[i].Picks[k]['P-pick'])
	else:
	  string1 = '%3sp%1s%1i%6.2f' % (ky,cat.events[i].Picks[k]['P-pol'],cat.events[i].Picks[k]['P-weight'],cat.events[i].Picks[k]['P-pick'])
      phase_strings.append(string1)
      if cat.events[i].Picks[k].has_key('S-pick'):
	try:
	  assert cat.events[i].Picks[k]['S-pick'] != '******'
	except AssertionError:
	  continue
	if not pol:
	  string2 = '%4ss%1i%6.2f' % (ky,int(cat.events[i].Picks[k]['S-weight']),float(cat.events[i].Picks[k]['S-pick']))
	else:
	  string2 = '%4ss %1i%6.2f' % (ky,cat.events[i].Picks[k]['S-weight'],cat.events[i].Picks[k]['S-pick'])
	phase_strings.append(string2)
    numb = len(phase_strings)
    for hu in range(len(phase_strings)):
      ofile.write(phase_strings[hu])
      if hu%6 == 5:
	ofile.write('\n')
    if numb%6 == 0:
      ofile.write('\n')
    else:
      ofile.write('\n\n')

#obsolete!
def find_resp(stat,chn,info_file):
  """
  extracts response information in obspy's dictionary format for given station + channel using information stored in ext.py and info_file
  stat: to be given as 2-to-4-digit string (has to be identical to station name contained in info_file)
  chn: string, either 'E', 'N' or 'Z'
  info_file: path to info_file
  """
  from resp import *
  if chn == 'E':
    cmd = 'awk \'/^/ { if ($1 == \"'+stat+'\") print $15}\' '+info_file
  elif chn == 'N':
    cmd = 'awk \'/^/ { if ($1 == \"'+stat+'\") print $16}\' '+info_file
  elif chn == 'Z':
    cmd = 'awk \'/^/ { if ($1 == \"'+stat+'\") print $17}\' '+info_file
  else:
    print 'Channel does not exist!!'
  proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
  code = proc.communicate()[0]
  cmd2 = 'awk \'/^/ { if ($1 == \"'+stat+'\") print $11}\' '+info_file
  proc2 = subprocess.Popen(cmd2,shell=True,stdout=subprocess.PIPE)
  code2 = proc2.communicate()[0]
  prg = float(code2)
  resp = {}
  resp['poles'] = eval(code).p
  resp['zeros'] = eval(code).z[1:]  # displacement --> velocity response
  resp['gain'] = eval(code).norm
  resp['sensitivity'] = eval(code).gain * prg #prg
  return resp

def correct_WA(data,resp,waterlevel):
  """
  simulate Wood-Anderson displacement sensor for input data
  data: Trace object of seismic (velocity) data 
  resp: response of sensor with which that data was recorded
  waterlevel: waterlevel for restitution calculations
  """
  samp_rate = data.stats.sampling_rate
  PAZ_WOOD_ANDERSON = {
      'poles': [-6.2832 - 4.7124j,- 6.2832 + 4.7124j],
      'zeros': [0.0 + 0.0j] * 1,
      'sensitivity': 2080., #a value of 2800 can also be found in literature
      'gain': 1.
  }

  data.data = data.data - data.data.mean()
  event_out = data.simulate(paz_remove=resp, paz_simulate=PAZ_WOOD_ANDERSON,remove_sensitivity=True,simulate_sensitivity=True)
  return event_out

def amplitude_finder(trace):
  """
  small helper function that selects maximum zero-to-peak amplitude for a given trace
  """
  max = trace.max()
  min = trace.min()

  if abs(max) > abs(min):
    amp = max
  else:
    amp = abs(min)
  return amp

def ml_calculator(hyp_dist,ep_dist,amp,depth,type=4):
  """
  calculate local magnitude from amplitude (Wood-Anderson trace) and distance to event
  hyp_dist: hypocentral distance from station to event (in km)
  ep_dist: epicentral distance from station to event (event depth ignored)
  type: quite a zoo of different formulae for ML exists in the literature; Four of these can currently be selected here
     1: presumably after Bullen and Bolt (1985); found in Shearer: Introduction to seismology; NOT RECOMMENDED!! Comparison with other three formulae and picked magnitudes from GIANT/PITSA shows that these magnitudes are systematically too high by roughly 0.5!!
     2: as derived by Bakun and Joyner (1984); for central California
     3: after Shin (1993), derived for Taiwan, with separate formula for deep earthquakes
     4: after IASPEI recommendation; attention: A in nm and static magnification of 1 for the Wood-Anderson!! (default)
  returns ML value;
  """  
  if type == 1:
    ML = log10(amp) + 2.56*log10(ep_dist) - 1.67
  elif type == 2:
    ML = log10(amp) + log10(hyp_dist) +0.00301*hyp_dist + 0.70
  elif type == 3:
    if depth < 35.:
      if ep_dist < 80.:
        ML = log10(amp) + log10(hyp_dist) + 0.00716*hyp_dist + 0.39
      else:
        ML = log10(amp) + 0.83*log10(hyp_dist) + 0.00261*hyp_dist + 1.07
    else:
      ML = log10(amp) + 0.83*log10(hyp_dist) + 0.00326*hyp_dist + 1.01
  elif type == 4:
    ML = log10(amp*(1e6/2080.)) + 1.11*log10(hyp_dist) + 0.00189*hyp_dist - 2.09
  return ML

def av_remove_outlier(inlist): 
  """
  removes all magnitudes from the input list that differ more than 2 standard deviation from its mean
  inlist: list of ML magnitudes
  """
  length = len(inlist)
  mean = array(inlist,dtype=float).mean()
  stand = std(inlist)
  for f in range(length):
    if abs(float(inlist[f]) - mean) > 2*stand:
      inlist[f] = 999.
  g = True
  while g == True:
    try:
      inlist.remove(999.)
    except:
      g = False
  return inlist

def mvavg(x,l,ab=False):
  """
  efficient computation of moving average over numpy array object
  x: numpy.array object
  l: length of moving window
  """
  n = len(x)
  a = np.zeros(n-l+1,float)
  for j in range(l):
    if ab:
      a += abs(x[j:n-l+1+j])
    else:
      a += x[j:n-l+1+j]
  a /= float(l)
  return a 

def fit_poly(arr_x,arr_y,deg):
  """
  get values of fitted polynomial for x array
  """
  coeffs = numpy.polyfit(arr_x,arr_y,deg)
  arr_return = numpy.polyval(coeffs,arr_x)
  return arr_return

def ms_cut(stt,ett,network,filepath):
  """
  efficient method to cut events from continuous MiniSEED data
  """
  starttime = UTC2time(stt)
  endtime = UTC2time(ett)
  try:
    data = read_ms(filepath,starttime,endtime)
    data.detrend()
    trace = extTrace2Trace(data,network)
    return trace
  except:
    return [] #better error handling needed (raise an exception or something??)

#needs adaption!!
def cut_events(cat,outfolder,data_path,info_file,restitute='None',dur=5,comp='all',mode='all'): #very slow...better: use record-wise cutting
  """
  function to cut events from raw data based on dictionary entries (as read in, for example, by read_cnv or read_catalogue)
  dict: input dictionary of events, format identical to the ones produced by the routines read_cnv or read_catalogue
  outfolder: directory to where cut events should be written (each event gets its own subfolder inside this folder)
  datapath: path to where the raw data is stored
  info_file: path to info file where station and response information is stored
  resititute: 'None' for simply cutting the data, if restitution to a certain sensor/digitizer combination is desired enter a station possessing this combination as a string here (has to be listed in info_file)
  dur: minutes of data to be cut (from origin time onwards)
  comp: 'all', 'N','E' or 'Z'
  mode: 'all' cuts event for all available stations, 'picked' only for those stations that have a P-pick in dict

  """
  for ev_name in cat.events.keys():
    orig_time = cat.events[ev_name].info['origin_time']
    doy = str(orig_time.julday)
    if len(doy) == 1:
      doy = '00'+doy
    if len(doy) == 2:
      doy = '0'+doy
    datestring = '%02i%02i%02i' %((orig_time.year)-2000.,orig_time.month,orig_time.day)
    datestring2 = '%02i%02i%02i' %(((orig_time-86400.).year)-2000.,(orig_time-86400.).month,(orig_time-86400.).day)
    if mode == 'all':
      if comp == 'all':
        data = glob.glob(data_path+'/'+str(orig_time.year)+'/*/*/?[H,N]?.D/*.'+doy)
      elif comp == 'N':
        data = glob.glob(data_path+'/'+str(orig_time.year)+'/*/*/?HN.D/*.'+doy)
      elif comp == 'E':
        data = glob.glob(data_path+'/'+str(orig_time.year)+'/*/*/?HE.D/*.'+doy)
      elif comp == 'Z':
        data = glob.glob(data_path+'/'+str(orig_time.year)+'/*/*/?HZ.D/*.'+doy)

      print data       
      stt = orig_time - 60.
      ett = orig_time + (dur*60)
      shour = '%02d' % (stt.hour)
      smin = '%02d' % (stt.minute)
      ssec = '%02d' % (stt.second)
      ehour = '%02d' % (ett.hour)
      emin = '%02d' % (ett.minute)
      esec = '%02d' % (ett.second)
      folderstring = '%04d%02d%02d_%02d%02d%02d' % (stt.year,stt.month,stt.day,stt.hour,stt.minute,stt.second)
      cmd1 = 'mkdir '+outfolder+'/'+folderstring
      os.system(cmd1)
      if restitute == 'None': #simply cut them out
        for file in data:
          #if stt.hour == 0 and stt.minute <= 20: 
          #  dd = read(file)
          #  stat = dd[0].stats.station
          #  cmp = dd[0].stats.channel
          #  try:
          #    dd2 = read(glob.glob(data_path+'/*/'+stat+datestring2+'*'+cmp)[0]) 
          #    for q in dd2:
          #      dd += q
          #  except:
          #    continue
          dd = read(file)
          stat = dd[0].stats.station
          cmp = dd[0].stats.channel
          aa = dd.slice(starttime=stt,endtime=ett)
          aa.detrend(type='demean')
          aa.detrend(type='linear')
          for d in aa:
            d.stats.mseed.dataquality = 'D'
            d.data = d.data.astype('int32')
          try:
            #aa.write(outfolder+'/'+folderstring+'/'+stat+datestring+'__'+shour+smin+ssec+'.BH'+aa[0].stats.channel[-1],format='MSEED',encoding='STEIM2',reclen=8192)
            aa.write(outfolder+'/'+folderstring+'/'+stat+'.BH'+aa[0].stats.channel[-1],format='MSEED',encoding='STEIM2',reclen=8192)
          except:
            continue        
      else: # restitution desired --> reading in, simulation, write temporary dayfile, then ms_extract
        for file in data:
          network = file.split('/')[5]
          #st_xx = ms_cut(stt,ett,network,file)

          st_new = read(file,format='MSEED')
          #st_new = st.slice(stt,ett)
#          try:
#            assert len(st_new) != 0
#          except AssertionError: #day file exists, but has gap at time of the event
#            continue
#         st_new = Stream(traces=[st_xx])
#         st_new.merge(method=0,fill_value='latest')
          st_new1 = st_new.slice(stt,ett)
          st_new1.merge(method=0,fill_value='latest')
          try:
            assert len(st_new1) != 0
          except AssertionError:
            continue
          paz_old = find_resp(str(st_new1[0].stats.station),str(st_new1[0].stats.channel)[2],info_file)
          paz_new = find_resp(restitute,'Z',info_file)
          st_new1.simulate(paz_remove=paz_old,paz_simulate=paz_new)
          st_new1[0].data *= 1000.
          st_new1[0].data = st_new1[0].data.astype('int32')
          namestring = st_new1[0].stats.station+'.'+st_new1[0].stats.network+'.'+st_new1[0].stats.channel+'.D.'+str(stt.year)+'.'+doy+'.'+shour+smin
          try:
            st_new1.write(outfolder+'/'+folderstring+'/'+namestring,format='MSEED',encoding='STEIM2')
          except:
            continue
  #mode = picked --> still to be implemented
  """
      if mode == 'picked':
      data = []
      for stat in dict[ev_nr]['Picks'].keys():
        if comp == 'all':
          dat = glob.glob(data_path+'/'+str(orig_time.year)+'/*/'+stat+'/?H?.D/*.'+doy)
        elif comp == 'N':
          dat = glob.glob(data_path+'/'+str(orig_time.year)+'/*/'+stat+'/?HN.D/*.'+doy)
        elif comp == 'E':
          dat = glob.glob(data_path+'/'+str(orig_time.year)+'/*/'+stat+'/?HE.D/*.'+doy)
        elif comp == 'Z':
          dat = glob.glob(data_path+'/'+str(orig_time.year)+'/*/'+stat+'/?HZ.D/*.'+doy)

        for uu in range(len(dat)):
          data.append(dat)
  """


def info_file2infodat(info_file,datfile,networks='all'):
  """
  routine for creating a valid info.dat file from info_file format
  input: paths to info file and desired info.dat-like output file
  networks: either choose 'all' or give a list of the network strings whose stations should go into the info.dat file
  """
  if networks == 'all':
    cmd = 'awk \'/^/ {print $0}\' '+info_file
    proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    output = proc.communicate()[0]
    lines = output.split('\n')[1:len(output.split('\n'))-1]
  
  else:
    lines = []
    for n in networks:
      cmd = 'awk \'/^/ {if ($14 == \"'+n+'\") print $0}\' '+info_file
      proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
      output = proc.communicate()[0]
      ll = output.split('\n')[0:len(output.split('\n'))-1]
      for xy in ll:
        lines.append(xy)


  outfile = open(datfile,'w')
  index1 = 1
  index2 = 2000
  index3 = 1000

  for l in lines:
    try:
      stat,sens,dum,dum,datl,dum,dum,lat,lon,elev,fac,dum,dum,network,dum,dum,dum = l.split(None)
    except:
      stat,sens,dum,dum,datl,dum,dum,lat,lon,elev,fac,dum,dum,network,dum,dum,dum,dum = l.split(None)
    outstring = '%4d  %-4s  1.  32  %-6s  %4d %03d  %-4s  %08.5f  %8.5f  %4i  2008.06.01  14:00  2010.12.31  23:00  TIPAGE' % (index2,datl,sens,index3,index1,stat,float(lat),float(lon),int(elev))
    outfile.write(outstring+'\n')
    index1 += 1
    index2 += 1
    index3 += 1

  outfile.close()



def write_pazfile(info_file,outfolder,networks='all'):
  """
  automatically generate GIANT/PITSA paz files from entries in info_file 
  """
  if networks == 'all':
    cmd = 'awk \'/^/ {print $1,$15,$16,$17}\' '+info_file
    proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    output = proc.communicate()[0]
    lines = output.split('\n')[1:len(output.split('\n'))-1]

  else:
    lines = []
    for n in networks:
      cmd = 'awk \'/^/ {if ($14 == \"'+n+'\") print $1,$15,$16,$17}\' '+info_file
      proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
      output = proc.communicate()[0]
      ll = output.split('\n')[0:len(output.split('\n'))-1]
      for xy in ll:
        lines.append(xy)

  index1 = 1
  index2 = 2000
  index3 = 1000

  for l in lines:
    stat,respE,respN,respZ = l.split(None)
  
    basename = '%04d_%03d' % (index3,index1)
    respfileE = open(basename+'.e','w')
    sig = 3
    outstring = 'CAL1    %-4s   %3d_%1i   PDAS PAZ 071231 1430 111010 2300\n2\n%8.5f%8.5f\n%8.5f%8.5f\n3\n0.00000 0.00000\n0.00000 0.00000\n0.00000 0.00000\n%6.4f' % (stat,index1,sig,find_resp(stat,'E',info_file)['poles'][0].real,find_resp(stat,'E',info_file)['poles'][0].imag,find_resp(stat,'E',info_file)['poles'][1].real,find_resp(stat,'E',info_file)['poles'][1].imag,find_resp(stat,'E',info_file)['sensitivity'])
    respfileE.write(outstring)
    respfileE.close()
    respfileN = open(basename+'.n','w')
    sig = 2
    outstring = 'CAL1    %-4s   %3d_%1i   PDAS PAZ 071231 1430 111010 2300\n2\n%8.5f%8.5f\n%8.5f%8.5f\n3\n0.00000 0.00000\n0.00000 0.00000\n0.00000 0.00000\n%6.4f' % (stat,index1,sig,find_resp(stat,'N',info_file)['poles'][0].real,find_resp(stat,'N',info_file)['poles'][0].imag,find_resp(stat,'N',info_file)['poles'][1].real,find_resp(stat,'N',info_file)['poles'][1].imag,find_resp(stat,'N',info_file)['sensitivity'])
    respfileN.write(outstring)
    respfileN.close()
    respfileZ = open(basename+'.z','w')
    sig = 1
    outstring = 'CAL1    %-4s   %3d_%1i   PDAS PAZ 071231 1430 111010 2300\n2\n%8.5f%8.5f\n%8.5f%8.5f\n3\n0.00000 0.00000\n0.00000 0.00000\n0.00000 0.00000\n%6.4f' % (stat,index1,sig,find_resp(stat,'Z',info_file)['poles'][0].real,find_resp(stat,'Z',info_file)['poles'][0].imag,find_resp(stat,'Z',info_file)['poles'][1].real,find_resp(stat,'Z',info_file)['poles'][1].imag,find_resp(stat,'Z',info_file)['sensitivity'])
    respfileZ.write(outstring)
    respfileZ.close()
    index1 += 1
    index2 += 1
    index3 += 1

  print 'response files written (hopefully correct...)'

#needs adaption!!
def read_cat_NEIC(cat,stt,ett,dist=True,sel_mag=[0.,10.],sel_az=[0.,360.],centr_lat=0,centr_lon=0,sel_lat=[-90.,90.],sel_lon=[-180.,180.],sel_dep=[0.,1000.],sel_dist=[0.,40000.]):
  """
  function for reading in (and filtering) NEIC catalogue and returning its contents as dictionary format
  for all kwargs enter two-element list of min and max desired value (always use floats)
  """
  infile = open(cat,'r')
  data = infile.readlines()
  infile.close()
  uha = 0
  #get rid of header
  for i in range(len(data)):
    if data[i][0:4] == ' PDE':
      uha = i
      break
  
  dict_out = {}
  num = 1
  for j in range(uha,len(data)):
    if data[j] == '\n':
      break
    if dist:
      cata,year,mon,day,origtime,lat,lon,dep,mag,magstring,infostring,xystring,dist = data[j].strip('\n').split(None)
    else:
      cata,year,mon,day,origtime,lat,lon,dep,mag,magstring,infostring,xystring = data[j].strip('\n').split(None)
    hr = origtime[0:2]
    min = origtime[2:4]
    sec = origtime[4:]
    orig = UTCDateTime(int(year),int(mon),int(day),int(hr),int(min),float(sec))
    az = gps2dist_azimuth(centr_lat,centr_lon,float(lat),float(lon))[1]
    if orig > stt and orig < ett: #right time period
      if float(mag) > sel_mag[0] and float(mag) < sel_mag[1] and float(lon) > sel_lon[0] and float(lon) < sel_lon[1] and float(lat) > sel_lat[0] and float(lat) < sel_lat[1] and int(dep) > sel_dep[0] and int(dep) < sel_dep[1] and float(az) > sel_az[0] and float(az) < sel_az[1]: #optional filtering parameters
        if dist:
          if float(dist) > sel_dist[0] and float(dist) < sel_dist[1]:
            dict_out[num] = {}
            dict_out[num]['origin_time'] = orig
            dict_out[num]['event_lat'] = lat
            dict_out[num]['event_lon'] = lon
            dict_out[num]['event_depth'] = dep
            dict_out[num]['magnitude'] = mag
            num += 1
        else:
          dict_out[num] = {}
          dict_out[num]['origin_time'] = orig
          dict_out[num]['event_lat'] = lat
          dict_out[num]['event_lon'] = lon
          dict_out[num]['event_depth'] = dep
          dict_out[num]['magnitude'] = mag
          num += 1
  return dict_out

#needs adaption!!
def cut_teleseis(dict,datadir,outfolder,info_file,rest=False,paz_new=0,phase='P',length=240,networks='all',rotate=False):
  for i in dict.keys():
    origtime = dict[i]['origin_time']
    lat_ev = dict[i]['event_lat']
    lon_ev = dict[i]['event_lon']
    dep_ev = dict[i]['event_depth']
    doy = '%03d' % (origtime.julday)
    #create subfolder for output
    fname = str(origtime.year)+'_'+str(origtime.month)+'_'+str(origtime.day)+'_'+str(origtime.hour)+'_'+str(origtime.minute)+'_'+str(origtime.second)
    cmdx = 'mkdir '+outfolder+'/'+fname
    os.system(cmdx)
    if networks == 'all':
      datalist = glob.glob(datadir+'/'+str(origtime.year)+'/*/*/*/*.'+doy)
      for stat in datalist:
        #filename = os.path.basename(stat)
        statID = os.path.basename(stat).split('.')[0]
        nw = os.path.basename(stat).split('.')[1] 
        comp = os.path.basename(stat).split('.')[2]
        stat_lat,stat_lon,stat_elev = get_stat_coords(statID,info_file)
        dist = gps2dist_azimuth(float(lat_ev),float(lon_ev),float(stat_lat),float(stat_lon))[0]/1000.
        #raytracing
        if phase == 'P':
          ttime = ttt.firstArrival(dist,float(dep_ev)).time
          centr_time = origtime + float(ttime)
        elif phase == 'S':
          try:
            ttime = ttt.compute(dist,float(dep_ev)).findPhase('S').time
          except:
            print 'No direct S phase exists for station '+statID
            break
          centr_time = origtime + float(ttime)
        early = centr_time - (length/2.) 
        late = centr_time + (length/2.)      
      #  trace = ms_cut(early,late,nw,stat) 
        stream = read(stat,format='MSEED',starttime=early,endtime=late,nearest_sample=True)
        if len(stream) > 1:
          stream = stream.merge()
        try:
          trace = stream[0]
        except:
          continue
        if phase == 'P':
          if comp[2] == 'Z':
            if rest:
              resp_old = find_resp(statID,comp[2],'/RAID/tipage/data/meta_data/info_file')
              trace.simulate(paz_remove=resp_old,paz_simulate=paz_new)
              #MSEED-switch...
              #for ii in range(len(trace)):
              #  trace[ii].data = trace[ii].data.astype('int32')*10000.
            trace.write(outfolder+'/'+fname+'/'+str(statID)+'_'+str(nw)+'_'+comp+'_'+fname+'_P.sac','SAC')
#       re-read,slice correctly, restitute if desired
          else:
            continue
        elif phase == 'S':
          if rest:
            resp_old = find_resp(statID,comp[2],'/RAID/tipage/data/meta_data/info_file')
            trace.simulate(paz_remove=resp_old,paz_simulate=paz_new)
            #MSEED-switch
            #for ii in range(len(trace)):
            #  trace[ii].data = trace[ii].data.astype('int32')*10000.
          trace.write(outfolder+'/'+fname+'/'+str(statID)+'_'+str(nw)+'_'+comp+'_'+fname+'_S.sac','SAC')
        #  data_new = read(outfolder+'/'+fname+'/'+str(statID)+'_'+str(nw)+'_'+comp+'_'+fname+'_S')
        #dnew = data_new.slice(early,late)
        #if rest:
         # resp_old = find_resp(statID,comp[2],'/RAID/tipage/data/meta_data/info_file')
         # dnew.simulate(paz_remove=resp_old,paz_simulate=paz_new)
          #for ii in range(len(dnew)):
           # dnew[ii].data = dnew[ii].data.astype('int32')
        #if phase == 'P':
        #  dnew.write(outfolder+'/'+fname+'/'+str(statID)+'_'+str(nw)+'_'+comp+'_'+fname+'_P','MSEED')
        #elif phase == 'S':
        #  dnew.write(outfolder+'/'+fname+'/'+str(statID)+'_'+str(nw)+'_'+comp+'_'+fname+'_S','SAC')  



def conf_matrix(normlist_1,normlist_2,normlist_3,normlist_4,normlist_1_delta,normlist_2_delta,normlist_3_delta,normlist_4_delta,list5_norm,path,name):
  """
  function for plotting confusion matrices...
  input parameters: normlist_1 to normlist_4: lists of picks for MPX classes 1 to 4, contains Hand quality classes of the same picks
                    list5_norm: same thing for the traces where MPX identified no pick at all
                    normlist_X_delta: standard deviation of the single picks for the different MPX quality classes
                    path, name: where the resulting pdf should be saved to

  """
  clf()
  bar(0,4,6,0,color='r')
  bar(6,4,6,0,color='r')
  bar(12,4,6,0,color='k')
  bar(18,4,6,0,color='k')
  bar(24,4,6,0,color='w')
  bar(0,4,6,4,color='r')
  bar(6,4,6,4,color='k')
  bar(12,4,6,4,color='k')
  bar(18,4,6,4,color='w')
  bar(24,4,6,4,color='b')
  bar(0,4,6,8,color='k')
  bar(6,4,6,8,color='k')
  bar(12,4,6,8,color='w')
  bar(18,4,6,8,color='b')
  bar(24,4,6,8,color='b')
  bar(0,4,6,12,color='k')
  bar(6,4,6,12,color='w')
  bar(12,4,6,12,color='b')
  bar(18,4,6,12,color='b')
  bar(24,4,6,12,color='b')
  bar(0,4,6,16,color='w')
  bar(6,4,6,16,color='b')
  bar(12,4,6,16,color='b')
  bar(18,4,6,16,color='b')
  bar(24,4,6,16,color='b')

  xlabel('classification MPX')
  ylabel('classification Hand')

  xticks([3,9,15,21,27],['class 0','class 1','class 2','class 3','not picked'])
  yticks([2,6,10,14,18],['class 4','class 3','class 2','class 1','class 0'])

  # cumulative pick number per field
  field_0_1_count = 0
  field_0_2_count = 0
  field_0_3_count = 0
  field_0_4_count = 0
  field_0_5_count = 0

  field_1_1_count = 0
  field_1_2_count = 0
  field_1_3_count = 0
  field_1_4_count = 0
  field_1_5_count = 0

  field_2_1_count = 0
  field_2_2_count = 0
  field_2_3_count = 0
  field_2_4_count = 0
  field_2_5_count = 0

  field_3_1_count = 0
  field_3_2_count = 0
  field_3_3_count = 0
  field_3_4_count = 0
  field_3_5_count = 0

  field_4_1_count = 0
  field_4_2_count = 0
  field_4_3_count = 0
  field_4_4_count = 0
  field_4_5_count = 0

  # all delta values per field...for subsequent mean calculation

  field_0_1_delta = []
  field_0_2_delta = []
  field_0_3_delta = []
  field_0_4_delta = []

  field_1_1_delta = []
  field_1_2_delta = []
  field_1_3_delta = []
  field_1_4_delta = []

  field_2_1_delta = []
  field_2_2_delta = []
  field_2_3_delta = []
  field_2_4_delta = []

  field_3_1_delta = []
  field_3_2_delta = []
  field_3_3_delta = []
  field_3_4_delta = []

  field_4_1_delta = []
  field_4_2_delta = []
  field_4_3_delta = []
  field_4_4_delta = []

  for a in range(len(normlist_1)):
    if normlist_1[a] == 0:
      field_0_1_count += 1
      field_0_1_delta.append(normlist_1_delta[a])
    elif normlist_1[a] == 1:
      field_1_1_count += 1
      field_1_1_delta.append(normlist_1_delta[a])
    elif normlist_1[a] == 2:
      field_2_1_count += 1
      field_2_1_delta.append(normlist_1_delta[a])
    elif normlist_1[a] == 3:
      field_3_1_count += 1
      field_3_1_delta.append(normlist_1_delta[a])
    elif normlist_1[a] == 4:
      field_4_1_count += 1
      field_4_1_delta.append(normlist_1_delta[a])


  for b in range(len(normlist_2)):
    if normlist_2[b] == 0:
      field_0_2_count += 1
      field_0_2_delta.append(normlist_2_delta[b])
    elif normlist_2[b] == 1:
      field_1_2_count += 1
      field_1_2_delta.append(normlist_2_delta[b])
    elif normlist_2[b] == 2:
      field_2_2_count += 1
      field_2_2_delta.append(normlist_2_delta[b])
    elif normlist_2[b] == 3:
      field_3_2_count += 1
      field_3_2_delta.append(normlist_2_delta[b])
    elif normlist_2[b] == 4:
      field_4_2_count += 1
      field_4_2_delta.append(normlist_2_delta[b])


  for c in range(len(normlist_3)):
    if normlist_3[c] == 0:
      field_0_3_count += 1
      field_0_3_delta.append(normlist_3_delta[c])
    elif normlist_3[c] == 1:
      field_1_3_count += 1
      field_1_3_delta.append(normlist_3_delta[c])
    elif normlist_3[c] == 2:
      field_2_3_count += 1
      field_2_3_delta.append(normlist_3_delta[c])
    elif normlist_3[c] == 3:
      field_3_3_count += 1
      field_3_3_delta.append(normlist_3_delta[c])
    elif normlist_3[c] == 4:
      field_4_3_count += 1
      field_4_3_delta.append(normlist_3_delta[c])

  for d in range(len(normlist_4)):
    if normlist_4[d] == 0:
      field_0_4_count += 1
      field_0_4_delta.append(normlist_4_delta[d])
    elif normlist_4[d] == 1:
      field_1_4_count += 1
      field_1_4_delta.append(normlist_4_delta[d])
    elif normlist_4[d] == 2:
      field_2_4_count += 1
      field_2_4_delta.append(normlist_4_delta[d])
    elif normlist_4[d] == 3:
      field_3_4_count += 1
      field_3_4_delta.append(normlist_4_delta[d])
    elif normlist_4[d] == 4:
      field_4_4_count += 1
      field_4_4_delta.append(normlist_4_delta[d])


  for e in range(len(list5_norm)):
    if list5_norm[e] == 0:
      field_0_5_count += 1
    elif list5_norm[e] == 1:
      field_1_5_count += 1
    elif list5_norm[e] == 2:
      field_2_5_count += 1
    elif list5_norm[e] == 3:
      field_3_5_count += 1
    elif list5_norm[e] == 4:
      field_4_5_count += 1


  text(1.5,2.0,'N = '+str(field_4_1_count)+'\n('+str(round(field_4_1_count/(field_4_1_count+field_4_2_count+field_4_3_count+field_4_4_count+field_4_5_count+0.0001)*100.,1))+'%)',color='w')
  text(7.5,2.0,'N = '+str(field_4_2_count)+'\n('+str(round(field_4_2_count/float(field_4_1_count+field_4_2_count+field_4_3_count+field_4_4_count+field_4_5_count+0.0001)*100.,1))+'%)',color='w')
  text(13.5,2.0,'N = '+str(field_4_3_count)+'\n('+str(round(field_4_3_count/float(field_4_1_count+field_4_2_count+field_4_3_count+field_4_4_count+field_4_5_count+0.0001)*100.,1))+'%)',color='w')
  text(19.5,2.0,'N = '+str(field_4_4_count)+'\n('+str(round(field_4_4_count/float(field_4_1_count+field_4_2_count+field_4_3_count+field_4_4_count+field_4_5_count+0.0001)*100.,1))+'%)',color='w')
  text(25.5,2.0,'N = '+str(field_4_5_count)+'\n('+str(round(field_4_5_count/float(field_4_1_count+field_4_2_count+field_4_3_count+field_4_4_count+field_4_5_count+0.0001)*100.,1))+'%)',color='k')

  text(1.5,6.0,'N = '+str(field_3_1_count)+'\n('+str(round(field_3_1_count/float(field_3_1_count+field_3_2_count+field_3_3_count+field_3_4_count+field_3_5_count+0.0001)*100.,1))+'%)',color='w')
  text(7.5,6.0,'N = '+str(field_3_2_count)+'\n('+str(round(field_3_2_count/float(field_3_1_count+field_3_2_count+field_3_3_count+field_3_4_count+field_3_5_count+0.0001)*100.,1))+'%)',color='w')
  text(13.5,6.0,'N = '+str(field_3_3_count)+'\n('+str(round(field_3_3_count/float(field_3_1_count+field_3_2_count+field_3_3_count+field_3_4_count+field_3_5_count+0.0001)*100.,1))+'%)',color='w')
  text(19.5,6.0,'N = '+str(field_3_4_count)+'\n('+str(round(field_3_4_count/float(field_3_1_count+field_3_2_count+field_3_3_count+field_3_4_count+field_3_5_count+0.0001)*100.,1))+'%)',color='k')
  text(25.5,6.0,'N = '+str(field_3_5_count)+'\n('+str(round(field_3_5_count/float(field_3_1_count+field_3_2_count+field_3_3_count+field_3_4_count+field_3_5_count+0.0001)*100.,1))+'%)',color='w')

  text(1.5,10.0,'N = '+str(field_2_1_count)+'\n('+str(round(field_2_1_count/float(field_2_1_count+field_2_2_count+field_2_3_count+field_2_4_count+field_2_5_count+0.0001)*100.,1))+'%)',color='w')
  text(7.5,10.0,'N = '+str(field_2_2_count)+'\n('+str(round(field_2_2_count/float(field_2_1_count+field_2_2_count+field_2_3_count+field_2_4_count+field_2_5_count+0.0001)*100.,1))+'%)',color='w')
  text(13.5,10.0,'N = '+str(field_2_3_count)+'\n('+str(round(field_2_3_count/float(field_2_1_count+field_2_2_count+field_2_3_count+field_2_4_count+field_2_5_count+0.0001)*100.,1))+'%)',color='k')
  text(19.5,10.0,'N = '+str(field_2_4_count)+'\n('+str(round(field_2_4_count/float(field_2_1_count+field_2_2_count+field_2_3_count+field_2_4_count+field_2_5_count+0.0001)*100.,1))+'%)',color='w')
  text(25.5,10.0,'N = '+str(field_2_5_count)+'\n('+str(round(field_2_5_count/float(field_2_1_count+field_2_2_count+field_2_3_count+field_2_4_count+field_2_5_count+0.0001)*100.,1))+'%)',color='w')

  text(1.5,14.0,'N = '+str(field_1_1_count)+'\n('+str(round(field_1_1_count/float(field_1_1_count+field_1_2_count+field_1_3_count+field_1_4_count+field_1_5_count+0.0001)*100.,1))+'%)',color='w')
  text(7.5,14.0,'N = '+str(field_1_2_count)+'\n('+str(round(field_1_2_count/float(field_1_1_count+field_1_2_count+field_1_3_count+field_1_4_count+field_1_5_count+0.0001)*100.,1))+'%)',color='k')
  text(13.5,14.0,'N = '+str(field_1_3_count)+'\n('+str(round(field_1_3_count/float(field_1_1_count+field_1_2_count+field_1_3_count+field_1_4_count+field_1_5_count+0.0001)*100.,1))+'%)',color='w')
  text(19.5,14.0,'N = '+str(field_1_4_count)+'\n('+str(round(field_1_4_count/float(field_1_1_count+field_1_2_count+field_1_3_count+field_1_4_count+field_1_5_count+0.0001)*100.,1))+'%)',color='w')
  text(25.5,14.0,'N = '+str(field_1_5_count)+'\n('+str(round(field_1_5_count/float(field_1_1_count+field_1_2_count+field_1_3_count+field_1_4_count+field_1_5_count+0.0001)*100.,1))+'%)',color='w')

  text(1.5,18.0,'N = '+str(field_0_1_count)+'\n('+str(round(field_0_1_count/float(field_0_1_count+field_0_2_count+field_0_3_count+field_0_4_count+field_0_5_count+0.0001)*100.,1))+'%)',color='k')
  text(7.5,18.0,'N = '+str(field_0_2_count)+'\n('+str(round(field_0_2_count/float(field_0_1_count+field_0_2_count+field_0_3_count+field_0_4_count+field_0_5_count+0.0001)*100.,1))+'%)',color='w')
  text(13.5,18.0,'N = '+str(field_0_3_count)+'\n('+str(round(field_0_3_count/float(field_0_1_count+field_0_2_count+field_0_3_count+field_0_4_count+field_0_5_count+0.0001)*100.,1))+'%)',color='w')
  text(19.5,18.0,'N = '+str(field_0_4_count)+'\n('+str(round(field_0_4_count/float(field_0_1_count+field_0_2_count+field_0_3_count+field_0_4_count+field_0_5_count+0.0001)*100.,1))+'%)',color='w')
  text(25.5,18.0,'N = '+str(field_0_5_count)+'\n('+str(round(field_0_5_count/float(field_0_1_count+field_0_2_count+field_0_3_count+field_0_4_count+field_0_5_count+0.0001)*100.,1))+'%)',color='w')

  text(1.5,1.0,'$\sigma$ = '+str(round(array(field_4_1_delta).std(),3)),color='w')
  text(7.5,1.0,'$\sigma$ = '+str(round(array(field_4_2_delta).std(),3)),color='w')
  text(13.5,1.0,'$\sigma$ = '+str(round(array(field_4_3_delta).std(),3)),color='w')
  text(19.5,1.0,'$\sigma$ = '+str(round(array(field_4_4_delta).std(),3)),color='w')

  text(1.5,5.0,'$\sigma$ = '+str(round(array(field_3_1_delta).std(),3)),color='w')
  text(7.5,5.0,'$\sigma$ = '+str(round(array(field_3_2_delta).std(),3)),color='w')
  text(13.5,5.0,'$\sigma$ = '+str(round(array(field_3_3_delta).std(),3)),color='w')
  text(19.5,5.0,'$\sigma$ = '+str(round(array(field_3_4_delta).std(),3)),color='k')

  text(1.5,9.0,'$\sigma$ = '+str(round(array(field_2_1_delta).std(),3)),color='w')
  text(7.5,9.0,'$\sigma$ = '+str(round(array(field_2_2_delta).std(),3)),color='w')
  text(13.5,9.0,'$\sigma$ = '+str(round(array(field_2_3_delta).std(),3)),color='k')
  text(19.5,9.0,'$\sigma$ = '+str(round(array(field_2_4_delta).std(),3)),color='w')

  text(1.5,13.0,'$\sigma$ = '+str(round(array(field_1_1_delta).std(),3)),color='w')
  text(7.5,13.0,'$\sigma$ = '+str(round(array(field_1_2_delta).std(),3)),color='k')
  text(13.5,13.0,'$\sigma$ = '+str(round(array(field_1_3_delta).std(),3)),color='w')
  text(19.5,13.0,'$\sigma$ = '+str(round(array(field_1_4_delta).std(),3)),color='w')

  text(1.5,17.0,'$\sigma$ = '+str(round(array(field_0_1_delta).std(),3)),color='k')
  text(7.5,17.0,'$\sigma$ = '+str(round(array(field_0_2_delta).std(),3)),color='w')
  text(13.5,17.0,'$\sigma$ = '+str(round(array(field_0_3_delta).std(),3)),color='w')
  text(19.5,17.0,'$\sigma$ = '+str(round(array(field_0_4_delta).std(),3)),color='w')

  savefig(path+'/'+name+'.pdf')


def quantile(x, q,  qtype = 7, issorted = False):
  """
Args:
x - input data
q - quantile
qtype - algorithm
issorted- True if x already sorted.
 
Compute quantiles from input array x given q.For median,
specify q=0.5.
 
References:
http://reference.wolfram.com/mathematica/ref/Quantile.html
http://wiki.r-project.org/rwiki/doku.php?id=rdoc:stats:quantile
 
Author:
Ernesto P.Adorio Ph.D.
UP Extension Program in Pampanga, Clark Field.
  """
  if not issorted:
    y = sorted(x)
  else:
    y = x
  if not (1 <= qtype <= 9):
    return None  # error!
    # Parameters for the Hyndman and Fan algorithm
  abcd = [(0,   0, 1, 0), # inverse empirical distrib.function., R type 1
         (0.5, 0, 1, 0), # similar to type 1, averaged, R type 2
         (0.5, 0, 0, 0), # nearest order statistic,(SAS) R type 3
         (0,   0, 0, 1), # California linear interpolation, R type 4
         (0.5, 0, 0, 1), # hydrologists method, R type 5
         (0,   1, 0, 1), # mean-based estimate(Weibull method), (SPSS,Minitab), type 6 
         (1,  -1, 0, 1), # mode-based method,(S, S-Plus), R type 7
         (1.0/3, 1.0/3, 0, 1), # median-unbiased ,  R type 8
         (3/8.0, 0.25, 0, 1)   # normal-unbiased, R type 9.
         ]

  a, b, c, d = abcd[qtype-1]
  n = len(x)
  g, j = modf( a + (n+b) * q -1)
  if j < 0:
    return y[0]
  elif j > n:
    return y[n]

  j = int(floor(j))
  if g ==  0:
    return y[j]
  else:
    return y[j] + (y[j+1]- y[j])* (c + d * g)

#needs adaption!!
def giant_phasecut(cat,stacorrfile,leng=120.,cnv=False,archive='/RAID/tipage/data'):
  """
  function for cutting out phases from a catalogue, creating a GIANT phase-file containing the actual pick to go with the cut data
  raytracing for stations not picked, listed as weight 4 picks
  """

  #if cnv:
  #  cat = read_cnv(catalogue,pol=True,statlen=4)

  #else:
  #  cat = read_catalogue(catalogue)

  # go through dict, cut out arrivals and write phase files
  for i in range(1,len(cat)+1):
    print i
    otime = cat[i]['origin_time']
    foldername = '%02d%02d%02d_%02d%02d%02d' % ((otime.year-2000),otime.month,otime.day,otime.hour,otime.minute,otime.second)
    os.system('mkdir '+foldername)
    phasefile = open(foldername+'.pha','w')
    ot_format = '%2d-%02d-%02d %02d:%02d:%5.2f' % ((otime.year-2000),otime.month,otime.day,otime.hour,otime.minute,(otime.second+(otime.microsecond/1e6)))
    phasefile.write('%22s %7.3fN %7.3fE %6.2f  0.0\n' % (ot_format,float(cat[i]['event_lat']),float(cat[i]['event_lon']),float(cat[i]['event_depth']))) #header line

    start = cat[i]['origin_time']
    stop =   start + leng 
    jday = '%03d' % (otime.julday)
    cuttime = '%02d%02d%02d' % (otime.hour,otime.minute,otime.second)

    #cut the P picks

    statlist = []
    for file in glob.glob(archive+'/'+str(otime.year)+'/*/*/??Z.D/*'+str(jday)):
      statlist.append(os.path.basename(file).split('.')[0])
    for stat in cat[i]['Picks'].keys():
      file = glob.glob(archive+'/'+str(otime.year)+'/*/'+stat+'/??Z.D/*.'+jday)[0]
      trace = read(file,format='MSEED',starttime=start,endtime=stop,nearest_sample=True)

      try:
        trace.write(foldername+'/'+stat+'.'+trace[0].stats.channel[2]+'.'+str(otime.year)[2:]+'.'+jday+'.'+str(cuttime),'MSEED') #trying to keep name short...

        phasefile.write('0   %4s        P I %1s %1d %08.2f -1.0 -1.0\n' % (stat,cat[i]['Picks'][stat]['P-pol'],int(cat[i]['Picks'][stat]['P-weight']),float(cat[i]['Picks'][stat]['P-pick'])))
      except:
        pass
      statlist.remove(stat)
  # raytracing for missing stations, cutting and phase file writing...

    #raytracing
    nlloc.initialize('/home/sippl/sandbox/nlloc')
    nlloc.change_input_line('/home/sippl/sandbox/nlloc',(float(cat[i]['event_lon'])-67.)*111.2*cos(float(cat[i]['event_lat'])*2*pi/360.),(float(cat[i]['event_lat'])-33.)*111.2,float(cat[i]['event_depth']))
    nlloc.call_nlloc('/home/sippl/sandbox/nlloc')
    raydict = nlloc.read_output('/home/sippl/sandbox/nlloc')

    #missing: add station corrections!!

    stacorrdict = simulps.read_stacorr(stacorrfile)

    for sta in statlist:
      if float(raydict[sta]) < 60.:
        file = glob.glob(archive+'/'+str(otime.year)+'/*/'+sta+'/??Z.D/*.'+jday)[0]
        trace = read(file,format='MSEED',starttime=start,endtime=stop,nearest_sample=True)
        try:
          trace.write(foldername+'/'+sta+'.'+trace[0].stats.channel[2]+'.'+str(otime.year)[2:]+'.'+jday+'.'+str(cuttime),'MSEED') #trying to keep name short...
          phasefile.write('0   %4s        P I %1s %1d %08.2f -1.0 -1.0\n' % (sta,'_',4,float(raydict[sta])+float(stacorrdict[sta]['corr_P'])))
        except:
          continue

    phasefile.close()

#needs checking!!
#input catalog can already have Picks, info, etc...then it's just filled further with Mw, centroid_depth, %DC etc.
def read_mts(cat,infile): 
  """
  read mts summary file into mechdict format (same as produced by HASH.read_hashout())
  """
  infi = open(infile,'r')
  data = infi.readlines()
  infi.close()

  for line in data:
    if line[0] == '#':
      continue #eliminate header(s)
    else:
      datestring,timestring,lat,lon,dep,mb,Mw,strike,dip,rake,strike2,dip2,rake2,pax_val,pax_strike,pax_dip,tax_val,tax_strike,tax_dip,dum,Mxx,Myy,Mzz,Mxy,Mxz,Myz,Nstat = line.strip('\n').split(None)
      orig_time = UTCDateTime(2000+int(datestring.split('/')[0]),int(datestring.split('/')[1]),int(datestring.split('/')[2]),int(timestring.split(':')[0]),int(timestring.split(':')[1]),int(timestring.split(':')[2]))
      for ev in cat.events.keys():
        if abs(cat.events[ev].info['orig_time'] - orig_time) < 5.: #arbitrary...
          cat.events[ev].info['centroid_depth'] = float(dep)
          cat.events[ev].info['DC%'] = int(dc)    
          cat.events[ev].info['pax'] = (int(pax_strike),int(pax_dip))
          cat.events[ev].info['tax'] = (int(tax_strike),int(tax_dip))
          cat.events[ev].info['Mw'] = float(Mw)
          cat.events[ev].info['strike/dip/rake'] = ((int(strike),int(dip),int(rake)),(int(strike2),int(dip2),int(rake2)))
          cat.events[ev].info['MomentTensor'] = array([[Mxx,Mxy,Mxz],[Mxy,Myy,Myz],[Mxz,Myz,Mzz]])
          break

  return cat    

#obsolete
def read_alllocs():
  """
  read events from summary format into (largely simplified) dictionary format
  """
  location='/home/sippl/sandbox/hypoDD/all_catalog/make_figures/all_locs'

  infile = open(location,'r')
  data = infile.readlines()
  infile.close()
  outdict = {}
  n = 1

  for line in data:
    outdict[n] = {}
    yr,mon,day,timestr,lon,lat,dep,mag,ident = line.strip('\n').split(None)
    outdict[n]['event_lat'] = float(lat)
    outdict[n]['event_lon'] = float(lon)
    outdict[n]['event_depth'] = float(dep)
    hr,min,sec = timestr.split(':')
    outdict[n]['origin_time'] = UTCDateTime(int(yr),int(mon),int(day),int(hr),int(min),float(sec))
    outdict[n]['magnitude'] = float(mag)
    outdict[n]['identifier'] = ident

    n += 1

  return outdict


def find_bin_val(num,bin_width,mirror=False):
  """
  helper function for binning of rose diagram(s)
  """
  bin_centres = arange(bin_width/2.,360.,bin_width)
  for i in bin_centres:
    if abs(num - i) <= bin_width/2.:
      num_new = i
      break

  if mirror:
    if num_new > 90 and num_new <= 270:
      num_new = (num_new+180)%360

  return num_new

def rose_diagram(mechdict,figure_inst='new',evalstr='tax_strike1',bin_width=20,lat_max=90.0,lat_min=-90.0,lon_max=180.0,lon_min=-180.,dep_max=750.,dep_min=0.,mag_max=10.0,mag_min=0.0,qualitylist=['A','B','C','D'],color='blue',output='None',norm=False,mirror=False):
  """
  evaluate earthquake mechanisms in rose diagrams for angles (strike,dip,rake,P-,T-axes)
  selection of events by location (lat,lon,dep) and magnitude possible (min/max values)
  give qualitylist to only evaluate certain quality classes (standard: all)
  evalstr defines which mechdict attribute is evaluated
  choose figure_inst = 'new' to open a new figure, figure_inst = 'old' to plot into an existing figure (color handling via kwarg color)
  if outfile is different from 'None', a pdf output will be generated (Name equal to output string)
  mirror - choose True for rose diagram only covering upper hemisphere (lower hemisphere values will receive += 180 deg)
  """
  #import from HASH --> read_hashout

  vallist = []

  for entry in mechdict.keys():
    if mechdict[entry]['event_lat'] < lat_max and mechdict[entry]['event_lat'] > lat_min and mechdict[entry]['event_lon'] > lon_min and mechdict[entry]['event_lon'] < lon_max and mechdict[entry]['event_depth'] < dep_max and mechdict[entry]['event_depth'] > dep_min:
      if mechdict[entry]['quality1'] in qualitylist:
        vallist.append(find_bin_val(mechdict[entry][evalstr],bin_width,mirror=mirror))

  if figure_inst == 'new':
    fig = figure(figsize=(8,8))
    ax = fig.add_axes([0.1,0.1,0.8,0.8],polar=True)

  #binning of vallist necessary 

  width= bin_width*pi/180.

  nbins = round(360/bin_width)

  for j in arange(nbins):
    num = 0
    for vall in vallist:
      if vall == ((bin_width/2.) + bin_width*j):
        num += 1
    if norm:
      num /= float(len(vallist))
    if figure_inst == 'new':
      ax.bar((360-(bin_width*j + bin_width)+90)*pi/180.,num,width=width,bottom=0.0,edgecolor=color,facecolor='none')
    elif figure_inst == 'old':
      bar((360-(bin_width*j + bin_width)+90)*pi/180.,num,width=width,bottom=0.0,edgecolor=color,facecolor='none')


  xticks([(90.*2.*pi)/360.,(45.*2.*pi)/360.,0.,(315.*2.*pi)/360.,(270.*2.*pi)/360.,(225.*2*pi)/360.,(180.*2*pi)/360.,(135.*2*pi)/360.],['N','45 deg','E','135 deg','S','225 deg','W','315 deg'])    
  if figure_inst == 'new':
    show()

  if not output == 'None':
    savefig(output+'.pdf')

def read_ermapper_grid(ifile):
    ers_index = ifile.find('.ers')
    if ers_index > 0:
        data_file = ifile[0:ers_index]
        header_file = ifile
    else:
        data_file = ifile
        header_file = ifile + '.ers'
        
    header = read_ermapper_header(header_file)

    nroflines = int(header['nroflines'])
    nrofcellsperlines = int(header['nrofcellsperline'])
    data = read_ermapper_data(data_file)
    data = reshape(data,(nroflines,nrofcellsperlines))
    return data

"""
next three functions taken from https://anuga.anu.edu.au/svn/anuga/branches/anuga_1_1/branches/numpy/anuga/abstract_2d_finite_volumes/ermapper_grids.py
"""

def read_ermapper_header(ifile):
    # function for reading an ERMapper header from file
    header = {}

    fid = open(ifile,'rt')
    header_string = fid.readlines()
    fid.close()

    for line in header_string:
        if line.find('=') > 0:
            tmp_string = line.strip().split('=')
            header[tmp_string[0].strip().lower()]= tmp_string[1].strip()

    return header

def read_ermapper_data(ifile, data_format = float32):
    # open input file in a binary format and read the input string
    fid = open(ifile,'rb')
    input_string = fid.read()
    fid.close()

    # convert input string to required format (Note default format is num.float32)
    grid_as_float = fromstring(input_string,data_format)
    return grid_as_float

#from _numeric import delandaz2coord

#def delandaz2coord(d, az, lat0, lon0):
#   """
#    -> lat, lon
#
#   Returns the coordinates of the point which is at an azimuth of 'az'
#   and a distance of 'd' as seen from the point (lat0,lon0).
#   """
#   d = d/111.195
#   # x in km, az,lat,lon in deg
#   _lat, az, d = pi*(90.-lat0)/180., az*pi/180., d*pi/180.
#   _latx = arccos(cos(_lat)*cos(d)+sin(_lat)*sin(d)*cos(az))
#   _lonx = arcsin(sin(d)*sin(az)/sin(_latx))
#   return 90. - 180.*_latx/pi, lon0 + 180.*_lonx/pi

def plot_recordsection(datapath,info_file='/home/sippl/info_file',eq_loc=(-30.,125.,35.),eq_time=UTCDateTime(2000,1,1,0,0,0),stations=[],component='Z',window=[30.,60.]):
  """
  plot record section for single earthquake, based on calculated traveltimes (ak135) and distances from event
  eq_loc - tuple for EQ's location (lat,lon,depth)
  eq_time - EQ's origin time (UTCDateTime object)
  stations - list of station names that should be included (provided there is data). Names have to appear in info_file
  """
  #event: retrieve station locations, calculate distances and theoretical travel times
  #example 1: Kalgoorlie; 30.705S, 121.192E, 1km depth; 26/2/2014, 00:00:17
  #example 2: Nepal; 28.147N, 84.708E, 15km depth; 25/4/2015; 06:11:26
  from obspy.taup import TauPyModel

  model = TauPyModel(model='ak135')

  ttdict = {}

  for stat in stations:
    #coordinate lookup
    #distance calculation
    #travel time calculation
    #cut data
    sta_lat,sta_lon,sta_elev = get_stat_coords(stat,info_file)

    #get P traveltime     
    dist_m,az,baz = obspy.core.util.geodetics.gps2DistAzimuth(sta_lat,sta_lon,eq_loc[0],eq_loc[1])
    dist_deg = dist_m / (1000. * 111.195)
    arr_P = model.get_travel_times(eq_loc[2],dist_deg,phase_list=['P'])[0].time 
   
    #define cut window (30 s before to 60 s after)
    bef = eq_time + arr_P - window[0]
    after = eq_time + arr_P + window[1]
    arr = eq_time + arr_P
    
    #load the correct trace(s)
    trcstrg = '%4s%02i%02i%02i' % (stat,(eq_time.year - 2000.),eq_time.month,eq_time.day)
    
    if eq_time.hour == 23:
      day_after = eq_time + 86400.
      trcstrg2 = '%4s%02i%02i%02i' % (stat,(day_after.year - 2000.),day_after.month,day_after.day)
      try:
        dat_str = read(glob.glob(datapath+'/'+stat+'/'+trcstrg+'*'+component)[0])
        dat_str += read(glob.glob(datapath+'/'+stat+'/'+trcstrg2+'*'+component)[0])
      except IndexError:
        continue
      #cat two traces, then cut
    elif eq_time.hour == 0:
      day_before = eq_time - 86400.
      trcstrg2 = '%4s%02i%02i%02i' % (stat,(day_before.year - 2000.),day_before.month,day_before.day)
      try:
        dat_str = read(glob.glob(datapath+'/'+stat+'/'+trcstrg+'*'+component)[0])
        dat_str += read(glob.glob(datapath+'/'+stat+'/'+trcstrg2+'*'+component)[0])      
      except IndexError:
        continue

    else:
      try:
        dat_str = read(glob.glob(datapath+'/'+stat+'/'+trcstrg+'*'+component)[0])
      except IndexError:
        continue

    #cut out data (if not empty, create entry in ttdict)
    dat_cut = dat_str.slice(bef,after)
    dat_cut.detrend('demean')
    dat_cut.detrend('linear')
    if not len(dat_cut) == 0:
      ttdict[stat] = {}
      ttdict[stat]['ttime'] = arr_P
      ttdict[stat]['arrival'] = arr
      ttdict[stat]['distance'] = dist_deg
      ttdict[stat]['trace'] = dat_cut
 
  #s, not samples
  figure(figsize=(8,25))
  maxpos = -9999.
  minpos = 9999.
  for j in ttdict.keys():
    ttdict[j]['trace'].merge(method=-1)
    ar = (ttdict[j]['trace'][0].data).astype('float')
    ar /= ar.max()
    xvals = arange(len(ar)) / ttdict[j]['trace'][0].stats.sampling_rate
    plot(xvals,ar - ttdict[j]['ttime'],'k')
    if (ttdict[j]['ttime'] * (-1)) >= maxpos:
      maxpos = (-1) * ttdict[j]['ttime']
    if (ttdict[j]['ttime'] * (-1) ) <= minpos:
      minpos = (-1) * ttdict[j]['ttime'] 

    text(0,(-1)*ttdict[j]['ttime'],j,color='red')

  plot([window[0],window[0]], [minpos-1,maxpos+1],'r--') 
  savefig('blap.pdf')
  xlabel('time [s]')
  ylabel('travel time [s*(-1)]')
  ylim([minpos-1,maxpos+1])
  xlim([0,xvals.max()])

  return ttdict

def write_evlines(cat,outfile):
  """
  writes event lines from catalog input (useful for creating additional MPX input after catalog.compare)
  """
  otf = open(outfile,'w')

  for j in sorted(cat.events.keys()):
    org = cat.events[j].info['origin_time']
    lat = cat.events[j].info['event_lat']
    lon = cat.events[j].info['event_lon']
    dep = cat.events[j].info['event_dep']
    #write dummies for nobs, gap and rms
 
    line = '%13.2f %04i %02i %02i %02i %02i %06.3f %8.4f %9.4f %7.2f 10 222 0.75\n' % (org.timestamp,org.year,org.month,org.day,org.hour,org.minute,(org.second + org.microsecond/1e6),lat,lon,(-1)*dep)
    otf.write(line)
  otf.write('\n')

def import_locking(infile):
  """

  """
  infl = open(infile,'r')
  dat =infl.readlines()
  infl.close()

  latlist = []
  lonlist = []
  locklist = []

  for j in dat:
    lon,lat,lock = j.split(None)
    if lock == 'NaN':
      continue
    latlist.append(float(lat))
    lonlist.append(float(lon))
    locklist.append(float(lock))


  Lats = (unique(sorted(latlist))).tolist()
  Lons = (unique(sorted(lonlist))).tolist()

  Locks = zeros([len(Lats),len(Lons)])

  for l in range(len(locklist)):
    la = latlist[l]
    lo = lonlist[l]
    pos_lon = Lons.index(lo)
    pos_lat = Lats.index(la)
    Locks[pos_lat,pos_lon] = locklist[l]


  return Lats,Lons,Locks

def import_slab(model='slab1.0'):
  """
  model - either 'slab1.0' or 'tassara' or 'ipoc'
  import dlab models (either SLAB1.0 or Tassara and Echaurren (2013)) into python format, can then be used for cross section plotting or for filtering Catalog by distance to slab surface
  """
  #get input file
  if model == 'slab1.0':
    infile = open('/home/sippl/jacob/kur_slab1.0_clip.xyz','r')
    dat = infile.readlines()
    infile.close()
  elif model == 'tassara':
    infile = open('/home/sippl/grids/slab_tasa.xyz','r')
    dat = infile.readlines()
    infile.close()
  elif model == 'ipoc':
    infile = open('/home/sippl/grids/IPOC_slab_grid.txt','r')
    dat = infile.readlines()
    infile.close()
  elif model == 'etopo':
    infile = open('/home/sippl/grids/etopo_NChile.xyz','r')
    dat = infile.readlines()
    infile.close()
  else: #model path is given
    infile = open(model,'r')
    dat = infile.readlines()
    infile.close()



  #save as matrix (adjust lon version of slab1.0)
  latlist = []
  lonlist = []
  deplist = []

  for j in dat:
    lon,lat,dep = j.split(None)
    if dep == 'NaN': #no need to store NaNs...
      continue
    latlist.append(float(lat))
    deplist.append(float(dep))
    if float(lon) < 180.:
      lonlist.append(float(lon))
    else:
      lonlist.append(float(lon) - 360.)


  #now reduce lat and lon to unique values
  Lats = (unique(sorted(latlist))).tolist()
  Lons = (unique(sorted(lonlist))).tolist() 

  Deps = zeros([len(Lats),len(Lons)])

  for k in range(len(deplist)):
    la = latlist[k]
    lo = lonlist[k]

    pos_lon = Lons.index(lo)
    pos_lat = Lats.index(la)

    Deps[pos_lat,pos_lon] = deplist[k]


  print model
  print Deps

  return Lats,Lons,Deps #OK


def get_closest_profile(lat,Lats,Lons,Deps):
  """
  extract closest E-W profile to given latitude from slab model
  """
  #find closest latitude to lat
  dist = 999.
  for j in range(len(Lats)):
    if abs(lat - Lats[j]) < dist:
      dist = abs(lat - Lats[j])
      pos = j

  #now extract that line from the matrix
  prof = Deps[pos,:]

  return prof


def compute_point_distance(lon,dep,prof,Lons):
  '''
  lon,dep - point (EQ) coordinates - lat omitted because structure is, to first order, 2.5D
  prof, Lons - slab coordinates along closest profile (output of get_closest_profile() )
 
  '''
  import shapely.geometry as geom

  #build matrix from prof, Lons
  matrx = zeros([len(prof),2])
  matrx[:,0] = (array(Lons) * 111.195 * cos(21*pi/180.))
  matrx[:,1] = prof

  line = geom.LineString(matrx)
  lon = lon * 111.195 * cos(21*pi/180.)
  point = geom.Point(lon,dep*(-1))
  #get sign (is it above or below?)
  lon_true = line.project(point)
  z_true = (line.interpolate(lon_true)).coords[0][1]
  if (-1)*z_true < dep:
    return point.distance(line) * (-1)  
  else:
    return point.distance(line)


def plot_topoprof(centr_lat,swath=0,step=1.):
  """
  swath - swath width in km
  step - only makes sense if swath != 0; step length of profiles for swath calculation 

 
  """
  figure(figsize=(12,3.5))
  if swath == 0:
    Lats, Lons, Deps = import_slab(model='etopo')
    prof = get_closest_profile(centr_lat,Lats,Lons,Deps)
    fac = cos(centr_lat * pi / 180.) * 111.195         
  else:
    # get central csec
    Lats, Lons, Deps = import_slab(model='etopo')
    prof = array(get_closest_profile(centr_lat,Lats,Lons,Deps))
    fac = cos(centr_lat * pi / 180.) * 111.195    

    #then loop over others, count how many (for normalization)
    s_min = centr_lat - (swath * (1/fac))
    s_max = centr_lat + (swath * (1/fac))
    lats = arange(s_min,s_max,(step * (1/fac)))  
    count = 1  
    for lat in lats:
      a = array(get_closest_profile(lat,Lats,Lons,Deps))
      prof += a
      count += 1 

    #normalize
    prof /= float(count)

  mnlon = -71.5

  plot((array(Lons) - mnlon)*fac,prof,'k-')
  plot((array(Lons) - mnlon)*fac,zeros(len(prof)),'k--')
  xlim([100,500])
  
def smooth_catalog(cat,grdspac=0.05,mom=True,cor=0.0,rad=0.25):
  """
  routine for spatial sampling of earthquake catalog 
  will return a grid, the value at each grid point is the sum of earthquake numbers (if mom=False) or seismic moment (if mom=True) of a circle with radius rad around the grid point
  cat - Catalog object
  grdspac - desired grid spacing in deg 
  mom - False -> earthquake numbers are summed up; True -> moments are summed
  cor - magnitude correction (for turning ML to Mw)
  rad - radius (in deg) around each grid point over which summation is applied
  """

  #get grid boundaries (2D grid only)
  lats = []
  lons = []
  for ev in cat.events.keys():
    lats.append(cat.events[ev].info['event_lat'])
    lons.append(cat.events[ev].info['event_lon'])

  mxlon = array(lons).max()
  mnlon = array(lons).min()
  mxlat = array(lats).max()
  mnlat = array(lats).min()
 
  #set up grid
  xs = arange(mnlon,mxlon,grdspac)
  ys = arange(mnlat,mxlat,grdspac)
  X,Y = meshgrid(xs,ys)

  matrx = zeros(shape(X))

  #loop: for each EQ, find nodes it contributes to, sum up value (convert to mom if necessary)
  for k in sorted(cat.events.keys()):
    for l in range(len(xs)):
      if abs(xs[l] - cat.events[k].info['event_lon']) <= rad:
        for m in range(len(ys)):
          if sqrt(abs(ys[m] - cat.events[k].info['event_lat'])**2 + abs(xs[l] - cat.events[k].info['event_lon'])**2 ) <= rad:
            #add to matrx position
            #inverse distance weighting
            weight = (rad - sqrt(abs(ys[m] - cat.events[k].info['event_lat'])**2 + abs(xs[l] - cat.events[k].info['event_lon'])**2 )) / rad
            if not mom: 
              matrx[m][l] += weight
            else:
              moment = 10**(1.5*(cat.events[k].info['magnitude'] + cor) + 9.1)
              matrx[m][l] += moment*(weight**5)


  return xs,ys,matrx


def read_JMA(infile):
  """
  read JMA Hypocenter format into Catalog object
  """
  ifl = open(infile,'r')
  data = ifl.readlines()
  ifl.close()

  ev_num = 1
  cat = basic.Catalog()
  for j in data:
    #origin time
    yr = j[1:5]
    mon = j[5:7]
    dy = j[7:9]
    hr = j[9:11]
    mn = j[11:13]
    sc = j[13:15]
    ms = j[15:17]
    print sc,ms
    try:
      sec = float(sc) + (float(ms)/100.)
    except ValueError:
      continue
    orig_time = UTCDateTime(int(yr),int(mon),int(dy),int(hr),int(mn),sec)
    #hypocenter location
    lat_deg = j[21:24]
    lat_min1 = j[24:26]
    lat_min2 = j[26:28]
    lat = int(lat_deg) + (int(lat_min1) + (float(lat_min2)/100.))/60.
    lon_deg = j[32:36]
    lon_min1 = j[36:38]
    lon_min2 = j[38:40]
    lon = int(lon_deg) + (int(lon_min1) + (float(lon_min2)/100.))/60.
 
    try:
      dep1 = int(j[44:47])
      dep2 = int(j[47:49])
      dep = dep1 + (float(dep2)/100.)
    except ValueError:
      continue
    #define event
    ev = basic.Event(orig_time,lat,lon,dep)
    ev.Picks = {}
    
    #magnitude
    mag1 = j[52]
    mag2 = j[53]
    try:
      mag = int(mag1) + 0.1*int(mag2)
    except ValueError:
      mag = -99.
    ev.info['magnitude'] = mag
    if ev.info['event_lat'] > 33. and ev.info['event_lat'] < 45. and ev.info['event_dep'] > 50.:
      cat.events[orig_time.isoformat()] = cp(ev)
      ev_num+=1
  cat.info['Nev'] = ev_num
  cat.info['source'] = infile

  return cat

def read_srtm(infile):
  """
  read srtm ASCII file and retun masked array of altitudes that can be plotted with basemap
  """
  ifile = open(infile,'r')
  data = ifile.readlines()
  ifile.close()

  sze = int(data[0].split(None)[1])
  lon_lleft = float(data[2].split(None)[1])
  lat_lleft = float(data[3].split(None)[1])
  inc = float(data[4].split(None)[1])

  matrx = zeros([sze,sze])
  xar = linspace(lon_lleft,(lon_lleft + inc*sze),num=sze)
  yar = linspace(lat_lleft,(lat_lleft + inc*sze),num=sze)

  count = 0
  for k in data[6:]:
    lst = k.split(None)
    matrx[count] = lst
    count += 1

  matrx_m = ma.masked_where(matrx == -9999,matrx)
  return xar,yar,matrx_m

