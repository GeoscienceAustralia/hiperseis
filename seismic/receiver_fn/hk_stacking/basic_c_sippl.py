#module basic
#coding=utf-8

from numpy import copy
from glob import glob
import os
from obspy.core import UTCDateTime
from copy import deepcopy as cp
from pylab import *
import util

def read_info(info_path,header=False):
  '''
  reads info_file, turns it into dictionary format
  '''
  infodict = {}

  fl = open(info_path,'r')
  dat = fl.readlines()
  fl.close()

  if header: 
    dat = dat[1:] #do not read header line

  for line in dat:
    sta,seis,sta3,seis_ser,dig,dig_ser,dum,lat,lon,elev,pregain,fsamp,name,net,resp_E,resp_N,resp_Z,path,pref = line.strip('\n').split(None)
    if sta[0] == '_': #station with seismometer change; go back to last line without _, get station name and archive from there; write different configurations
      sta_new = sta[1:]
      if not sta_new in infodict.keys():
        print 'Warning: station name beginning with _ only allowed if indicating a configuration change to a peviously (!) defined station - station ignored!'
        continue
      if not type(infodict[sta_new]['responses']) == dict:
        tmp = tuple(copy(infodict[sta_new]['responses']))
        infodict[sta_new]['responses'] = {}
        infodict[sta_new]['responses']['config1'] = tmp
        infodict[sta_new]['responses']['switch1_2'] = path
        infodict[sta_new]['responses']['config2'] = (pregain,resp_E,resp_N,resp_Z) 
        tmp2 = tuple(copy(infodict[sta_new]['instruments']))
        infodict[sta_new]['instruments'] = {}          
        infodict[sta_new]['instruments']['config1'] = tmp2
        infodict[sta_new]['instruments']['config2'] = (seis,dig)
      else: #there already is more than one config
        num = 1 + ((len(infodict[sta_new]['responses']) + 1) / 2)
        infodict[sta_new]['responses']['config'+str(num)] = (pregain,resp_E,resp_N,resp_Z)  
        infodict[sta_new]['instruments']['config'+str(num)] = (seis,dig)
        infodict[sta_new]['responses']['switch'+str(num-1)+'_'+str(num)] = path
      continue
    if not infodict.has_key(sta):
      infodict[sta] = {}
    else:
      sta_before = sta + '_' + infodict[sta]['network']
      infodict[sta_before] = infodict[sta].copy()
      del infodict[sta]
      sta = sta+'_'+net
      infodict[sta] = {}
    infodict[sta]['lat'] = lat
    infodict[sta]['lon'] = lon
    infodict[sta]['elev'] = elev
    infodict[sta]['network'] = net
    infodict[sta]['responses'] = (pregain,resp_E,resp_N,resp_Z)
    infodict[sta]['fsamp'] = fsamp
    infodict[sta]['abbrev'] = sta3
    infodict[sta]['archive'] = path
    infodict[sta]['instruments'] = (seis,dig)
    infodict[sta]['name'] = sta
    infodict[sta]['locality'] = name
    infodict[sta]['pref'] = pref

  return infodict  

def enter_4names(infodict):
  '''
  adds unique 4-letter names to infodict
  '''
  list4 = []
  for j in infodict.keys():
    if len(j) <= 4:
      infodict[j]['name4'] = j
    else:
      #truncation
      name_trunc = j[0:4]
      #now check whether that already exists
      if not name_trunc in infodict.keys() and not name_trunc in list4:
        list4.append(name_trunc)
        infodict[j]['name4'] = name_trunc
      else:
        count = 1
        while count < 10:
          name_new = j[0:3] + str(count)
          if not name_new in infodict.keys() and not name_new in list4:
            list4.append(name_new)
            infodict[j]['name4'] = name_new
            count = 99
          else:
            count += 1
         
  return infodict     


def find_resp(stat,resps):
  """
  extracts response information in obspy's dictionary format for given station + channel using information stored in ext.py and info_file
  stat: to be given as 2-to-5-digit string (has to be identical to station name contained in info_file)
  """
  from resp import *

  if resps[1] == resps[2] == resps[3] == 'XXX': #dummy entered
    #print 'Warning: no valid response information available for station '+stat
    return {},{},{}

  resp_E = {}
  resp_E['poles'] = eval(resps[1]).p
  resp_E['zeros'] = eval(resps[1]).z[1:]
  resp_E['gain'] = eval(resps[1]).norm
  resp_E['sensitivity'] = eval(resps[1]).gain * float(resps[0])

  resp_N = {}
  resp_N['poles'] = eval(resps[2]).p
  resp_N['zeros'] = eval(resps[2]).z[1:]
  resp_N['gain'] = eval(resps[2]).norm
  resp_N['sensitivity'] = eval(resps[2]).gain * float(resps[0])

  resp_Z = {}
  resp_Z['poles'] = eval(resps[3]).p
  resp_Z['zeros'] = eval(resps[3]).z[1:]
  resp_Z['gain'] = eval(resps[3]).norm
  resp_Z['sensitivity'] = eval(resps[3]).gain * float(resps[0])

  return resp_E,resp_N,resp_Z


def time_checker(time,changetimes):
  '''
  helper function to determine which set of responses to use for stations that have experienced instrument changes
  '''
  ln = len(changetimes)
  #case 1: before first changetime
  count = 0
  try:
    while time > changetimes[count]:
      count += 1
    return count
  except IndexError: #ran out of times
    return count



class Station:
  #Station class, filled from info_file format 

  def __init__(self,sta,info_path,network=None,stacorr=0,headerline=False):
    '''
    name - string; station name (has to be identical to name given in info_file)
    info_path - string; exact path of info_file to be used
    network - string; network name (use in case there's two stations of the same name in info_file)
    '''
    infodict_old = read_info(info_path,header=headerline)
    infodict = enter_4names(infodict_old)    

    if not sta in infodict.keys():
      if network == None:
        #print 'Station name not unique! Please provide network as input kwarg!'
	#return
        pass
      else:
        if sta+'_'+network in infodict.keys():      
          sta = sta+'_'+network
        else:
          print 'Station does not exist in provided info file!'
	  return

    try:
      self.lat = float(infodict[sta]['lat'])
    except ValueError:
      print 'Error: Provided Latitude value should be of type float!!'
      self.lat = 0.

    try:
      self.lon = float(infodict[sta]['lon'])
    except ValueError:
      print 'Error: Provided Longitude value should be of type float!!'
      self.lon = 0.

    try:
      self.elev = float(infodict[sta]['elev'])
    except ValueError:
      print 'Error: Provided Elevation value should be of type float!!'
      self.elev = 0.

    self.name = infodict[sta]['name']
    self.name4 = infodict[sta]['name4']
    self.network = infodict[sta]['network']
    self.archive = infodict[sta]['archive']
    self.pref = infodict[sta]['pref']
    self.stacorr = stacorr
    
    if not type(infodict[sta]['responses']) == dict: #simple (usual) case
      self.pazE,self.pazN,self.pazZ = find_resp(sta,infodict[sta]['responses'])
      self.instr = infodict[sta]['instruments']
      self.changetimes = []
    else: #different digitizer/seismometer configurations...additionally configure switch time (UTCDateTime)
      num = (len(infodict[sta]['responses']) + 1) / 2 
      self.pazE = []
      self.pazN = []
      self.pazZ = []
      self.instr = []
      self.changetimes = []
      for i in range(1,num+1):
        dum1,dum2,dum3 = find_resp(sta,infodict[sta]['responses']['config'+str(i)])
        xyz = infodict[sta]['instruments']['config'+str(i)]
        self.pazE.append(dum1)
        self.pazN.append(dum2)
        self.pazZ.append(dum3)
        self.instr.append(xyz)
        if i > 1:
          tm = infodict[sta]['responses']['switch'+str(i-1)+'_'+str(i)]
          dy,mon,yr = tm.split('/')
          utctm = UTCDateTime(int(yr),int(mon),int(dy),0,0,0)
          self.changetimes.append(utctm)

  def lat(self):
    return self.lat
  
  def lon(self):
    return self.lon
  
  def elev(self):
    return self.elev

  def name(self):
    return self.name

  #def rot(self):
  #return rotation angle (in case of misoriented sensors)

  def response(self,channel='Z',time='None'):
    '''
    return station response...in case of response differences between channels, enter channel as keyword argument; default is vertical channel
    beware of obspy's nomenclature: 'gain' refers to the A0 normalization factor, 'sensitivity' sums up the sensitivities of digitizer (incl. pre-gain) and sensor
    '''
    if len(self.changetimes) == 0:
      if channel == 'Z':
        return self.pazZ
      elif channel == 'N':
        return self.pazN
      elif channel == 'E':
        return self.pazE
      
    else:
      if time == 'None':
        print 'Error: Seismometer type cannot be uniquely determined (i.e. a device change occurred at some point in time); please provide a timestamp (UTCDateTime object) with the keyword argument time!'
      else:
        #determine which set to use
        indx = time_checker(time,self.changetimes)
        if channel == 'Z':
          return self.pazZ[indx]
        elif channel == 'N':
          return self.pazN[indx]
        elif channel == 'E':
          return self.pazE[indx]

  def archive(self):
    return self.archive

  def coords(self):
    return (self.lon,self.lat)

  def network(self):
    return self.network

  def uptime(self,channel='pref',comp='Z'):
    '''
    return dictionary of station running time
    channel: give channel identifier to be checked (e.g. HH, BH, SH)
    comp: add component to be checked (default: Z)
    '''
    from glob import glob
    from obspy.core import read
    if channel == 'pref':
      channel = self.pref
    yearlist = glob(self.archive+'/*')
    yrs = []
    for fls in yearlist:
      yr = int(fls.split('/')[-1])
      if not yr in yrs:
        yrs.append(yr)

    updict = {}
    for y in yrs:
      updict[y] = []
      llist = glob(self.archive+'/'+str(y)+'/'+self.network+'/'+self.name+'/'+channel+comp+'.D/*')
      for fl in llist:
        jday = fl.split('/')[-1].split('.')[-1]
        if not jday in updict[y]:
          try:
            updict[y].append(int(jday))
          except ValueError:
            continue

    for j in updict.keys():
      updict[j].sort()

    return updict

  def scan(self,ident='HH'):
    '''
    call to obspy scan
    add channel identifier (HH,BH,...) as kwarg ident
    '''
    print 'This calls obspy-scan for station '+self.name+'...this will take several minutes or up to an hour, depending on data volume. Results will be saved in file '+self.name+'_scan.pdf!'
    os.system('obspy-scan -o '+self.name+'_scan.pdf '+self.archive+'/*/'+self.network+'/'+self.name+'/'+ident+'*')
    #add '-o '+self.name+'_scan.pdf


  def seismometer(self,time='None'):
    if len(self.changetimes) == 0:
      return self.instr[0]
    else:
      if time == 'None':
        print 'Error: Seismometer type cannot be uniquely determined (i.e. a device change occurred at some point in time); please provide a timestamp (UTCDateTime object) with the keyword argument time!'
      else:
        #determine which set to use
        indx = time_checker(time,self.changetimes)
        return self.instr[indx][0]  

  def digitizer(self,time='None'):
    if len(self.changetimes) == 0:
      return self.instr[1]
    else:
      if time == 'None':
        print 'Error: Digitizer type cannot be uniquely determined (i.e. a device change occurred at some point in time); please provide a timestamp (UTCDateTime object) with the keyword argument time!'
      else:
        #determine which set to use
        indx = time_checker(time,self.changetimes)
        return self.instr[indx][1]


  def channels(self):
    ''' gives a list of all available channels'''

    from glob import glob
    #go through year folders, network and station set...list all channel folders, if anything new add to list
    channellist = []
    globlist = glob(self.archive+'/*/'+self.network+'/'+self.name+'/*')
    for entry in globlist:
      chn = entry.split('/')[-1].split('.')[0]
      if not chn in channellist:
        channellist.append(chn)

    return channellist

  def name4(self):
    return self.name4

  def stacorr(self):
   return self.stacorr


  def datafiles(self,startday,endday,channel='pref',comp='Z'):
    '''
    improved glob - provide list of file paths for station, channel and component between day startday and day endday (both UTCDateTime objects)
    choice of wildcard '?' for comp is allowed
    '''

    if channel == 'pref':
      channel = self.pref
    from glob import glob   

    jday1 = startday.julday
    year1 = startday.year

    jday2 = endday.julday
    year2 = endday.year

    outlist = []
    yearlist = glob(self.archive+'/*')
    for ntr in yearlist:
      yr = int(ntr.split('/')[-1])
      if year1 == year2: #special case...
        llist = glob(self.archive+'/'+str(yr)+'/'+self.network+'/'+self.name+'/'+channel+comp+'.D/*.[0-9][0-9][0-9]')
        for entry in llist:
          jday = int(entry.split('/')[-1].split('.')[-1])
          if yr == year1 and jday >= jday1 and jday <= jday2:
            outlist.append(entry)
      else:
        if yr >= year1 and yr <= year2:
          llist = glob(self.archive+'/'+str(yr)+'/'+self.network+'/'+self.name+'/'+channel+comp+'.D/*.[0-9][0-9][0-9]')
          for entry in llist:
            jday = int(entry.split('/')[-1].split('.')[-1])
            if yr == year1 and jday >= jday1:
              outlist.append(entry)
            elif yr != year1 and yr != year2:
              outlist.append(entry)
            elif yr == year2 and jday <= jday2:
              outlist.append(entry)

    outlist.sort()
    return outlist


class Event:
  """
  dic with keys: lat,lon,dep,nobs,origin_time,magnitude,gap,RMS_residual,cumsum,cumsumS,Picks,moment_tensor,DC%,Mw,centroid_depth,
  """
  #MUST have: location (lat,lon,dep,orig_time)
  #everything initialized with dummies

  def __init__(self,origin_time,lat,lon,dep,info={},Picks={}):
    
    info['origin_time'] = origin_time #check whether it's a valid UTCDateTime object!!
    info['event_lon'] = float(lon)
    info['event_lat'] = float(lat)
    info['event_dep'] = float(dep)
    info.setdefault('cumsum',-99)
    info.setdefault('cumsumS',-99)
    info.setdefault('nobs',-99)
    info.setdefault('magnitude',0.)
    info.setdefault('gap',199)
    info.setdefault('RMS_residual',-1.)
    info.setdefault('MomentTensor','None')
    info.setdefault('strike/dip/rake','None')
    info.setdefault('DC%',-99.)
    info.setdefault('Mw',-99.)
    info.setdefault('centroid_depth',-99.)

    self.info = info
    self.Picks = Picks

  def picks(self):
    return self.Picks



class Catalog:
  """
  time-sorted dictionary of Events
  """

  def __init__(self,events={},info={}):
    """
    format - select file format for reading of event list/catalog
             calls routines from util
    """
    info.setdefault('source','not specified')    
    info.setdefault('M_c','not specified') #completeness magnitude
    info.setdefault('Nev',0) #number of events

    self.events = events
    self.info = info

  def nev(self):
    return self.info['Nev']

  def filter(self,**kwargs):
    """
    allowed kwargs:
    gap_max, nobs_min, cumsum_min, cumsumS_min - int values
    mag_min, mag_max, rms_max, lat_min, lat_max, lon_min, lon_max, dep_min, dep_max  - float values
    starttime, endtime - UTCDateTime objects 
    ident - str; possible: 'NN', 'ID', 'P1', 'P2', 'P3', 'IF', 'UP'; works only if labels are set
    filter catalog by different criteria...modifies Catalog object in-place by deleting Events that do not fit criteria 
    also throws out events where the parameter in question is not defined/set to dummy

    """
    #kwargs with both-way selection -> numeric flag
    latflag=0 #1 -> only lower end, 2 -> only upper end, 3 -> both ends
    depflag=0
    lonflag=0
    magflag=0
    timeflag=0

    #kwargs with one-way selection -> boolean flag
    gapflag=False
    nobsflag=False
    cumflag=False
    cumSflag=False
    rmsflag=False
    identflag=False

    for key,value in kwargs.iteritems(): #other kwargs than the ones listed here are ignored
      if key == 'starttime':
        timeflag += 1
      elif key == 'endtime':
        timeflag += 2
      elif key == 'lat_max':
        latflag += 2
      elif key == 'lat_min':
        latflag += 1       
      elif key == 'lon_max':       
        lonflag += 2
      elif key == 'lon_min':       
        lonflag += 1
      elif key == 'dep_max':       
        depflag += 2
      elif key == 'dep_min':       
        depflag += 1
      elif key == 'mag_max':       
        magflag += 2
      elif key == 'mag_min':       
        magflag += 1
      elif key == 'nobs_min':       
        nobsflag = True
      elif key == 'gap_max':       
        gapflag = True
      elif key == 'rms_max':       
        rmsflag = True
      elif key == 'cumsum_min':       
        cumflag = True
      elif key == 'cumsumS_min':
        cumSflag = True
      elif key == 'ident':
        identflag = True
    
    if identflag:
      for j in self.events.keys():
        if self.events[j].info['ident'] != kwargs['ident']:
          del self.events[j]

    for j in self.events.keys(): #now go through catalog
      #booleans first
      if cumSflag:
        if self.events[j].info['cumsumS'] == -99 or self.events[j].info['cumsumS'] < kwargs['cumsumS_min']:
          del self.events[j]
          continue  
      if cumflag:
        if self.events[j].info['cumsum'] == -99 or self.events[j].info['cumsum'] < kwargs['cumsum_min']:
          del self.events[j]
          continue
      if rmsflag:
        if self.events[j].info['RMS_residual'] == -99. or self.events[j].info['RMS_residual'] > kwargs['rms_max']:
          del self.events[j]
          continue
      if gapflag:
        if self.events[j].info['gap'] == -99 or self.events[j].info['gap'] > kwargs['gap_max']:
          del self.events[j]
          continue
      if nobsflag:
        if self.events[j].info['nobs'] == -99 or self.events[j].info['nobs'] < kwargs['nobs_min']:
          del self.events[j]
          continue

      #now the numeric ones...
      if timeflag == 1: #only minimum set
        if (self.events[j].info['origin_time'] - kwargs['starttime']) < 0.:
          del self.events[j]
          continue
      elif timeflag == 2: #only maximum set
        if (self.events[j].info['origin_time'] - kwargs['endtime']) > 0.:
          del self.events[j]
          continue
      elif timeflag == 3: #both set
        if (self.events[j].info['origin_time'] - kwargs['endtime']) > 0. or (self.events[j].info['origin_time'] - kwargs['starttime']) < 0.:
          del self.events[j]
          continue

      if latflag == 1: #only minimum set
        if self.events[j].info['event_lat'] < kwargs['lat_min']:
          del self.events[j]
          continue
      elif latflag == 2: #only maximum set
        if self.events[j].info['event_lat'] > kwargs['lat_max']:
          del self.events[j]
          continue
      elif latflag == 3: #both set
        if self.events[j].info['event_lat'] > kwargs['lat_max'] or self.events[j].info['event_lat'] < kwargs['lat_min']:
          del self.events[j]
          continue

      if lonflag == 1: #only minimum set
        if self.events[j].info['event_lon'] < kwargs['lon_min']:
          del self.events[j]
          continue
      elif lonflag == 2: #only maximum set
        if self.events[j].info['event_lon'] > kwargs['lon_max']:
          del self.events[j]
          continue
      elif lonflag == 3: #both set
        if self.events[j].info['event_lon'] > kwargs['lon_max'] or self.events[j].info['event_lon'] < kwargs['lon_min']:
          del self.events[j]
          continue

      if depflag == 1: #only minimum set
        if self.events[j].info['event_dep'] < kwargs['dep_min']:
          del self.events[j]
          continue
      elif depflag == 2: #only maximum set
        if self.events[j].info['event_dep'] > kwargs['dep_max']:
          del self.events[j]
          continue
      elif depflag == 3: #both set
        if self.events[j].info['event_dep'] > kwargs['dep_max'] or self.events[j].info['event_dep'] < kwargs['dep_min']:
          del self.events[j]
          continue

      if magflag == 1: #only minimum set
        if self.events[j].info['magnitude'] < kwargs['mag_min']:
          del self.events[j]
          continue
      elif magflag == 2: #only maximum set
        if self.events[j].info['magnitude'] > kwargs['mag_max']:
          del self.events[j]
          continue
      elif magflag == 3: #both set
        if self.events[j].info['magnitude'] > kwargs['mag_max'] or self.events[j].info['magnitude'] < kwargs['mag_min']:
          del self.events[j]
          continue

      
  def filter_slabdist(self,model='slab1.0',dist_min=0.,dist_max=999.):
    """
    filter catalog by hypocenter distance from slab surface
    """
    Lats,Lons,Deps = util.import_slab(model=model)
    for j in self.events.keys():
      prof = util.get_closest_profile(self.events[j].info['event_lat'],Lats,Lons,Deps)
      dist = util.compute_point_distance(self.events[j].info['event_lon'],self.events[j].info['event_dep'],prof,Lons)
      if dist < dist_min:
        del self.events[j]
        continue
      if dist > dist_max:
        del self.events[j]
        continue

  #plotting tools as functions here!!!
  
  #comparison between two catalogs (e.g. extract common events and events that are unique to either)
  def compare(self,other,output='unique',window=5.):
    """
    other - a second catalog object
    output - 'unique' gives back all events that are in catalog but not in other
             'other' gives back all events that are in other but not in this catalog
             'both1', 'both2' gives back all events that are in both catalogs (with metadata from catalog 1 or 2, respectively)
    window: search window around origin times in seconds

    """
    blup = 0
    #search for events with similar origin time (window!); somehow make use of fact that Catalog is always time-sorted
    #once candidate event pair exists, compare location/depth
    #three output dicts (return the one the user wants -> output variable)

    cat_both1 = Catalog() #events that appear in both catalogs, metadata from catalog 1 is kept
    cat_both1.events = {}

    cat_both2 = Catalog() #events that appear in both catalogs, metadata from catalog 2 is kept
    cat_both2.events = {}

    cat_other = Catalog() #events that only show up in catalog 2
    cat_other.events = {}

    cat_unique = Catalog() #events that only show up in catalog 1
    cat_unique.events = {}


    lendic1 = len(self.events.keys())
    lendic2 = len(other.events.keys())

    key1 = 0
    key2 = 0 
    bkp = 0
    flag = False

    orig_end = other.events[sorted(other.events.keys())[lendic2-1]].info['origin_time']

    while key1 < lendic1:
      if not flag:
        key2 = bkp
      orig = self.events[sorted(self.events.keys())[key1]].info['origin_time']
      while key2 < lendic2:
        orig2 = other.events[sorted(other.events.keys())[key2]].info['origin_time']
        if orig2 > orig + window:
          key1 += 1
          flag = True
          key2 -= 1
          break #already past
        elif abs(orig - orig2) <= window: #found one
          #check if location is at least roughly similar
          cat_both1.events[sorted(self.events.keys())[key1]] = cp(self.events[sorted(self.events.keys())[key1]])
          cat_both2.events[sorted(other.events.keys())[key2]] = cp(other.events[sorted(other.events.keys())[key2]])
          key1 += 1
          bkp = cp(key2)
          flag = True
          break
        elif orig > orig_end: #gone beyond other catalog
          key1 = lendic1
          break
        key2 += 1
      flag = False

    #if output was set to both1 or both2, return these catalogs
    if output == 'both1':
      return cat_both1
    elif output == 'both2':
      return cat_both2
    elif output == 'other': #get events that are unique to catalog 2
      for j in other.events.keys():
        if not j in cat_both2.events.keys():
          cat_other.events[j] = cp(other.events[j])
      return cat_other
    elif output == 'unique':
      for k in self.events.keys():
        if not k in cat_both1.events.keys():
          cat_unique.events[k] = cp(self.events[k])
      return cat_unique

  def plot_map(self,llcrnrlat='none',llcrnrlon='none',urcrnrlat='none',urcrnrlon='none',save=False,colorcode=True,dots=True,mines=False,lock=False,trench=True):
    """
    catalog map plotting
    """
    
    lats = []
    lons = []
    deps = []
    mags = []
    for j in sorted(self.events.keys()):
      lats.append(self.events[j].info['event_lat'])
      lons.append(self.events[j].info['event_lon'])
      deps.append(self.events[j].info['event_dep'])
      mags.append(self.events[j].info['magnitude'])

    mxlat = (array(lats)).max()
    mnlat = (array(lats)).min()
    mxlon = (array(lons)).max()
    mnlon = (array(lons)).min()

    figure(figsize=(10,15))

    from mpl_toolkits import basemap
    if llcrnrlat == 'none':
      llcrnrlat = mnlat-0.2
    if llcrnrlon == 'none':
      llcrnrlon = mnlon - 0.2
    if urcrnrlat == 'none':
      urcrnrlat = mxlat + 0.2
    if urcrnrlon == 'none':
      urcrnrlon = mxlon + 0.2

    m = basemap.Basemap(projection='merc',llcrnrlat=llcrnrlat,llcrnrlon=llcrnrlon,urcrnrlat=urcrnrlat,urcrnrlon=urcrnrlon,resolution='i')
    m.drawcoastlines()
    m.drawcountries() 

    mers = range(int(mnlon)-1,int(mxlon)+1,2)
    pars = range(int(mnlat)-1,int(mxlat)+1,2)

    m.drawmeridians(mers,labels=[0,0,0,1])

    m.drawparallels(pars,labels=[1,0,0,0])
    if lock:
      Lat,Lon,Lock = util.import_locking('/home/sippl/temp/Locking_NChile.xyz')
      Xs,Ys = meshgrid(Lon,Lat)
      X,Y = m(Xs,Ys)
      m.pcolor(X,Y,Lock,cmap='hot_r')
      colorbar()
      m.contour(X,Y,Lock,levels=[0.05,0.25,0.5,0.7,0.85],colors='blue')

    if trench:
      infl = open('/home/sippl/temp/trench_Chile.txt','r')
      dat = infl.readlines()
      infl.close()
    
      latl = []
      lonl = []
      for j in dat:
        lon,lat = j.split(None)
        lonl.append(float(lon))
        latl.append(float(lat))
      XX,YY = m(lonl,latl)
      m.plot(XX,YY,'k-')
    '''
    #temp
    stalats = [-21.04323,-21.31973,-22.04847,-22.33369,-22.85283,-22.70580,-21.72667,-20.14112,-21.79638,-23.51343,-19.76096,-18.61406,-18.33585,-24.62597,-23.20833,-18.33510,-21.53414,-17.58594,-20.27822,-19.13108,-19.59717,-20.82071,-23.51149,-21.3948,-20.93345,-20.65325,-20.60909,-20.24393,-20.7613,-22.6182,-19.6685,-18.3708,-20.5656]
    stalons = [-69.48740,-69.89603,-69.75310,-70.14918,-70.20235,-69.57188,-69.88618,-69.15340,-69.24192,-70.55408,-69.65582,-70.32808,-69.50160,-70.40379,-69.47092,-69.50767,-68.73017,-68.48000,-69.88791,-69.59553,-70.12305,-70.15288,-70.24961,-69.6101,-68.87026,-69.11923,-68.79609,-70.14041,-69.3973,-68.9113,-69.1942,-70.3419,-70.1807]
    for j in range(len(stalats)):
      CorX,CorY = m(stalons[j],stalats[j])
      plot(CorX,CorY,'^',markersize=11,markerfacecolor='y',markeredgecolor='k')
    '''
    x,y = m(lons,lats)
    if colorcode:
      if dots:
        m.scatter(x,y,c=deps,s=2,alpha=0.7,linewidths=0)
      else:
        m.scatter(x,y,s=5,c=deps,facecolors='none')
      colorbar(shrink=0.5,label='depth [km]')
    else:
      if dots:
        m.scatter(x,y,c='black',s=2,alpha=0.7,linewidths=0)
      else:
        m.scatter(x,y,s=array(mags)*4,facecolors='none',edgecolors='black')

    if mines:
      minelats = [-22.310,-22.332,-21.004,-20.9665,-20.9766,-20.886,-22.3844,-20.0504,-20.4077,-23.4303,-22.7908,-21.986,-21.919,-24.269,-20.1952,-22.6819]
      minelons = [-68.908,-69.653,-68.808,-68.721,-69.6824,-70.032,-70.2118,-69.257,-69.868,-70.0568,-69.2552,-68.7008,-68.8307,-69.066,-69.3362,-70.1766]

      nilslats = [-22.2097,-22.3806,-22.5210,-22.8198,-22.6359,-22.5680,-22.6785,-23.4333,-23.4312,-22.9742,-22.7916,-23.1227,-22.8596,-22.1297,-20.8879,-20.9497,-20.9896,-20.9806,-20.8943,-20.4869,-20.4322,-20.2546,-20.2171,-20.0604,-22.6072]
      nilslons = [-69.6510,-70.2116,-69.8789,-69.7205,-69.8829,-69.8125,-70.1691,-70.0635,-69.5130,-69.0622,-69.2572,-69.8556,-69.3253,-69.8607,-70.0320,-69.9924,-69.6557,-69.5920,-69.6541,-69.9489,-69.8669,-69.8339,-69.7951,-69.2697,-69.6721]

      for j in range(len(minelons)):
        x_mine,y_mine = m(minelons[j],minelats[j])
        m.plot(x_mine,y_mine,'rv',markersize=12)

      for k in range(len(nilslats)):
        x_nils,y_nils = m(nilslons[k],nilslats[k])
        m.plot(x_nils,y_nils,'gv',markersize=10)
    
    if save:
      savefig('output.pdf')


  def plot_ewsection(self,centr_lat='none',width='none',maxlon='none',minlon='none',save=False,mag=True,hist=False,slab=False,slabmod='tassara',ident=False):
    """
    plot catalog onto E-W section
    """

        #everything is 'none' --> get coordinate max/min from locations
    lats = []
    lons = []
    deps = []
    mags = []
    idents = [] 

    if centr_lat == 'none':
      print 'Error: You have to specify a central latitude for the profile'
      return

    cat_new = cp(self)
    if not maxlon == 'none' and not minlon == 'none':
      cat_new.filter(lon_max=maxlon,lon_min=minlon)

    if not width == 'none':
      mxlat = centr_lat + (width/111.195)
      mnlat = centr_lat - (width/111.195)

      cat_new.filter(lat_max=mxlat,lat_min=mnlat)

    for j in sorted(cat_new.events.keys()):
      lats.append(cat_new.events[j].info['event_lat'])
      lons.append(cat_new.events[j].info['event_lon'])
      deps.append(cat_new.events[j].info['event_dep'])
      mags.append(cat_new.events[j].info['magnitude'])
      if ident:
        idents.append(cat_new.events[j].info['ident'])

    mxlon = (array(lons)).max()
    mnlon = (array(lons)).min()
    mxdep = (array(deps)).max()

    mxlon= -65.
    mnlon= -73.
    mxdep=300.

    cat_new.filter(lon_max=mxlon,lon_min=mnlon,dep_max=mxdep)

    lat_mean = (array(lats)).mean()
    fac = cos(lat_mean * pi / 180.) * 111.195

    lonrange = ((mxlon - mnlon)) * fac

    ysize = 16 * (mxdep/lonrange) 

    figure(figsize=(16,ysize))

    if slab:
      Lats, Lons, Deps = util.import_slab(model='tassara')
      prof = util.get_closest_profile(centr_lat,Lats,Lons,Deps)
      plot((array(Lons)-mnlon)*fac,prof,'r-')

      Lats,Lons,Deps = util.import_slab(model='slab1.0')
      prof = util.get_closest_profile(centr_lat,Lats,Lons,Deps)
      plot((array(Lons)-mnlon)*fac,prof,'g-')

      Lats,Lons,Deps = util.import_slab(model='ipoc')
      prof = util.get_closest_profile(centr_lat,Lats,Lons,Deps)
      plot((array(Lons)-mnlon)*fac,prof,'b-')


    import matplotlib as mpl
    if hist:
      hist2d(((array(lons)-mnlon)*fac),(array(deps)*(-1)),(225,125),norm=mpl.colors.LogNorm(),cmap=cm.jet)
    
    elif mag:
      if not ident:
        scatter(((array(lons)-mnlon)*fac),(array(deps)*(-1)),s=array(mags)*3.5,facecolor='none',edgecolors='black')
      else:
        for j in range(len(lons)):
          if idents[j] == 'ID':
            plot((lons[j]-mnlon)*fac,deps[j]*(-1),'o',markersize=4,markerfacecolor='none',markeredgecolor='cyan')
          elif idents[j] == 'NN':
            plot((lons[j]-mnlon)*fac,deps[j]*(-1),'o',markersize=4,markerfacecolor='none',markeredgecolor='black')
          elif idents[j] == 'P1':
            plot((lons[j]-mnlon)*fac,deps[j]*(-1),'o',markersize=4,markerfacecolor='none',markeredgecolor='blue')
          elif idents[j] == 'P2':
            plot((lons[j]-mnlon)*fac,deps[j]*(-1),'o',markersize=4,markerfacecolor='none',markeredgecolor='green')
          elif idents[j] == 'P3':
            plot((lons[j]-mnlon)*fac,deps[j]*(-1),'o',markersize=4,markerfacecolor='none',markeredgecolor='red')
          elif idents[j] == 'UP':
            plot((lons[j]-mnlon)*fac,deps[j]*(-1),'o',markersize=4,markerfacecolor='none',markeredgecolor='magenta')
    
    else:
      if not ident:
        scatter(((array(lons)-mnlon)*fac),(array(deps)*(-1)),s=0.5,facecolors='none',edgecolors='black')
      else:
        for j in range(len(lons)):
          if idents[j] == 'ID':
            plot((lons[j]-mnlon)*fac,deps[j]*(-1),'o',markersize=4,markerfacecolor='none',markeredgecolor='cyan')
          elif idents[j] == 'NN':
            plot((lons[j]-mnlon)*fac,deps[j]*(-1),'o',markersize=4,markerfacecolor='none',markeredgecolor='black')
          elif idents[j] == 'P1':
            plot((lons[j]-mnlon)*fac,deps[j]*(-1),'o',markersize=4,markerfacecolor='none',markeredgecolor='blue')
          elif idents[j] == 'P2':
            plot((lons[j]-mnlon)*fac,deps[j]*(-1),'o',markersize=4,markerfacecolor='none',markeredgecolor='green')
          elif idents[j] == 'P3':
            plot((lons[j]-mnlon)*fac,deps[j]*(-1),'o',markersize=4,markerfacecolor='none',markeredgecolor='red')
          elif idents[j] == 'UP':
            plot((lons[j]-mnlon)*fac,deps[j]*(-1),'o',markersize=4,markerfacecolor='none',markeredgecolor='magenta')
    


    #for j in range(len(lons)):
    #  plot((lons[j]-mnlon)*fac,deps[j]*(-1),'o',markersize=4,markerfacecolor='none',markeredgecolor='black')
    axis('equal')
    xlabel('distance [km]')
    ylabel('depth [km]')

    #plot([0,350],[-100,-100],'r--')
    #plot([0,350],[-110,-110],'r--')
    #plot([0,350],[-120,-120],'r--')
    

    ylim([-300,0])
    xlim([0,830])
    #axis('off')

    if save:
      savefig(save+'.pdf')

  def plot_topoprof(centr_lat):
    """

    """
    Lats, Lons, Deps = util.import_slab(model='etopo')
    prof = util.get_closest_profile(centr_lat,Lats,Lons,Deps)
    fac = cos(centr_lat * pi / 180.) * 111.195   

    plot(array(Lons)*fac,prof,'k--') 



  def plot_anysection(self,startpoint,endpoint,width='none',save=False):
    """
    plot arbitrarily oriented cross section, input: two points (start, end)
    startpoint, endpoint: [lat,lon]

    """
    lats = []
    lons = []
    deps = []
    mags = []
    for j in sorted(self.events.keys()):
      lats.append(self.events[j].info['event_lat'])
      lons.append(self.events[j].info['event_lon'])
      deps.append(self.events[j].info['event_dep'])
      mags.append(self.events[j].info['magnitude'])

    mxdep = (array(deps)).max()

    #now configure profile line
    ybig = max(startpoint[0],endpoint[0])
    ysmall = min(startpoint[0],endpoint[0])
    xbig = max(startpoint[1],endpoint[1])
    xsmall = min(startpoint[1],endpoint[1])
    dx = (xbig - xsmall) * 111.195 * cos((ybig+ysmall)/2. * pi / 180.)
    dy = (ybig - ysmall) * 111.195
    try:
      az = (arctan(dy/(float(dx)))*180.)/pi
    except ZeroDivisionError:
      az = 0.
    length = sqrt(dx**2 + dy**2)

    ysize = 12 * (mxdep/length)
    figure(figsize=(12,ysize))

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
    angle_norm = az%90

    for j in range(len(lats)):
      ev_offset_n = ((startpoint[0] - lats[j]) * 111.195) / sin(angle_norm*pi/180.)
      ev_offset_e = (-1) * ((startpoint[1] - lons[j]) * 111.195 * cos(startpoint[0] * pi /180.)) / sin((90.-angle_norm)*pi/180.)

      #check whether event is within profile line
      xstart = 0
      ystart = 0
      xend = (endpoint[1] - startpoint[1]) * 111.195 * cos((endpoint[1] + startpoint[1])/2. * pi / 180.) 
      yend = (endpoint[0] - startpoint[0]) * 111.195
      xev = (lons[j] - startpoint[1]) * 111.195 * cos((lons[j] + startpoint[1])/2. * pi / 180.)
      yev = (lats[j] - startpoint[0]) * 111.195

      dist = ((yend - ystart) * xev - (xend - xstart) * yev + xend*ystart - yend*xstart) / sqrt((yend - ystart)**2 + (xend - xstart)**2)

      if abs(dist) <= width:

        if abs(az-90.) > 30. and abs(az-270.) > 30.:
          n_correction = dist / tan(angle_norm*pi/180.)
          ev_offset = ev_offset_n + n_correction

        else:
          e_correction = dist / tan((90 - angle_norm)*pi/180.)
          ev_offset = ev_offset_e - e_correction

        if ev_offset > 0 and ev_offset < length:
          #add station to CCP stack
          plot(ev_offset,deps[j]*(-1),marker='.',markerfacecolor='None',markeredgecolor='green')

    axis('equal')
    xlabel('distance [km]')
    ylabel('depth [km]')
    ylim([(-1)*mxdep,0])
    xlim([0,length])
    if save:
      savefig('output.pdf')



  def plot_nssection(self,centr_lon='none',width='none',maxlat='none',minlat='none',save=False,mag=False):
    """
    plot catalog onto N-S section
    """

            #everything is 'none' --> get coordinate max/min from locations
    lats = []
    lons = []
    deps = []
    mags = []

    if centr_lon == 'none':
      print 'Error: You have to specify a central longitude for the profile'
      return
    for i in sorted(self.events.keys()):
      lats.append(self.events[i].info['event_lat'])

    mean_lat = (array(lats)).mean()
    fac = 111.195 * cos(mean_lat*pi/180.)

    cat_new = cp(self)
    if not maxlat == 'none' and not minlat == 'none':
      cat_new.filter(lat_max=maxlat,lat_min=minlat)

    if not width == 'none':
      mxlon = centr_lon + (width/fac)
      mnlon = centr_lon - (width/fac)

      cat_new.filter(lon_max=mxlon,lon_min=mnlon)

    lats = []
    for j in sorted(cat_new.events.keys()):
      lats.append(cat_new.events[j].info['event_lat'])
      lons.append(cat_new.events[j].info['event_lon'])
      deps.append(cat_new.events[j].info['event_dep'])
      mags.append(cat_new.events[j].info['magnitude'])

    mxdep = (array(deps)).max()
    mxlat = (array(lats)).max()
    mnlat = (array(lats)).min()

    latrange = (mxlat - mnlat) * 111.195    
    ysize = 20 * (mxdep/latrange)

    figure(figsize=(16,6))
    #figure()

    if mag:
      scatter(((array(lats)-mnlat)*111.195),(array(deps)*(-1)),s=array(mags)*4,facecolor='none',edgecolors='black')
    else:
      scatter(((array(lats)-mnlat)*111.195),(array(deps)*(-1)),s=0.5,facecolors='none',edgecolors='black')

    axis('equal')

    ylabel('depth [km]')
    xlabel('distance [km]')

    Send = (-24.75 - mnlat)*111.195
    Nend = (-18. - mnlat)*111.195

    print Send, Nend

    xlim([Send,Nend])
    ylim([-280,0])

    if save:
      savefig('section.pdf')


  def Npicks(self):
    """
    return total number of picks (P/S/combined) contained in the catalog
    """
    np = 0
    ns = 0
    for u in self.events.keys():
      for stat in self.events[u].Picks.keys():
        if self.events[u].Picks[stat].has_key('P-pick'):
          np += 1
        if self.events[u].Picks[stat].has_key('S-pick'):
          ns += 1

    ntot = np + ns

    return np,ns,ntot

  def plot_histog_mag(self,bins=[],xlimits=[],xticks=[],**kwargs):
    """
    plot histogram of magnitudes of events contained in catalog.
    bins = [min_bin, max_bin, step_bin]: list to set bins shape.
    xlim = [xmin, xmax]: list of x axis limits to plot.
    xticks = [min, max, step]: list to set ticks shape.
    allowed kwargs (most of them based on plt.hist documentation):
      hist_log, hist_cumul, hist_normed, hist_combined (this adds a logarithmic histogram to original histogram) - boolean values
      hist_cumul - integer values (it can also be int, see plt.hist documentation)
      hist_type - string
    """

    import numpy as np
    import subprocess

    outpath = "test_magnitud"
    cmd1 = "mkdir -p " + outpath
    subprocess.call(cmd1,shell=True,stdout=subprocess.PIPE)

    #kwargs by default
    hist_logflag = False
    hist_cumulflag = False
    hist_typeflag = 'bar'
    hist_combinedflag = False
    hist_normedflag = False

    for key,value in kwargs.iteritems(): #other kwargs than the ones listed here are ignored
      if key == 'hist_log':
        hist_logflag = kwargs['hist_log']
      elif key == 'hist_cumul':
        hist_cumulflag = kwargs['hist_cumul']
      elif key == 'hist_type':
        hist_typeflag = kwargs['hist_type']
      elif key == 'hist_combined':
        hist_combinedflag = kwargs['hist_combined']
      elif key == 'hist_normed':
        hist_normedflag = kwargs['hist_normed']

    mag_arr = np.empty([len(self.events.keys())])
    n = 0
    for k in self.events.keys():
      ev = self.events[k]
      mag_arr[n] = ev.info['magnitude']
      n += 1 

    print str(len(self.events.keys())) + ' events found in catalog'
    print 'Max magnitude = ' + str(mag_arr.max())
    print 'Min magnitude = ' + str(mag_arr.min())

    max_weight = 0.1
    fac=10
    max_bin = int(max_weight*(int(np.amax(mag_arr)/max_weight)+1))
    if len(bins) != 3:
      bins = [0, max_bin, max_weight]
    nbin = np.arange(bins[0], bins[1], bins[2])       

    if len(xlimits) != 2:
      xlimits = [cumsum_arr.min(), max_bin]

    if len(xticks) != 3:
      xticks = [xlimits[0], xlimits[1], max_weight*fac]

    ticks = np.arange(xticks[0], xticks[1], xticks[2])
    hrange=(xlimits[0], xlimits[1])

    hy, hx, patches = plt.hist(mag_arr,nbin,range=hrange,normed=hist_normedflag,log=hist_logflag,cumulative=hist_cumulflag,histtype=hist_typeflag,facecolor='green',alpha=0.75)
    if hist_normedflag:
      ymin = 0
      ymax = hy.max()
      ylabel_str='Density'
      hist_str='normalized'      
    if hist_logflag:
      ymin = 1
      ymax = hy.max() + 1000
      ylabel_str = 'Log(Nev)'
      hist_str = 'counts_log'
    if hist_logflag and (hist_cumulflag or hist_cumulflag==-1):
      ymin = 1
      ymax = hy.max() + 1000
      ylabel_str = 'Log(Nev)'
      hist_str = 'counts_log_cumul'
    if not hist_normedflag and not hist_logflag and not hist_cumulflag and hist_cumulflag!=-1: 
      ymin = 0
      ymax = hy.max() + 100
      ylabel_str = 'Number of events'
      hist_str = 'counts'
    if hist_combinedflag: 
      plt.hist(mag_arr,nbin,range=hrange,log=True,cumulative=False,histtype='bar',facecolor='green',alpha=0.75)
      ymin = 1
      ymax = hy.max() + 1000
      ylabel_str='Log(Nev)'
      hist_str='combined'    

    plt.axis([xlimits[0], xlimits[1], ymin, ymax])
    plt.xticks(ticks)
    title_str = 'Histogram of magnitudes'
    plt.title(title_str)
    plt.xlabel('Magnitude')
    plt.ylabel(ylabel_str)
    plt.grid(True)
    f_out= outpath + '/hist_mag_' + hist_str + '.pdf'
    plt.savefig(f_out)
    plt.close()


  def plot_histog_cumsum(self,bins=[],xlimits=[],xticks=[],**kwargs):
    """
    plot histogram of cumulative sum of class points of P,S picks contained in catalog.
    bins = [min_bin, max_bin, step_bin]: list to set bins shape.
    xlim = [xmin, xmax]: list of x axis limits to plot.
    allowed kwargs (most of them based on plt.hist documentation):
      hist_log, hist_cumul, hist_normed, hist_combined (this adds a logarithmic histogram to original histogram) - boolean values
      hist_cumul - integer values (it can also be int, see plt.hist documentation)
      hist_type, hist_pick('P'/'S') - string
    xticks = [min, max, step]: list to set ticks shape.
    """
    import numpy as np
    import subprocess
    from util import get_cumsum,get_cumsumP,get_cumsumS

    outpath = "test_pick_class"
    cmd1 = "mkdir -p " + outpath
    subprocess.call(cmd1,shell=True,stdout=subprocess.PIPE)

    #kwargs by default
    hist_pickflag = 'P'
    hist_logflag = False
    hist_cumulflag = False
    hist_typeflag = 'bar'
    hist_combinedflag = False
    hist_normedflag = False

    for key,value in kwargs.iteritems(): #other kwargs than the ones listed here are ignored
      if key == 'hist_pick':
        hist_pickflag = kwargs['hist_pick']
      elif key == 'hist_log':
        hist_logflag = kwargs['hist_log']
      elif key == 'hist_cumul':
        hist_cumulflag = kwargs['hist_cumul']
      elif key == 'hist_type':
        hist_typeflag = kwargs['hist_type']
      elif key == 'hist_combined':
        hist_combinedflag = kwargs['hist_combined']
      elif key == 'hist_normed':
        hist_normedflag = kwargs['hist_normed']

    cumsum_arr = np.empty([len(self.events.keys())])
    n = 0
    for k in self.events.keys():
      ev = self.events[k]
      if hist_pickflag == 'P':
        cumsum_arr[n] = get_cumsumP(ev)
      elif hist_pickflag == 'S':
        cumsum_arr[n] = get_cumsumS(ev)     
      n += 1 

    print str(len(self.events.keys())) + ' events found in catalog'
    print 'Cumsum_max = ' + str(np.amax(cumsum_arr))
    print 'Cumsum_min = ' + str(np.amin(cumsum_arr))

    if hist_pickflag == 'P':
      max_weight = 4.
      fac=5
    elif hist_pickflag == 'S':
      max_weight = 8.
      fac=5
    max_bin = int(max_weight*(int(np.amax(cumsum_arr)/max_weight)+1))

    if len(bins) != 3:
      bins = [0, max_bin, max_weight]
    nbin = np.arange(bins[0], bins[1], bins[2])       

    if len(xlimits) != 2:
      xlimits = [cumsum_arr.min(), max_bin]

    if len(xticks) != 3:
      xticks = [xlimits[0], xlimits[1], max_weight*fac]

    ticks = np.arange(xticks[0], xticks[1], xticks[2])
    hrange=(xlimits[0], xlimits[1])

    hy, hx, patches = plt.hist(cumsum_arr,nbin,range=hrange,normed=hist_normedflag,log=hist_logflag,cumulative=hist_cumulflag,histtype=hist_typeflag,facecolor='green',alpha=0.75)
    if hist_normedflag:
      ymin = 0
      ymax = hy.max()
      ylabel_str='Density'
      hist_str='normalized'      
    if hist_logflag:
      ymin = 1
      ymax = hy.max() + 1000
      ylabel_str = 'Log(Nev)'
      hist_str = 'counts_log'
    if hist_logflag and (hist_cumulflag or hist_cumulflag==-1):
      ymin = 1
      ymax = hy.max() + 1000
      ylabel_str = 'Log(Nev)'
      hist_str = 'counts_log_cumul'
    if not hist_normedflag and not hist_logflag and not hist_cumulflag and hist_cumulflag!=-1: 
      ymin = 0
      ymax = hy.max() + 100
      ylabel_str = 'Number of events'
      hist_str = 'counts'
    if hist_combinedflag: 
      plt.hist(cumsum_arr,nbin,range=hrange,log=True,cumulative=False,histtype='bar',facecolor='green',alpha=0.75)
      ymin = 1
      ymax = hy.max() + 1000
      ylabel_str='Log(Nev)'
      hist_str='combined'    

    plt.axis([xlimits[0], xlimits[1], ymin, ymax])
    plt.xticks(ticks)
    title_str = 'Histogram of class points per event: ' + hist_pickflag + '-picks'
    plt.title(title_str)
    plt.xlabel('Class points')
    plt.ylabel(ylabel_str)
    plt.grid(True)
    f_out= outpath + '/hist_cumsum_' + hist_pickflag + '-picks_' + hist_str + '.pdf'
    plt.savefig(f_out)
    plt.close()


  def plot_histog_npicks(self,bins=[],xlimits=[],xticks=[],**kwargs):
    """
    plot histogram of the number of P,S picks per event in catalog.
    bins = [min_bin, max_bin, step_bin]: list to set bins shape.
    xlim = [xmin, xmax]: list of x axis limits to plot.
    xticks = [min, max, step]: list to set ticks shape.
    allowed kwargs (most of them based on plt.hist documentation):
      hist_log, hist_cumul, hist_normed, hist_combined (this adds a logarithmic histogram to original histogram) - boolean values
      hist_cumul - integer values (it can also be int, see plt.hist documentation)
      hist_type, hist_pick('P'/'S') - string
    """
    import numpy as np
    import subprocess

    outpath = "test_pick_number"
    cmd1 = "mkdir -p " + outpath
    subprocess.call(cmd1,shell=True,stdout=subprocess.PIPE)

    #kwargs by default
    hist_pickflag = 'P'
    hist_logflag = False
    hist_cumulflag = False
    hist_typeflag = 'bar'
    hist_combinedflag = False
    hist_normedflag = False

    for key,value in kwargs.iteritems(): #other kwargs than the ones listed here are ignored
      if key == 'hist_pick':
        hist_pickflag = kwargs['hist_pick']
      elif key == 'hist_log':
        hist_logflag = kwargs['hist_log']
      elif key == 'hist_cumul':
        hist_cumulflag = kwargs['hist_cumul']
      elif key == 'hist_type':
        hist_typeflag = kwargs['hist_type']
      elif key == 'hist_combined':
        hist_combinedflag = kwargs['hist_combined']
      elif key == 'hist_normed':
        hist_normedflag = kwargs['hist_normed']

    npick_arr = np.empty([len(self.events.keys())])
    n = 0
    for k in self.events.keys():
      ev = self.events[k]
      n_ph,n_p,n_s = ev.get_sum_picks()
      if hist_pickflag == 'P':
        npick_arr[n] = n_p
      elif hist_pickflag == 'S':
        npick_arr[n] = n_s     
      n += 1 

    print str(len(self.events.keys())) + ' events found in catalog'
    print 'Max number of picks= ' + str(np.amax(npick_arr))
    print 'Min number of picks = ' + str(np.amin(npick_arr))

    if hist_pickflag == 'P':
      max_weight = 2
      fac=2
    elif hist_pickflag == 'S':
      max_weight = 2
      fac=2
    max_bin = int(max_weight*(int(np.amax(npick_arr)/max_weight)+1))

    if len(bins) != 3:
      bins = [0, max_bin, max_weight]
    nbin = np.arange(bins[0], bins[1], bins[2])       

    if len(xlimits) != 2:
      xlimits = [npick_arr.min(), max_bin]

    if len(xticks) != 3:
      xticks = [xlimits[0], xlimits[1], max_weight*fac]

    ticks = np.arange(xticks[0], xticks[1], xticks[2])
    hrange=(xlimits[0], xlimits[1])

    hy, hx, patches = plt.hist(npick_arr,nbin,range=hrange,normed=hist_normedflag,log=hist_logflag,cumulative=hist_cumulflag,histtype=hist_typeflag,facecolor='green',alpha=0.75)
    if hist_normedflag:
      ymin = 0
      ymax = hy.max()
      ylabel_str='Density'
      hist_str='normalized'      
    if hist_logflag:
      ymin = 1
      ymax = hy.max() + 1000
      ylabel_str = 'Log(Nev)'
      hist_str = 'counts_log'
    if hist_logflag and (hist_cumulflag or hist_cumulflag==-1):
      ymin = 1
      ymax = hy.max() + 1000
      ylabel_str = 'Log(Nev)'
      hist_str = 'counts_log_cumul'
    if not hist_normedflag and not hist_logflag and not hist_cumulflag and hist_cumulflag!=-1: 
      ymin = 0
      ymax = hy.max() + 100
      ylabel_str = 'Number of events'
      hist_str = 'counts'
    if hist_combinedflag: 
      plt.hist(npick_arr,nbin,range=hrange,log=True,cumulative=False,histtype='bar',facecolor='green',alpha=0.75)
      ymin = 1
      ymax = hy.max() + 1000
      ylabel_str='Log(Nev)'
      hist_str='combined'    

    plt.axis([xlimits[0], xlimits[1], ymin, ymax])
    plt.xticks(ticks)
    title_str = 'Histogram of picks number per event: ' + hist_pickflag + '-picks'
    plt.title(title_str)
    plt.xlabel('Number of picks')
    plt.ylabel(ylabel_str)
    plt.grid(True)
    f_out= outpath + '/hist_npicks_' + hist_pickflag + '-picks_' + hist_str + '.pdf'
    plt.savefig(f_out)
    plt.close()


