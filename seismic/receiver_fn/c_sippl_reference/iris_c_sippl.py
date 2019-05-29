#module iris
#coding=utf-8

import obspy, os, glob, basic_c_sippl
from obspy.core import UTCDateTime,Stream,read
from obspy.clients.fdsn import Client
#from obspy.taup import taup
from obspy.geodetics.base import gps2dist_azimuth
from copy import deepcopy as cp

def read_ISC(infile,cat=False):
  """
  read ISC catalog file, hand over to raytracer (dict)
  other catalog format might have to be covered as well..
  """
  inf = open(infile,'r')
  data = inf.readlines()
  inf.close()

  if cat:
    bigdict = basic_c_sippl.Catalog() 
    for line in data:
      try:
        date,time,lat,lon,dep,mag = line.strip('\n').split(None)
        magflag = True
      except:
        date,time,lat,lon,dep = line.strip('\n').split(None)
        magflag = False
      yr,mn,dy = date.split('/')
      hr,min,sec = time.split(':')
      orig_time = UTCDateTime(int(yr),int(mn),int(dy),int(hr),int(min),float(sec))
      ev = basic_c_sippl.Event(orig_time,float(lat),float(lon),float(dep))
      ev.Picks = {}
      ev.info['magnitude'] = float(mag)
      ev.info['ident'] = 'NN'
      bigdict.events[orig_time.isoformat()] = cp(ev)


  else:
    counter = 1
    bigdict = {}

    for line in data:
      try:
        date,time,lat,lon,dep,mag = line.strip('\n').split(None)
        magflag = True
      except:
        date,time,lat,lon,dep = line.strip('\n').split(None)
        magflag = False
      yr,mn,dy = date.split('/')
      hr,min,sec = time.split(':')
      bigdict[counter] = {}
      bigdict[counter]['orig_time'] = UTCDateTime(int(yr),int(mn),int(dy),int(hr),int(min),float(sec))
      bigdict[counter]['event_lon'] = float(lon)
      bigdict[counter]['event_lat'] = float(lat)
      bigdict[counter]['event_depth'] = float(dep.strip('df'))
      if not magflag:
        bigdict[counter]['magnitude'] = 3.0
      else:
        bigdict[counter]['magnitude'] = float(mag)

      counter += 1

  return bigdict

def get_ttimes(indict,outdir,stat_lat,stat_lon,mode='IRIS',files='day_mseed',phase='P',stat='KMBL',netwrk='AU',interval=[-60,120]):
  """
  calculate travel times to phase, define time interval around these
  adds entries t1 and t2 to dict, also writes out new catalog with baz and takeoff angle

  stat_lat,stat_lon: station coordinates
  mode - IRIS for data retrieval from IRIS, for local data enter datapath in here
  files - only valid if mode = local, enter 'hr_mseed; for miniSEED hour files, 'day_mseed' for miniSEED day files and 'day_sac' for daily SAC files
  phase: phase name that is to be isolated
  interval: seconds before (negative!) and after the theoretical arrival that should be cut out
  writes data into outdir (event-wise folders) + generates a file with event information there
  mode - if SAC, set some SAC header variables
  """

  cutdict = {}

  outf = open(outdir+'/eventinfo','w')
  eid = 1

  for j in indict.keys():
    #calculate delta  
    distm = gps2dist_azimuth(indict[j]['event_lat'],indict[j]['event_lon'],stat_lat,stat_lon)[0]
    delta = distm / 111195. #XXX #seis.geo.delazi(indict[j]['event_lat'],indict[j]['event_lon'],stat_lat,stat_lon)[0]    
    if delta < 30. or delta>95.:
      continue
    try:
      ttdic = 0#obspy.taup.taup.getTravelTimes(delta,indict[j]['event_depth'],model='ak135')
    except AttributeError:
      print 'Error in taup'
      continue
    phaseflag = False
    for entry in ttdic:
      if entry['phase_name'] == phase:
        ttime = entry['time']
        tof = entry['take-off angle']
        phaseflag = True
    if not phaseflag:
      continue
    arrival = indict[j]['orig_time'] + ttime
    t1 = arrival + interval[0]
    t2 = arrival + interval[1]
    
    line = '%5i %04i %02i %02i %02i:%02i:%05.2f %8.4f %9.4f %6.2f %4s %07.2f %05.2f\n' % (int(eid),indict[j]['orig_time'].year, indict[j]['orig_time'].month, indict[j]['orig_time'].day, indict[j]['orig_time'].hour, indict[j]['orig_time'].minute, round(indict[j]['orig_time'].second + (indict[j]['orig_time'].microsecond / 10e6),3), indict[j]['event_lat'], indict[j]['event_lon'], indict[j]['event_depth'], phase, ttime, tof)
    #create directory, download waveform data
    dirname = outdir+'/'+str(eid)+'__'+str(indict[j]['orig_time'].year)+str(indict[j]['orig_time'].month)+str(indict[j]['orig_time'].day)+'_'+str(indict[j]['orig_time'].hour)+str(indict[j]['orig_time'].minute)+str(indict[j]['orig_time'].second)
    os.system("mkdir "+dirname)	
    os.chdir(dirname)
    
    if mode == 'IRIS':
      client = Client("IRIS")
      
      try:
        tr_z = client.get_waveforms(netwrk, stat,'10', 'BHZ', t1, t2) 
        tr_e = client.get_waveforms(netwrk, stat,'10', 'BH1', t1, t2) 
        tr_n = client.get_waveforms(netwrk, stat,'10', 'BH2', t1, t2) 
        tr_z.write(netwrk+'.'+stat+'.BHZ.mseed','MSEED') 
        tr_n.write(netwrk+'.'+stat+'.BHN.mseed','MSEED')
        tr_e.write(netwrk+'.'+stat+'.BHE.mseed','MSEED')

      except Exception:
        os.chdir("..")
        os.system("rm -rf "+dirname)
        continue
    else: #local data
      #find data and cut it out
      if files == 'hr_mseed':
        tracename_1 = '%4s%02i%02i%02i%02i' % (stat,indict[j]['orig_time'].year-2000,indict[j]['orig_time'].month,indict[j]['orig_time'].day,indict[j]['orig_time'].hour)
        print tracename_1
        if not indict[j]['orig_time'].hour == 23:
          tracename_2 = '%4s%02i%02i%02i%02i' % (stat,indict[j]['orig_time'].year-2000,indict[j]['orig_time'].month,indict[j]['orig_time'].day,indict[j]['orig_time'].hour+1)
        else:
          tracename_2 = '%4s%02i%02i%02i00' % (stat,indict[j]['orig_time'].year-2000,indict[j]['orig_time'].month,indict[j]['orig_time'].day+1)
        try:
          julday = '%03i' % (indict[j]['orig_time'].julday)
          yr = str(indict[j]['orig_time'].year)
          print mode+'/'+yr+'/'+julday+'/'+tracename_1
          print glob.glob(mode+'/'+yr+'/'+julday+'/'+tracename_1+'*.[H,B,S]HE')[0]
          #in case there's more than one trace merge them
          etrace = Stream()
          for i in glob.glob(mode+'/'+yr+'/'+julday+'/'+tracename_1+'*.[H,B,S]HE'):
            etrace += read(i)
          for j in glob.glob(mode+'/'+yr+'/'+julday+'/'+tracename_2+'*.[H,B,S]HE'):
            etrace += read(j)
          etrace.merge(method=-1,fill_value='interpolate')
          ntrace = Stream()
          for i in glob.glob(mode+'/'+yr+'/'+julday+'/'+tracename_1+'*.[H,B,S]HN'):
            ntrace += read(i)
          for j in glob.glob(mode+'/'+yr+'/'+julday+'/'+tracename_2+'*.[H,B,S]HN'):
            ntrace += read(j)
          ntrace.merge(method=-1,fill_value='interpolate')
          ztrace = Stream()
          for i in glob.glob(mode+'/'+yr+'/'+julday+'/'+tracename_1+'*.[H,B,S]HZ'):
            ztrace += read(i)
          for j in glob.glob(mode+'/'+yr+'/'+julday+'/'+tracename_2+'*.[H,B,S]HZ'):
            ztrace += read(j)
          ztrace.merge(method=-1,fill_value='interpolate')

        except:
          os.chdir("..")
          os.system("rm -rf "+dirname)
          continue
      elif files == 'day_mseed':
        tracename_base = '%02i%02i%02i' % (int(str(indict[j]['orig_time'].year)[2:]),indict[j]['orig_time'].month,indict[j]['orig_time'].day)
        print mode+'/'+stat+'/'+stat+tracename_base+'*.HHE'
        print glob.glob(mode+'/'+stat+'/'+stat+tracename_base+'*.[H,B]HE')
        
        try:
          etr = read(glob.glob(mode+'/'+stat+'/'+stat+tracename_base+'*.[H,B]HE')[0])
          etr.merge(fill_value='interpolate')
          print len(etr)
          etrace = etr[0]
          ntr = read(glob.glob(mode+'/'+stat+'/'+stat+tracename_base+'*.[H,B]HN')[0])
          ntr.merge(fill_value='interpolate')
          ntrace = ntr[0]
          ztr = read(glob.glob(mode+'/'+stat+'/'+stat+tracename_base+'*.[H,B]HZ')[0])
          ztr.merge(fill_value='interpolate')
          ztrace = ztr[0]
        except Exception:
          os.chdir("..")
          os.system("rm -rf "+dirname)
          continue

      elif files == 'geofon':

        jday = '%03i' % (indict[j]['orig_time'].julday)

        print mode+'/'+str(indict[j]['orig_time'].year)+'/'+netwrk+'/'+stat+'/*/*.'+jday

        try:
          etr = read(glob.glob(mode+'/'+str(indict[j]['orig_time'].year)+'/'+netwrk+'/'+stat+'/HHE.D/*.'+jday)[0])
          etr.merge(fill_value='interpolate')
          print len(etr)
          etrace = etr[0]
          ntr = read(glob.glob(mode+'/'+str(indict[j]['orig_time'].year)+'/'+netwrk+'/'+stat+'/HHN.D/*.'+jday)[0])
          ntr.merge(fill_value='interpolate')
          ntrace = ntr[0]
          ztr = read(glob.glob(mode+'/'+str(indict[j]['orig_time'].year)+'/'+netwrk+'/'+stat+'/HHZ.D/*.'+jday)[0])
          ztr.merge(fill_value='interpolate')
          ztrace = ztr[0]
        except Exception:
          os.chdir("..")
          os.system("rm -rf "+dirname)
          continue

      elif files == 'day_sac':

        tracename_base = '%04i.%03i' % (indict[j]['orig_time'].year,indict[j]['orig_time'].julday)
        print glob.glob(mode+'/'+stat+'/'+tracename_base+'*.BHE.SAC')

        try:
          etrace = read(glob.glob(mode+'/'+stat+'/'+tracename_base+'*.BHE.SAC')[0],format='SAC')[0]
          ntrace = read(glob.glob(mode+'/'+stat+'/'+tracename_base+'*.BHN.SAC')[0],format='SAC')[0] 
          ztrace = read(glob.glob(mode+'/'+stat+'/'+tracename_base+'*.BHZ.SAC')[0],format='SAC')[0]

        except Exception:
          os.chdir("..")
          os.system("rm -rf "+dirname)
          continue

        etrace.stats.starttime = UTCDateTime(indict[j]['orig_time'].year,indict[j]['orig_time'].month,indict[j]['orig_time'].day,0,0,0)
        ntrace.stats.starttime = UTCDateTime(indict[j]['orig_time'].year,indict[j]['orig_time'].month,indict[j]['orig_time'].day,0,0,0)
        ztrace.stats.starttime = UTCDateTime(indict[j]['orig_time'].year,indict[j]['orig_time'].month,indict[j]['orig_time'].day,0,0,0)
      
      e_event = etrace.slice(t1,t2)
      n_event = ntrace.slice(t1,t2)
      z_event = ztrace.slice(t1,t2)      

      try:
        e_event.write(netwrk+'.'+stat+'.BHE.mseed','MSEED')
        n_event.write(netwrk+'.'+stat+'.BHN.mseed','MSEED')
        z_event.write(netwrk+'.'+stat+'.BHZ.mseed','MSEED')
      except:
        os.chdir("..")
        os.system("rm -rf "+dirname)
        continue
      """
        #if len(e_event) == 0:
        #try:
        #  t1.year -= 1
        # t2.year -= 1
        #  e_event = etrace.slice(t1,t2)
        #  e_event.write(netwrk+'.'+stat+'.BHE.mseed','MSEED')
        # n_event = ntrace.slice(t1,t2)
        # n_event.write(netwrk+'.'+stat+'.BHN.mseed','MSEED')
        # z_event = ztrace.slice(t1,t2)
        # z_event.write(netwrk+'.'+stat+'.BHZ.mseed','MSEED')
        #except:
        #  continue
      """
    outf.write(line)
    os.chdir("..")
    eid += 1 

def get_bulkdata(stat,netwrk,loc_code,t1,t2,outdir,ident='BH'):
  """
  bulk data download from IRIS fdsn server
  ident: BH, HH, SH, LH etc.
  """

  client = Client("IRIS")
  while t1 < t2:
    print t1
    try:
      tr_e = client.get_waveforms(netwrk,stat,loc_code,ident+'E',t1,t1+86400.)
      tr_n = client.get_waveforms(netwrk,stat,loc_code,ident+'N',t1,t1+86400.)
      tr_z = client.get_waveforms(netwrk,stat,loc_code,ident+'Z',t1,t1+86400.)
      print 'Data found'
    except:
      print 'No data'
      t1 += 86400.
      continue
     

    string = '%02i%02i%02i' % (int(str(t1.year)[2:]),t1.month,t1.day)
    #write to outdir
    tr_e.write(outdir+'/'+stat+'/'+stat+string+'000000.'+ident+'E','MSEED')
    tr_n.write(outdir+'/'+stat+'/'+stat+string+'000000.'+ident+'N','MSEED')
    tr_z.write(outdir+'/'+stat+'/'+stat+string+'000000.'+ident+'Z','MSEED')

    t1 += 86400.    


def ttime_dict(ev_lat,ev_lon,ev_dep,s_lat,s_lon,phase='P'):
  """
  queries raytracing results from IRIS' traveltime webservice
  """
  from obspy.clients.iris.client import Client as Cl

  client = Cl()
  result = client.traveltime(evloc=(ev_lat,ev_lon),staloc=(s_lat,s_lon),evdepth=ev_dep)
  lines = result.split('\n')
  phasedict = {}
  for j in lines:
    try:
      if j.split(None)[2] == phase:
        phasedict['name'] = phase
        phasedict['dist'] = float(j.split(None)[0])
        phasedict['ttime'] = float(j.split(None)[3])
        phasedict['rayp'] = float(j.split(None)[4])
        phasedict['tof'] = float(j.split(None)[5])
        phasedict['incidence'] = float(j.split(None)[6])
    except IndexError:
      continue
  if len(phasedict) == 0:
    print "Phase not found"
  return phasedict
