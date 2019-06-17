#!/usr/bin/env python
# coding=utf-8
# module HK_stacking
"""
Convenience and plotting routines for HK-stacking code.  Requires custom built binary executable `hk2d`. (Provenance: TBD)

This code used for this publication by Christian Sippl:
  https://www.sciencedirect.com/science/article/pii/S0040195116300245
"""

# pylint: skip-file

from numpy import *
import glob, os, ausrem, util, cPickle, iris
from pylab import *
#from obspy.sac import SacIO
from obspy.io.sac.sactrace import SACTrace
from obspy.core.util.geodetics import gps2DistAzimuth

def HK_stack(stat,vp=6.5,H=[15.,100.,0.05],K=[1.5,2.0,0.001],weight=[0.33,0.33,0.33],smooth=0.1,n_stack=4,X=0,baz=[0,360],indir='.'):
  """
  call hk2d with input parameters
  """

  wgt = str(weight[0])+'_'+str(weight[1])+'_'+str(weight[2])
  outdir = get_folder_struct(stat,wgt)

  #calculate HK stack

  cmd = 'hk2d -P '+str(vp)+' -Z '+str(H[0])+'/'+str(H[1])+'/'+str(H[2])+' -K '+str(K[0])+'/'+str(K[1])+'/'+str(K[2])+' -W '+str(weight[0])+'/'+str(weight[1])+'/'+str(weight[2])+' -T '+str(smooth)+' -N '+str(n_stack)+' -X '+str(X)+' -R '+str(baz[0])+'/'+str(baz[1])+' -I '+indir+' -L '+indir+'/'+stat+'__allRFs.list -D '+outdir+' -O sfr2d.AH.ANQ/hkr2d.AH.ANQ -o'
  os.system(cmd)

  #plot results into folder
  hkdict = read_outfiles(outdir)
  plot_sfunc(hkdict,wgt,stat,outf_dir=outdir+'/plots')  


def get_folder_struct(stat,weight,rootf='.'):
  """
  set up output structure before running hk2d
  """
 
  if not rootf == '.':
    cmd = 'mkdir -p '+rootf+'/'+stat+'/'+weight
  else:
    cmd = 'mkdir -p '+stat+'/'+weight
  os.system(cmd)
  os.system(cmd+'/plots')

  outdir = cmd.split(None)[2]
  return outdir


def read_outfiles(outpath):
  """
  read output files into dictionary format
  returns: hkdict (keys = settings for X)
  """
  hkdict = {}  

  for j in range(10):
    hkdict[j] = {}
    hkfiles = glob.glob(outpath+'/hkr2d*.x'+str(j)+'.??') 
    if len(hkfiles) == 0:
      continue
    else:
      for q in hkfiles:
        num = q.split('.')[-1]
        fl = open(q,'r')
        data = fl.readlines()
        fl.close()  
    
        hkdict[j]['vp_avg'] = data[0].split(None)[-2]
    
        hkdict[j]['H_inc'] = float(data[1].split(None)[-4])
        hkdict[j]['H_min'] = float(data[1].split(None)[3])
        hkdict[j]['H_max'] = float(data[1].split(None)[5])
        hkdict[j]['H_steps'] = int(round((float(hkdict[j]['H_max']) - float(hkdict[j]['H_min']))/float(hkdict[j]['H_inc'])+1,1))  
    
        hkdict[j]['K_inc'] = float(data[2].split(None)[-4])
        hkdict[j]['K_min'] = float(data[2].split(None)[3])
        hkdict[j]['K_max'] = float(data[2].split(None)[5])
        hkdict[j]['K_steps'] = int(round((float(hkdict[j]['K_max']) - float(hkdict[j]['K_min']))/float(hkdict[j]['K_inc'])+1,1))
    
        hkdict[j]['weights'] = [float(data[3].split(None)[3].strip(',')), float(data[3].split(None)[6].strip(',')), float(data[3].split(None)[-1].strip('\n'))]
        hkdict[j]['smooth_win'] = data[4].split(None)[-1].strip('\n')
        hkdict[j]['N_root'] = data[5].split(None)[-1].strip('\n')
        hkdict[j]['baz'] = [float(data[7].split(None)[3]),float(data[7].split(None)[5].strip('\n'))]
        hkdict[j]['nrf'] = int(data[8].split(None)[1])
    
        hkdict[j]['H_opt_'+num] = float(data[16].split(None)[3])
        hkdict[j]['H_uncert_'+num] = float(data[16].split(None)[5])
        hkdict[j]['K_opt_'+num] = float(data[17].split(None)[3])
        hkdict[j]['K_uncert_'+num] = float(data[17].split(None)[5])
    
        hkdict[j]['sigma_s'] = float(data[15].split(None)[2].strip('\n')) 
        hkdict[j]['a_'+num] = float(data[18].split(None)[2].strip(','))
        hkdict[j]['b_'+num] = float(data[18].split(None)[5].strip(','))
        hkdict[j]['alpha_'+num] = float(data[18].split(None)[-1].strip('\n'))
    
        hkdict[j]['matrix'] = zeros([hkdict[j]['K_steps'],hkdict[j]['H_steps']])
    
    try:
      sfl = open(glob.glob(outpath+'/sfr2d*.x'+str(j))[0],'r')
    except:
      continue
    sdata = sfl.readlines()
    sfl.close()
    for k in sdata[3:]:
      K,H,sval = k.strip('\n').split(None) 
      index_K = int(round((float(K) - hkdict[j]['K_min']) / hkdict[j]['K_inc'],1))
      index_H = hkdict[j]['H_steps'] - (int(round((float(H) - hkdict[j]['H_min']) / hkdict[j]['H_inc'],1))) - 1
      hkdict[j]['matrix'][index_K][index_H] = float(sval)

  return hkdict


def get_hist(matrix):
  """
  get vertical and horizontal histogram from sfunction matrix
  """

  #rows first (K)
  ln = shape(matrix)[0]
  wd = shape(matrix)[1]
  ln_ar = zeros(ln)
  wd_ar = zeros(wd)
  for k in range(len(ln_ar)):
    ln_ar[k] = matrix[k].sum()
   
  #now the columns
  wd_ar = zeros(wd)
  for u in range(len(wd_ar)):
    wd_ar[u] = matrix[:,u].sum()   
   
  return ln_ar,wd_ar


def plot_sfunc(hkdict,wgt,stat,outf_dir='.'):
  """
  plots sfunc, optional: maxima indication, side histograms,...
  """

  for j in hkdict.keys():

    if not hkdict[j].has_key('matrix'):
      continue
    figure(figsize=(12,12))

    #get axes configuration
    left, width = 0.07, 0.75
    bottom, height = 0.07, 0.75
    bottom_h = left_h = left + width + 0.01

    main_ax = [left,bottom,width,height]
    histx_ax = [left,bottom_h,width,0.12]
    histy_ax = [left_h,bottom,0.12,height]

    axMain = axes(main_ax)
  
    x_ar_depth = arange(hkdict[j]['K_min'],hkdict[j]['K_max']+hkdict[j]['K_inc'],0.1)
    y_ar_depth = arange(hkdict[j]['H_max'],hkdict[j]['H_min']-hkdict[j]['H_inc'],-5)
 
    x_ar_num = arange(0,hkdict[j]['K_steps']+1,(1/float(hkdict[j]['K_inc']))*0.1)
    y_ar_num = arange(0,hkdict[j]['H_steps']+1,(1/float(hkdict[j]['H_inc']))*5)

    axMain.imshow(hkdict[j]['matrix'].transpose(),aspect='auto',interpolation='none')

    xticks(x_ar_num,x_ar_depth)
    xlabel('Vp/Vs ratio')
    yticks(y_ar_num,y_ar_depth)
    ylabel('Moho depth [km]')  
 
    #contouring levels
    V = [0.1,0.25,0.5,0.75,0.85,0.9,0.95]
 
    axMain.contour(hkdict[j]['matrix'].transpose(),V,colors='black',linewidths=[0.8,0.8,2,0.8,2,0.8,2,0.8])
    
    x_opt_list = []
    y_opt_list = []
    for a in range(9):
      if hkdict[j].has_key('K_opt_0'+str(a)):
        x_opt_list.append(hkdict[j]['K_opt_0'+str(a)])   
        y_opt_list.append(hkdict[j]['H_opt_0'+str(a)])
        

    for qua in range(len(x_opt_list)):
      x_now = (x_opt_list[qua] - hkdict[j]['K_min']) / hkdict[j]['K_inc']
      y_now = hkdict[j]['H_steps'] - ((y_opt_list[qua] - hkdict[j]['H_min']) / hkdict[j]['H_inc'])
      axMain.plot(x_now,y_now,'ks',markersize=7)
      axMain.plot([0,x_now],[y_now,y_now],'w--')
      axMain.plot([x_now,x_now],[hkdict[j]['H_steps'],y_now],'w--') 

      axMain.text(0,y_now-4,'H = '+str(y_opt_list[qua]),color='red',fontsize=14)    
      axMain.text(x_now+4,hkdict[j]['H_steps'] - 160,'K = '+str(x_opt_list[qua]),color='red',fontsize=14,rotation='vertical')

    #title('Contour plot sfunction, X='+str(j))

    #side histograms??
    ln_ar, wd_ar = get_hist(hkdict[j]['matrix'])
  
    ln_ar /= ln_ar.max() #normalization
    wd_ar /= wd_ar.max()

    axHistX = axes(histx_ax)
    axHistX.plot(ln_ar,'r--')
    yticks([0,0.25,0.5,0.75,1.])
    xticks([])
    ylabel('normalized amplitude')
    ylim([0,1.03])

    y_ar = arange(1100,-1,-1)
    axHistY = axes(histy_ax)
    axHistY.plot(wd_ar,y_ar,'r--')
    ylim([0,1100])
    xlim([0,1.03])
    yticks([])
    xticks([0,0.25,0.5,0.75,1],rotation=45)
    xlabel('normalized amplitude')

    outf_name = outf_dir+'/'+stat+'_'+wgt+'_'+str(j)+'.pdf'
    savefig(outf_name)


def get_ausrem_avg(stat,info_file='/home/christian/info_WOMBAT'):
  """
  retrieve average vp from ausrem
  """
  lat,lon,elev = util.get_stat_coords(stat,info_file)
  vp = ausrem.read_infiles(typ='vp')
  v,d = ausrem.get_1D_prof(vp,lat,lon,dep=[0,60],interpolate=False)

  #sum up values until vp>7.5
  vp_all = 0
  count = 0
  for i in range(len(v)):
    if float(v[i]) < 7.5:
      vp_all += v[i] 
      count += 1

  avg_vp = vp_all / float(count) 
  return avg_vp


def automat_HK(flder):
  """
  piece of advice: use absolute path for flder!!
  """

  diclist = glob.glob(flder+'/*[rR][fF]_sel*.bin')
  #turn dicts to (temporary) SAC files
  """
  try:
    os.mkdir(flder+'/temp')
  except OSError: #already exists
    pass
  """
  write_files_from_dict(diclist,outfolder=flder+'/HK_stacks_new')
  
  os.chdir(flder+'/HK_stacks_new')
  #works so far...
  #now call HK stacking for each station (initial run...determination of raw matrix)
  for item in diclist:
    stat = item.split('/')[-1].split('_')[0]  
    print(stat)
    call_HK(stat,dirct=flder+'/HK_stacks_new',X=7,info_file='/home/christian/info_file')
  #read everything in and set up giant dict
  dic = dict_setup(info_file='/home/christian/info_file')
  cPickle.dump(dic,open('te','wb'))
  """
  dic = cPickle.load(open('te','rb'))
  dic_new = pick_max(dic)
  nit = range(5)
  for i in nit:
    if i == 0:
      #do initial weighting, pick max again
      dic_2 = mask_primary(dic_new,H_coarse=[25.,55.],H_fine=[30.,45.],K_coarse=[1.6,1.95],K_fine=[1.65,1.85])
      dic_2_new = pick_max(dic_2,apply_mask=True)
  return dic,dic_new,dic_2_new
  """
  #read_outfiles()
  #now the magic has to happen


def dict_setup(flder='.',X=7,info_file='/home/christian/info_file'):

  #get list of stations
  fld_lst = glob.glob(flder+'/????')
  fld_lst.append(flder+'/AUKAL')
  bigdict = {}

  statlst = []
  for i in fld_lst:
    stat = i.split('/')[-1]
    statlst.append(stat) #list of all utilized stations
  for j in fld_lst:
    mtrx = read_outfiles(j+'/0.33_0.33_0.33')[X]
    sta = j.split('/')[-1]
    print(sta)
    bigdict[sta] = {}
    bigdict[sta]['matrix'] = mtrx['matrix']
    bigdict[sta]['header'] = {} #store stuff like dimensions, translation to vp/vs and thickness values etc.
    bigdict[sta]['header']['K_min'] = mtrx['K_min']
    bigdict[sta]['header']['K_max'] = mtrx['K_max']
    bigdict[sta]['header']['H_max'] = mtrx['H_max']
    bigdict[sta]['header']['H_min'] = mtrx['H_min']
    bigdict[sta]['header']['K_steps'] = mtrx['K_steps']
    bigdict[sta]['header']['H_steps'] = mtrx['H_steps']
    bigdict[sta]['header']['H_inc'] = mtrx['H_inc']
    bigdict[sta]['header']['K_inc'] = mtrx['K_inc']
    bigdict[sta]['mask'] = ones(shape(bigdict[sta]['matrix']))
    #inter-station distances
    bigdict[sta]['distances'] = {}
    lat_s,lon_s,el_s = util.get_stat_coords(sta,info_file)
    for k in statlst:
      if not k == sta: #that's obviously 0
        lat1,lon1,el1 = util.get_stat_coords(k,info_file)
        distkm = round(gps2DistAzimuth(lat1,lon1,lat_s,lon_s)[0] / 1000.,3)
        bigdict[sta]['distances'][k] = distkm 

  return bigdict


def pick_max(dic,apply_mask=False):
  """
  """
  dic_n = dic.copy()
  for sta in dic_n.keys():
    maxlist = []
    hklist = []  

    if apply_mask:
      dic_n[sta]['matrix'] *= dic_n[sta]['mask']
    mtx = dic_n[sta]['matrix']
    
    mx = mtx.max()
    k,h = shape(mtx)
    for i in range(k):
      for j in range(h):
        if mtx[i,j] == mx:
          maxlist.append([i,j]) 
          hklist.append([dic_n[sta]['header']['K_min'] + i*dic_n[sta]['header']['K_inc'], dic_n[sta]['header']['H_max'] - j*dic_n[sta]['header']['H_inc']])
    dic_n[sta]['header']['max_HK'] = hklist
    dic_n[sta]['header']['max_ij'] = maxlist   
  return dic_n


def mask_primary(bigdict,H_coarse=[25.,55.],H_fine=[30.,45.],K_coarse=[1.6,1.95],K_fine=[1.65,1.85]):
  """
  initialization of mask
  1 inside fine, 0.25 outside coarse, linear between
  """
  #get matrix positions of these  

  bigdict2 = bigdict.copy()

  for i in bigdict2.keys():

    mtrx = bigdict2[i]['mask'].copy()

    H_lst = array([bigdict2[i]['header']['H_max'] - H_coarse[1], H_coarse[1] - H_fine[1], H_fine[1] - H_fine[0], H_fine[0] - H_coarse[0], H_coarse[0] - bigdict2[i]['header']['H_min']])
    H_lst *= 1/bigdict2[i]['header']['H_inc'] 

    #same for K
    K_lst = array([K_coarse[0] - bigdict2[i]['header']['K_min'],K_fine[0] - K_coarse[0], K_fine[1] - K_fine[0], K_coarse[1] - K_fine[1], bigdict2[i]['header']['K_max'] - K_coarse[1]])
    K_lst *= 1/bigdict2[i]['header']['K_inc']
 
    #now first the outer shell
    mtrx_mod = mtrx.copy()
    mtrx_mod[0:K_lst[0]] = 0.25
    mtrx_mod[(K_lst.sum() - K_lst[-1]):] = 0.25
   
    mtrx_mod[:,0:H_lst[0]] = 0.25
    mtrx_mod[:,H_lst.sum()-H_lst[-1]:] = 0.25
   
    #now the linear interpolative one
    #start at H, upper end
    #calculate single step
    H_step_upper = (1-0.25) / float(H_lst[1]) 
    H_step_lower = (1-0.25) / float(H_lst[-2])

    K_step_lower = (1-0.25) / float(K_lst[1])
    K_step_upper = (1-0.25) / float(K_lst[-2])

    for ln in range(int(H_lst[0]),int(H_lst[0]+H_lst[1])):
      for zl in range(int(K_lst[0]),int(K_lst[0]+K_lst[1]+K_lst[2]+K_lst[3])):
        mtrx_mod[zl,ln] = (ln - H_lst[0])*H_step_upper + 0.25

    #now lower end of H
 
    for ln in range(int(H_lst[0] + H_lst[1] + H_lst[2]), int(H_lst[0]+H_lst[1]+H_lst[2]+H_lst[3])):
      for zl in range(int(K_lst[0]),int(K_lst[0]+K_lst[1]+K_lst[2]+K_lst[3])):
        mtrx_mod[zl,ln] = (int(H_lst[0]+H_lst[1]+H_lst[2]+H_lst[3]) - ln) * H_step_lower + 0.25

    #now upper end of K
    ql = range(int(H_lst[0]), int(H_lst[0]+H_lst[1]+H_lst[2]+H_lst[3]))
    for ln in ql:
      al = range(int(K_lst[0]+K_lst[1]+K_lst[2]),int(K_lst[0]+K_lst[1]+K_lst[2]+K_lst[3]))
      for zl in al:
        if mtrx_mod[zl,ln] == 1.: #not already modified
          mtrx_mod[zl,ln] = (int(K_lst[0]+K_lst[1]+K_lst[2]+K_lst[3]) - zl) * K_step_upper + 0.25
        #else:
        #  if (zl - int(K_lst[0]+K_lst[1]+K_lst[2])) / float(len(al)) > (ln - int(H_lst[0])) / float(len(ql)):
        #    mtrx_mod[zl,ln] = (int(K_lst[0]+K_lst[1]+K_lst[2]+K_lst[3]) - zl) * K_step_upper + 0.25

    #and the lower one...
    for ln in range(int(H_lst[0]), int(H_lst[0]+H_lst[1]+H_lst[2]+H_lst[3])):
      for zl in range(int(K_lst[0]),int(K_lst[0]+K_lst[1])):
        if mtrx_mod[zl,ln] == 1.: #not already modified
          mtrx_mod[zl,ln] = (zl - int(K_lst[0])) * K_step_lower + 0.25
        #else:
        #  mtrx_mod[zl,ln] = (mtrx_mod[zl,ln] + (zl - int(K_lst[0])) * K_step_lower + 0.25)/2.

    bigdict2[i]['mask'] = mtrx_mod
    
  return bigdict2


def change_mask(bigdict,rad_H=8.):
  """
  modify masking array
  """
  #get 4 closest stations
  for i in bigdict.keys():
    distlist = sorted(bigdict[i]['distances'].values())[0:4]
    klist = []
    hlist = []
    #retrieve their values for H and K
    for val in distlist:
      #get associated station
      for sta in bigdict[i]['distances'].keys():
        if bigdict[i]['distances'][sta] == val:
          k,h = bigdict[sta]['header']['max_ij'][0]
          klist.append(k)
          hlist.append(h)
          break

    #define region around these values, modify mask
    for u in range(len(hlist)):
      dist = rad_H / 0.05  
      k,h = shape(bigdict[i]['mask'])    
      bigdict[i]['mask'][:,0:(hlist[u]-dist)] *= 0.5
      bigdict[i]['mask'][:,(hlist[u]+dist):] *= 0.5
      for a in range(k):
        for b in range(h):
          if b > (hlist[u]-dist) and b < (hlist[u]+dist):
            if a > klist[u]+dist:
              bigdict[i]['mask'][a][b] *= 0.5
            elif a < klist[u]-dist:
              bigdict[i]['mask'][a][b] *= 0.5 

  return bigdict


def draw_nRFs(indict,stat,reps=100):
  """

  """
  from random import shuffle
  from copy import deepcopy as cp
  n = len(indict)
  for j in range(1,reps+1):
    os.system('mkdir '+str(j))
    dict_new = {}
    for k in range(n):
      lst = indict.keys()
      shuffle(lst)
      ky_chosen = lst[0]
      dict_new[k+1] = cp(indict[ky_chosen])

    write_files_from_dict(dict_new,stat,outfolder=str(j))


def write_files_from_dict(dat,stat,outfolder='.'):
  """
  write SAC files for RFs in indict list
  """
  num = 1
#  for dic in indict:
#    dat = cPickle.load(open(dic,'rb')) 
#    stat = dic.split('/')[-1].split('_')[0]
  outf = open(outfolder+'/'+stat+'__allRFs.list','w')
  for entry in dat.keys():
    name = '%03i__%04i%02i%02i__%02i%02i%02i' % (num,dat[entry]['orig_time'].year,dat[entry]['orig_time'].month,dat[entry]['orig_time'].day,dat[entry]['orig_time'].hour,dat[entry]['orig_time'].minute,dat[entry]['orig_time'].second)
    name_full = name+stat+'__RF.SAC'
    trc = dat[entry]['Traces']['RF']
    trc.write(outfolder+'/'+name_full,format='SAC')
    #now set header values
    hd = SACTrace.read(outfolder+'/'+name_full,headonly=True)
    hd.dist = dat[entry]['Distance']
    hd.az = dat[entry]['Azimuth']
    hd.baz = dat[entry]['Backazimuth']
    hd.evla = dat[entry]['event_lat']
    hd.evlo = dat[entry]['event_lon']
    hd.evdp = dat[entry]['event_depth']
    hd.write(outfolder+'/'+name_full,headonly=True)
    #and write list file for each station
    raytr = iris.ttime_dict(dat[entry]['event_lat'],dat[entry]['event_lon'],dat[entry]['event_depth'],dat[entry]['station_lat'],dat[entry]['station_lon'],phase='P')
    outf.write(name_full+'  '+str(round(raytr['rayp'],5))+'  0.00  '+str(round(dat[entry]['Backazimuth'],3))+'  '+str(round(raytr['dist'],3))+'\n')
    num += 1

  outf.close()
    
def call_HK(stat,vp,H=[30.,60.,0.05],K=[1.65,1.85,0.002],weight=[0.33,0.33,0.33],smooth=0.1,n_stack=4,X=0,baz=[0,360],dirct='.',info_file='/home/christian/info_WOMBAT'):
  """
  modified HK main code for automated usage
  """

  #get vp
  #vp = round(get_ausrem_avg(stat,info_file=info_file),3)

  wgt = str(weight[0])+'_'+str(weight[1])+'_'+str(weight[2])
  #outdir = get_folder_struct(stat,wgt)
  outdir = '.'

  print(outdir)
  #calculate HK stack

  cmd = 'hk2d -P '+str(vp)+' -Z '+str(H[0])+'/'+str(H[1])+'/'+str(H[2])+' -K '+str(K[0])+'/'+str(K[1])+'/'+str(K[2])+' -W '+str(weight[0])+'/'+str(weight[1])+'/'+str(weight[2])+' -T '+str(smooth)+' -N '+str(n_stack)+' -X '+str(X)+' -R '+str(baz[0])+'/'+str(baz[1])+' -I '+dirct+' -L '+dirct+'/'+stat+'__allRFs.list -D '+dirct+'/'+outdir+' -O sfr2d.AH.ANQ/hkr2d.AH.ANQ -o'
  os.system(cmd)


def plot_HK_map(flder='.',lonrange=[120.,125.],latrange=[-34.,-30.],info_file='/home/christian/info_file',z='H',
                gravity=False,smooth=0.1,eps=0.33,Ps=False,pval=False,save=False,mode='color'):
  """
  plots an interpolated map of crustal thickness (or vp/vs) from a bunch of HK-stacking derived values
  Ps=True only works with z='H'
  mode = 'color' or 'contour' or 'both' (the last only works if gravity=False)
  """
  #get dictionary, pick maxima from matrices
  if Ps:
    #read values from file, multiply with 8
    xvals = []
    yvals = []
    zvals = []
    if pval:
      infile = open('all_checked','r')
    else:
      infile = open('Ps_all','r')
    dtt = infile.readlines()
    infile.close()

    for j in dtt:
      if pval:
        sta,stat_lat,stat_lon,H,K,dum = j.strip('\n').split(None)
        print(sta)
      else:
        sta,Ps = j.strip('\n').split(None)
      #stat_lat,stat_lon,stat_elev = util.get_stat_coords(sta,info_file)
      
      if float(stat_lat) >= latrange[0] and float(stat_lat) <= latrange[1] and float(stat_lon) >= lonrange[0] and float(stat_lon) <= lonrange[1]:
        if z == 'K' and K == '----':
          continue
        xvals.append(float(stat_lon))
        yvals.append(float(stat_lat))
        if pval:
          if z == 'H':
            zvals.append(float(H))
          elif z == 'K':
            zvals.append(float(K))
        else:
          zvals.append(float(Ps)*8)

  else:
    dic = dict_setup(flder=flder)
    dic_picked = pick_max(dic)

    #cPickle.dump(dic_picked,open('storage_dic','wb'))
    #dic_picked = cPickle.load(open('storage_dic','rb'))

    #now station-wise...get coordinates, store in 4-tuples (lat,lon,H,K)
    xvals = []
    yvals = []
    zvals = []
    for stat in dic_picked.keys():
      lat,lon,elev = util.get_stat_coords(stat,info_file)
      H = dic_picked[stat]['header']['max_HK'][0][1]
      K = dic_picked[stat]['header']['max_HK'][0][0]
      if float(lat) >= latrange[0] and float(lat) <= latrange[1] and float(lon) >= lonrange[0] and float(lon) <= lonrange[1]:
        xvals.append(lon)
        yvals.append(lat)
        if z == 'H':
          zvals.append(H)
        elif z == 'K':
          zvals.append(K)  


  #set control nodes (at 40 km)
  for x in arange(lonrange[0],lonrange[1],0.25):
    for y in [latrange[0]-0.25,latrange[1]+0.25]:
      xvals.append(x)
      yvals.append(y)
      if z =='H':
        zvals.append(37.5)
      elif z == 'K':
        zvals.append(1.73)

  for y in arange(latrange[0],latrange[1],0.25):
    for x in [lonrange[0] - 0.25, lonrange[1]+0.25]:
      xvals.append(x)
      yvals.append(y)
      if z =='H':
        zvals.append(37.5)
      elif z == 'K':
        zvals.append(1.73)   

  print(xvals, yvals)


  yi = linspace(latrange[0],latrange[1],500)
  xi = linspace(lonrange[0],lonrange[1],500)
  XI,YI = meshgrid(xi,yi)

  from mpl_toolkits import basemap
  figure(figsize=(12,12))
  m = basemap.Basemap(projection='merc',llcrnrlon=lonrange[0],llcrnrlat=latrange[0],urcrnrlon=lonrange[1],urcrnrlat=latrange[1],resolution='h')

  XImap,YImap = m(XI,YI)

  from scipy.interpolate import Rbf

  xvals_map,yvals_map = m(xvals,yvals)
  rbf = Rbf(xvals,yvals,zvals,epsilon=eps,smooth=smooth)
  #rbf = Rbf(xvals_map,yvals_map,zvals)
  ZI = rbf(XI,YI)
 
  outfile = open('Moho_lowres.txt','w')
  for i in range(len(ZI)):
    for j in range(len(ZI[0])):
      outfile.write(str(round(XI[i][j],3))+' '+str(round(YI[i][j],3))+' '+str(round(ZI[i][j],3))+'\n')
 
  m.drawcoastlines()
  #ylim([latrange[0],latrange[1]])
  #xlim([lonrange[0],lonrange[1]])
  if gravity:
    #read in gravity grid, plot contours on top
    gravdat = loadtxt('/home/sippl/sandbox/gravity.txt')

    #convert to regular grid
    gravxi = arange(120.,128.01,0.01)
    gravyi = arange(-34.,-30.01,0.01)
    gravzi = griddata(gravdat[:,0],gravdat[:,1],gravdat[:,2],gravxi,gravyi,interp='linear')

    gravlon,gravlat = np.meshgrid(gravxi,gravyi)
    gravx,gravy = m(gravlon,gravlat)

    m.contour(gravx,gravy,gravzi,[-800,-400,-200,0,150,300],colors=['b','b','r','r','r','r'],linestyles='-')
  if mode == 'color':
    m.pcolor(XImap,YImap,ZI,cmap='jet',vmin=32,vmax=52)
  
  elif mode == 'contour':
    m.contour(XImap,YImap,ZI)

  elif mode == 'both':
    m.pcolor(XImap,YImap,ZI,cmap='jet')
    m.contour(XImap,YImap,ZI)

  m.scatter(xvals_map,yvals_map,120,zvals,cmap='jet',vmin=32,vmax=52)

  lat_draw = arange(latrange[0],latrange[1],1)
  lon_draw = arange(lonrange[0],lonrange[1],1)

  m.drawmeridians(lon_draw,labels=[0,0,0,1])
  m.drawparallels(lat_draw,labels=[1,0,0,0])

  if z == 'H':
    colorbar(shrink=0.7,label='Moho depth [km]')
  elif z == 'K':
    colorbar(shrink=0.7,label='vp/vs')

  if save:
    savefig('Mohomap_Vers1.pdf')

"""
def likelihood_filters()
"""




"""
def plot_RF_phases():

def get_arrivals():
"""

