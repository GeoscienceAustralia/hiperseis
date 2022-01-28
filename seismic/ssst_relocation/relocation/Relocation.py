# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 14:28:34 2021

@author: U37509
"""

import os, subprocess
import numpy as np
from datetime import datetime

def IsP(phase):
    """
    Determine if a wave arrives in P phase or S phase, looking at the last lag.
    
    
    Parameters
    ----------
    phase : string
        Phase
        
    
    Returns
    -------
    isp : integer
        1 if P wave, -1 if S wave, 0 if undetermined.
        
        
    """
    isp = 0
    if phase == 'p':
        isp = 1
    elif phase == 's':
        isp = -1
    else:
        l = len(phase)
        for i in np.arange(l, 0, -1):
            char = phase[i-1]
            if char == 'P': 
                isp = 1
                break
            elif char == 'S' or char == 'L':
                isp = -1
                break
            #end if
        #end for
    #end if
    return isp
#end func
    
def ang_dist(lon1, colat1, lon2, colat2, units='degrees'):
    """
    Function to calculate the angular distance from (lon1, colat1) to 
    (lon2, colat2).
    
    
    Parameters
    ----------
    lon1 : float
        Longitude of point 1.
        
    colat1 : float
        Colatitude of point 1.
        
    lon2 : float
        Longitude of point 2.
    
    colat2 : float
        Colatitude of point 2.
        
    units : string (optional)
        'degrees' or 'radians' describing the units in which lon1, lat1, lon2,
        lat2 are input. Default is 'degrees'.
        
        
    Returns
    -------
    value : float
        Angular distance from (lon1, colat1) to (lon2, colat2).
        
        
    """
    if units == 'degrees':
        lon1 = lon1*np.pi/180
        colat1 = colat1*np.pi/180
        lon2 = lon2*np.pi/180
        colat2 = colat2*np.pi/180
    else:
        colat1 = np.pi/2 - colat1
        colat2 = np.pi/2 - colat2
    #end if
    value = 2*np.arcsin(np.sqrt(np.sin((colat1 - colat2)/2)**2 + \
                                np.sin((lon1 - lon2)/2)**2* \
                                np.sin(colat1)*np.sin(colat2)))
    if units == 'degrees': value = value*180/np.pi
    return value
#end func

def push_time_corrections_to_database(picks):
    """
    Pushes time correction to SeisComp3 SQL database before relocation is 
    performed.
    
    
    Parameters
    ----------
    picks : numpy.ndarray
        Structured array with the following data type.
        dtype = [('event_id', 'U20'), ('pick_id', 'U30'), ('stat', 'U10'), 
                 ('net', 'U5'), ('cha', 'U10'), ('elon', 'single'), 
                 ('ecolat', 'single'), ('edepth', 'single'), 
                 ('origin_time', 'double'), ('mag', 'half'), 
                 ('slon', 'single'), ('scolat', 'single'), ('selev', 'half'), 
                 ('phase', 'U8'), ('arrival_time', 'double'), 
                 ('ptt', 'single'), ('tcor', 'half'), ('residual', 'half'), 
                 ('snr', 'half'), ('qualityMeasureCWT', 'half'), 
                 ('domFreq', 'half'), ('qualityMeasureSlope', 'half'), 
                 ('bandIndex', 'uint8'), ('nSigma', 'uint8')]
        
        
    """
    import MySQLdb
    
    print('Pushing time corrections to database')
    pick_ids = picks['pick_id']
    arrival_times = picks['arrival_time']
    tcors = picks['tcor']
    phases = picks['phase']
    
    t = arrival_times - tcors
    time_value = [datetime.fromtimestamp(int(ti)). \
                  strftime('%Y-%m-%d %H:%M:%S') for ti in t]
    time_value_ms = [int((ti - int(ti))*1e6) for ti in t]
    
    temp = [(pick_ids[i], time_value[i], time_value_ms[i], phases[i]) \
            for i in range(len(picks))]
    
    batch_size = 1000
    num_batches = int(np.ceil(len(temp)/batch_size))
    temp_split = list(np.array_split(temp, num_batches))        
    
    db = MySQLdb.connect(host="localhost", user="sysop", passwd="sysop", 
                         db="seiscomp3") 
    c = db.cursor()
    
    print('Removing old tables')
    c.execute('drop table if exists temp')
    db.commit()
    
    c.execute('drop table if exists picks_temp')
    db.commit()
    
    print('Adding new arrival times to database')
    c.execute(str('create table temp (pickID varchar(255), ' + \
                  'time_value datetime, time_value_ms int(11), ' + \
                  'phase char(32))'))
    db.commit()
    
    for lst in temp_split:
        query = str('insert into temp (pickID, time_value, time_value_ms, ' + \
                    'phase) values (%s, %s, %s, %s)')
        c.executemany(query, [tuple(row) for row in lst])
        db.commit()
    #end for
    
    print('Matching arrival times to picks')
    c.execute(str('create table picks_temp(index (pickID), index (_oid)) ' + \
                  'select * from temp join PublicObject on ' + \
                  'temp.pickID=PublicObject.publicID'))
    db.commit()
    
    print('Updating good picks')
    c.execute('alter table Pick order by _oid')
    db.commit()
    
    c.execute(str('update Pick right join picks_temp on ' + \
                  'Pick._oid=picks_temp._oid set Pick.time_value_ms=' + \
                  'picks_temp.time_value_ms, Pick.time_value=' + \
                  'picks_temp.time_value, Pick.phaseHint_code=' + \
                  'picks_temp.phase'))
    db.commit()
    
    c.close()
    db.close()
    
    return
#end func
    
def update_hypocentres_from_database(events, picks, hypo_dict, config, 
                                     unstable_events=list()):    
    """
    Attaches new hypocentres to their corresponding picks, or if the hypocentre 
    has changed by more than 'thr_dist' or 'thr_time', event ID is added to 
    'unstable_events'.
    
    
    Parameters
    ----------
    events : list
        List of event IDs for which to retrieve updated hypocentres.
        
    picks : numpy.ndarray
        Structured array with the following data type.
        dtype = [('event_id', 'U20'), ('pick_id', 'U30'), ('stat', 'U10'), 
                 ('net', 'U5'), ('cha', 'U10'), ('elon', 'single'), 
                 ('ecolat', 'single'), ('edepth', 'single'), 
                 ('origin_time', 'double'), ('mag', 'half'), 
                 ('slon', 'single'), ('scolat', 'single'), ('selev', 'half'), 
                 ('phase', 'U8'), ('arrival_time', 'double'), 
                 ('ptt', 'single'), ('tcor', 'half'), ('residual', 'half'), 
                 ('snr', 'half'), ('qualityMeasureCWT', 'half'), 
                 ('domFreq', 'half'), ('qualityMeasureSlope', 'half'), 
                 ('bandIndex', 'uint8'), ('nSigma', 'uint8')]
    
    hypo_dict : dictionary
        Dictionary of updated event hypocentres.
    
    config : configparser.SectionProxy object
        Information from config file.
    
    unstable_events : list
        List of event IDs which have a poor origin quality.
        
    
    Returns
    -------
    picks : numpy.ndarray
        Input 'picks' array with updated hypocentres.
        
    unstable_events : list
        Input 'unstable_events' list with added event IDs where hypocentre is
        unstable.
        
        
    """
    
    thr_dist = float(config['hypo_thr_dist_deg'])
    thr_time = float(config['hypo_thr_time_sec'])
    
    for event in events:
        inds = picks['event_id'] == event
        picks_temp = picks[inds]
        lon1 = picks_temp['elon'][0]
        lat1 = 90.0 - picks_temp['ecolat'][0]
        time1 = picks_temp['origin_time'][0]
        
        lon2, lat2, depth2, time2 = hypo_dict[event]
        
        if np.abs(lon1 - lon2) > thr_dist or np.abs(lat1 - lat2) > thr_dist \
            or np.abs(time1 - time2) > thr_time: 
            if event not in unstable_events:
                unstable_events.append(event)
            #end if
            continue
        #end if
        
        picks['ecolat'][inds] = 90.0 - lat2
        picks['elon'][inds] = lon2
        picks['edepth'][inds] = depth2
        picks['origin_time'][inds] = time2
    #end for
    return picks, unstable_events
#end func
    
def extract_hypocentres_from_database(events):
    """
    Retrieves updated hypocentres from SeisComp3 database after relocation has
    been performed.
    
    
    Parameters
    ----------
    events : list
        List of event IDs for which to retrieve updated hypocentres.
        
    
    Returns
    -------
    hypo_dict : dictionary
        Dictionary of updated event hypocentres.
        
        
    """
    import MySQLdb
    from obspy import UTCDateTime
    
    db = MySQLdb.connect(host="localhost", user="sysop", passwd="sysop", 
                         db="seiscomp3")
    c = db.cursor()
    c.execute(str('select ep.publicID, o.longitude_value, ' + \
                  'o.latitude_value, o.depth_value, o.time_value, ' + \
                  'o.time_value_ms, o._oid from Origin o, ' + \
                  'PublicObject op, Event e, PublicObject ep where ' + \
                  'o._oid=op._oid and e._oid=ep._oid and ' + \
                  'e.preferredOriginID=op.publicID'))
    rows = c.fetchall()
    hypo_dict = {row[0]: (float(row[1]), float(row[2]), float(row[3])*1e3,
                          UTCDateTime(row[4]).timestamp + int(row[5])/1e6) \
                          for row in rows if row[0] in events}
    c.close()
    db.close()
    
    return hypo_dict
#end func    
    
def compute_new_hypocentre(events, picks, TT_dict, ellipcorr_dict, output_path, 
                           rank, config, unstable_events=list(), fast=False):    
    """
    Uses a fortran subroutine to perform relocation of events using travel
    time corrections.
    
    
    Parameters
    ----------
    events : list
        List of IDs of events to be relocated.
        
    picks : numpy.ndarray
        Structured array with the following data type.
        dtype = [('event_id', 'U20'), ('pick_id', 'U30'), ('stat', 'U10'), 
                 ('net', 'U5'), ('cha', 'U10'), ('elon', 'single'), 
                 ('ecolat', 'single'), ('edepth', 'single'), 
                 ('origin_time', 'double'), ('mag', 'half'), 
                 ('slon', 'single'), ('scolat', 'single'), ('selev', 'half'), 
                 ('phase', 'U8'), ('arrival_time', 'double'), 
                 ('ptt', 'single'), ('tcor', 'half'), ('residual', 'half'), 
                 ('snr', 'half'), ('qualityMeasureCWT', 'half'), 
                 ('domFreq', 'half'), ('qualityMeasureSlope', 'half'), 
                 ('bandIndex', 'uint8'), ('nSigma', 'uint8')]
        
    TT_dict : dictionary
        Dictionary of phase_TT_table objects, containing travel time tables to
        use for relocation.
    
    ellipcorr_dict : dictionary
        Dictionary of ellipcorr_object objects, containing ellipticity 
        correction tables to use for relocation.
        
    output_path : string
        Output directory.
        
    rank : integer
        Rank of in-use processor.
        
    config : configparser.SectionProxy object
        Information from config file.
        
    unstable_events : list
        List of event IDs with a poor origin quality. 
        
    fast : boolean
        If True, only the first 100 picks for an event are used for relocation.
        
    
    Returns
    -------
    picks : numpy.ndarray
        Input pick array with updated hypocentres and travel time residuals.
        
    unstable_events : list
        Input unstable event list with added event IDs for events which move
        more than thr_dist degrees or change origin time by more than thr_time
        seconds after relocation.
        
        
    """
    
    from relocation_helper import relocation
    
    temp_networks = np.array(config['temp_networks'].split(', '))
    phases = np.array(config['phases'].split(', '))
    dcolat = float(config['reloc_dlat'])
    ddep = float(config['reloc_ddep'])
    niter = int(config['reloc_nit'])
    norm = int(config['reloc_norm'])
    frac = float(config['reloc_sfrac'])
    thr_dist = float(config['hypo_thr_dist_deg'])
    thr_time = float(config['hypo_thr_time_sec'])
    
    degrad = np.pi/180.0
    phase_alias = {'P': 'P', 'Pg': 'Pup', 'Pn': 'P', 'Pb': 'P', 'S': 'S', 
                   'Sg': 'Sup', 'Sn': 'S', 'Sb': 'S'}
    phase_ind = {phases[i]:i for i in range(len(phases))}
    
    tt_tables = np.zeros((500, 100, len(phases)))
    dtdd_tables = np.zeros((500, 100, len(phases)))
    ecdists_1 = np.zeros((500, len(phases)))
    depths_1 = np.zeros((100, len(phases)))
    
    for i in range(len(phases)):
        phase = phases[i]
        ecdists_temp = TT_dict[phase].ecdists
        depths_temp = TT_dict[phase].depths
        nx = len(ecdists_temp)
        nd = len(depths_temp)
        ecdists_1[:nx,i] = ecdists_temp
        depths_1[:nd,i] = depths_temp
        tt_tables[:nx,:nd,i] = TT_dict[phase].tt
        dtdd_tables[:nx,:nd,i] = TT_dict[phase].dtdd
    #end for
    
    tau0 = np.zeros((50, 6, len(phases)))
    tau1 = np.zeros((50, 6, len(phases)))
    tau2 = np.zeros((50, 6, len(phases)))
    ecdists_2 = np.zeros((50, len(phases)))
    depths_2 = np.zeros((6, len(phases)))
    
    for i in range(len(phases)):
        phase = phase_alias[phases[i]]
        ecdists_temp = ellipcorr_dict[phase].ecdist
        depths_temp = ellipcorr_dict[phase].depth
        nx = len(ecdists_temp)
        nd = len(depths_temp)
        ecdists_2[:nx,i] = ecdists_temp
        depths_2[:nd,i] = depths_temp
        tau0[:nx,:nd,i] = ellipcorr_dict[phase].tau0
        tau1[:nx,:nd,i] = ellipcorr_dict[phase].tau1
        tau2[:nx,:nd,i] = ellipcorr_dict[phase].tau2
    #end for
    
    filename = os.path.join(output_path, 'out%s.txt'%str(rank).zfill(3))
    with open(filename, 'w') as file:
        file.write(str('Relocating ' + str(len(events)) + ' events \n'))
    #end with
    
    good_pick_ind = np.logical_and(np.isin(picks['phase'], phases),
                                   ~np.isin(picks['net'], temp_networks))
    
    for event in events:
        with open(filename, 'a') as file:
            file.write(str('Relocating event ' + event + '\n'))
        #end with
        inds = np.logical_and(picks['event_id'] == event, good_pick_ind)
        picks_temp = picks[inds]
        
        if len(picks_temp) == 0: 
            with open(filename, 'a') as file:
                file.write(str('Finished relocating event ' + event + '\n'))
            #end with
            continue
        elif len(picks_temp) > 100 and fast == True:
            ecdist = ang_dist(picks_temp['elon'], picks_temp['ecolat'],
                              picks_temp['slon'], picks_temp['scolat'])
            picks_temp = picks_temp[ecdist.argsort()]
            npick_use = 100
        else:
            npick_use = len(picks_temp)
        #end if
        
        with open(filename, 'a') as file:
            file.write(str(str(len(picks_temp)) + ' picks \n'))
        #end with
        
        ecolat0 = picks_temp['ecolat'][0]
        elon0 = picks_temp['elon'][0]
        edep0 = picks_temp['edepth'][0]/1e3
        scolat0 = picks_temp['scolat']
        slon0 = picks_temp['slon']
        selev0 = picks_temp['selev']
        iph0 = np.array([phase_ind[pick['phase']] + 1 for pick in picks_temp])
        wt0 = np.array([IsP(pick['phase']) for pick in picks_temp])
        tt = picks_temp['arrival_time'] - picks_temp['origin_time']
        term = picks_temp['tcor']        
    
        phases_temp = np.full((len(picks_temp), 8), ' ', dtype='U1')
        for i in range(len(picks_temp)):
            for j in range(len(picks_temp['phase'][i])):
                phases_temp[i, j] = picks_temp['phase'][i][j]
            #end for
        #end for
        
        dlon = dcolat/np.sin(ecolat0*degrad)
        
        elon, ecolat, edep, ot, resid, qual = \
            relocation(elon0, ecolat0, edep0, slon0, scolat0, selev0, iph0, 
                       phases_temp, wt0, tt, term, dlon, dcolat, ddep, niter, 
                       frac, norm, tt_tables, dtdd_tables, ecdists_1, 
                       depths_1, ecdists_2, depths_2, tau0, tau1, tau2, 
                       npick_use)
        
        if qual == -1 or np.abs(elon - elon0) > thr_dist \
            or np.abs(ecolat - ecolat0) > thr_dist or np.abs(ot) > thr_time:
            if event not in unstable_events:
                unstable_events.append(event)
            #end if
            with open(filename, 'a') as file:
                file.write(str('Finished relocating event ' + event + '\n'))
            #end with
            continue
        #end if
        
        #picks['residual'][inds] = resid
        
        inds = picks['event_id'] == event
        
        picks['ecolat'][inds] = ecolat
        picks['elon'][inds] = elon
        picks['edepth'][inds] = edep*1e3
        picks['origin_time'][inds] = picks['origin_time'][inds] + ot
        
        with open(filename, 'a') as file:
            file.write(str('Finished relocating event ' + event + '\n'))
        #end with
    #end for
    with open(filename, 'a') as file:
        file.write(str('Finished relocating events \n'))
    #end with
    return picks, unstable_events
#end func
    
def compute_new_hypocentre_iloc(events, output_path, config, rank):
    """
    Passes the time corrections included in 'picks' into a SeisComp3 database 
    for iLoc to use, then executes commands in the command line to perform 
    relocation with iLoc.
    
    
    Parameters
    ----------
    events : list
        List of event IDs for which to perform relocation.
        
    output_path : string
        Output directory.
        
    config : configparser.SectionProxy object
        Information from config file.
        
    rank : integer
        Rank of in-use processor.
        
        
    """
    from tqdm import tqdm
    
    iloc_redefine_phases = config['iloc_redefine_phases'] == 'True'
    iloc_use_rstt = config['iloc_use_rstt'] == 'True'
    
    filename = os.path.join(output_path, 'out%s.txt'%str(rank).zfill(3))
    with open(filename, 'w') as file:
        file.write(str('Relocating ' + str(len(events)) + ' events \n'))
    #end with
    
    if iloc_redefine_phases:
        DoNotRenamePhase = '0'
    else:
        DoNotRenamePhase = '1'
    #end if
    if iloc_use_rstt:
        UseRSTT = '1'
    else:
        UseRSTT = '0'
    #end if
    
    devnull = open(os.devnull, 'wb')
    for event in tqdm(events):
        with open(filename, 'a') as file:
            file.write(str('Relocating event ' + event + '\n'))
        #end with
        
        cmd = 'echo "' + str(event) + ' UpdateDB=1 DoGridSearch=1 ' + \
              'DoNotRenamePhase=' + DoNotRenamePhase + \
              ' UseRSTTPnSn=' + UseRSTT + ' UseRSTTPgLg=' + UseRSTT + \
              ' Verbose=0 NAsearchRadius=1" | iLocSC seiscomp'
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, 
                             stderr=devnull)
        results = list()
        for line in p.stdout:
            results.append(line.strip())
        #end for
        p.wait()
        if p.returncode:
            print('Event', event, 'encountered an error during relocation!')
            with open(os.path.join(event.replace('/', '-'), 
                                   '.txt'), 'a') as file:
                for line in results:
                    file.write(str(line + '\n'))
                #end for
            #end with
        #end if
        with open(filename, 'a') as file:
            file.write(str('Finished relocating event ' + event + '\n'))
        #end with
    #end for
    with open(filename, 'a') as file:
        file.write(str('Finished relocating events \n'))
    #end with
    return
#end func