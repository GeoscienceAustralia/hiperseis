"""
Description
-----------
This script is used to read event catalogues in either .pkl binary format as 
obspy.core.event Catalog objects, or QuakeML/SC3ML .xml format, then merge and 
filter the catalogues.


Criteria
--------
- Events within the Australian bounding box [105, -45, 160, -13] are filtered 
to have azimuthal gap less than 180 degrees.

- Events in the Indian Ocean in bounding boxes [15, -60, 90, 30] and [90, -60, 
180, 45] are filtered to have a magnitude larger than 5.

- Events in the rest of the world are filtered to have magnitude larger than 
5.5.

- Events with +/- 15 seconds and within 50km of one another are counted as 
duplicates, and only one is kept.

- Hypocentres sourced from Gary Gibson are used to replace those in the GA 
catalogues if a duplicate is found.


Developer: Lachlan Adams 
Contact: lachlan.adams@ga.gov.au or lachlan.adams.1996@outlook.com

"""

from obspy.core.event import Catalog, Event, Origin, Pick, WaveformStreamID, \
                             Arrival, Magnitude, CreationInfo, QuantityError, \
                             OriginQuality, OriginUncertainty, read_events
from obspy import UTCDateTime
import glob, os, time, pickle, argparse
import pandas as pd
import numpy as np

def azimuth(lon1, lat1, lon2, lat2, units='degrees'):
    """
    Function to calculate the azimuth from (lon1, lat1) to (lon2, lat2).
    
    
    Parameters
    ----------
    lon1 : float
        Longitude of point 1.
        
    lat1 : float
        Latitude of point 1.
        
    lon2 : float
        Longitude of point 2.
    
    lat2 : float
        Latitude of point 2.
        
    units : string (optional)
        'degrees' or 'radians' describing the units in which lon1, lat1, lon2,
        lat2 are input. Default is 'degrees'.
        
        
    Returns
    -------
    azim : float
        Azimuth from (lon1, colat1) to (lon2, colat2).
        
        
    """
    if units == 'degrees':
        degrad = np.pi/180.0
        a = lon1*degrad
        b = (90.0-lat1)*degrad
        x = lon2*degrad
        y = (90.0-lat2)*degrad
    else:
        a = lon1
        b = np.pi/2.0 - lat1
        x = lon2
        y = mp.pi/2.0 - lat2
    #end if
    azim = np.arctan(np.sin(x - a)/(np.sin(b)/np.tan(y) - \
                                    np.cos(b)*np.cos(x - a)))
    
    offset1 = np.pi*(lat1 > lat2)
    offset2 = np.pi*(np.logical_and(lat1 == lat2, b > np.pi/2.0))
    
    azim = (azim + offset1 + offset2 + 2*np.pi) % (2*np.pi)
    if units == 'degrees': 
        azim = azim/degrad
    return azim
#end func
    
def ang_dist(lon1, lat1, lon2, lat2, units='degrees'):
    """
    Function to calculate the angular distance from (lon1, lat1) to 
    (lon2, lat2).
    
    
    Parameters
    ----------
    lon1 : float
        Longitude of point 1.
        
    lat1 : float
        Latitude of point 1.
        
    lon2 : float
        Longitude of point 2.
    
    lat2 : float
        Latitude of point 2.
        
    units : string (optional)
        'degrees' or 'radians' describing the units in which lon1, lat1, lon2,
        lat2 are input. Default is 'degrees'.
        
        
    Returns
    -------
    value : float
        Angular distance from (lon1, lat1) to (lon2, lat2).
        
        
    """
    if units == 'degrees':
        colat1 = 90 - lat1
        colat2 = 90 - lat2
        lon1 = lon1*np.pi/180
        colat1 = colat1*np.pi/180
        lon2 = lon2*np.pi/180
        colat2 = colat2*np.pi/180
    else:
        colat1 = np.pi/2 - lat1
        colat2 = np.pi/2 - lat2
    #end if
    value = 2*np.arcsin(np.sqrt(np.sin((colat1 - colat2)/2)**2 + \
                                np.sin((lon1 - lon2)/2)**2* \
                                np.sin(colat1)*np.sin(colat2)))
    if units == 'degrees': value = value*180/np.pi
    return value
#end func
    
def read_preprocessed(path):
    """
    Function to read event information already ingested from xml files. 
    
    
    Parameters
    ----------
    path : string
        Directory in which all files are stored. The function will also look in
        all sub-directories of 'path'.
        
        
    Returns
    -------
    catalogue : obspy.core.event Catalog object
        Catalogue of all events read in from files in directory 'path' and its
        sub-directories.
        
    
    """
    catalogue = Catalog(events = None)
    
    files = glob.glob(os.path.join(path, '**/*.pkl'), recursive=True)
    
    if len(files) != 0:
        print('Number of files =', len(files))
        with open('out.txt', 'a') as temp:
            temp.write(str('Number of files = ' + str(len(files)) + '\n'))
        i = 0
        t0 = time.time()
        for file in files:
            i = i + 1
            #print('Reading file', i, 'of', len(files), ', time =', 
            #      time.time()-t0)
            with open('out.txt', 'a') as temp:
                temp.write(str('Reading file ' + str(i) + ' of ' + \
                               str(len(files)) + ', time = ' + \
                               str(time.time()-t0) + '\n'))
            with open(file, 'rb') as inp:
                newcat = pickle.load(inp)
            #end with
            for event in newcat:
                if event.creation_info is None:
                    event.creation_info = CreationInfo(agency_id = 'X', 
                                                       auhtor = 'X')
                #end if
            #end for
            catalogue.extend(list(newcat))
        #end for
    #end if
    
    return catalogue
#end func
    
def read_event_xml_files(path):
    """
    Function to read event information from xml files. First tries to read in
    SC3ML format, and if this fails, attempts to read in QuakeML format.
    
    
    Parameters
    ----------
    path : string
        Directory in which all files are stored. The function will also look in
        all sub-directories of 'path'.
        
        
    Returns
    -------
    catalogue : obspy.core.event Catalog object
        Catalogue of all events read in from files in directory 'path' and its
        sub-directories.
        
    
    """
    catalogue = Catalog(events = None)
    
    if path == "":
        files = list()
    else:
        files = glob.glob(os.path.join(path, '**/*.xml'), recursive=True)
    #end if
    
    if len(files) != 0:
        for file in files:
            try:
                newcat = read_events(file, format = 'SC3ML')
            except:
                try:
                    newcat = read_events(file, format = 'QuakeML')
                except:
                    continue
                #end try
            #end try
            for event in newcat:
                if event.preferred_origin() is None:
                    event.preferred_origin_id = event.origins[0].resource_id
                #end if
                if event.preferred_magnitude() is None:
                    station_count_list = list()
                    for magnitude in event.magnitudes:
                        station_count_list.append(magnitude.station_count)
                    #end for
                    max_count = max(station_count_list)
                    max_ind = station_count_list.index(max_count)
                    event.preferred_magnitude_id = \
                        event.magnitudes[max_ind].resource_id
                #end if
                if event.creation_info is None:
                    event.creation_info = CreationInfo(agency_id = 'X', 
                                                       auhtor = 'X')
                #end if
            #end for
            catalogue.extend(newcat)
        #end for
    #end if
    
    return catalogue
#end func
    
def read_GG_csv_files(path):
    """
    A function to read event hypocentre information from .csv files in the 
    format used by Gary Gibson.
    
    
    Parameters
    ----------
    path : string
        Directory in which all files are stored.
        
        
    Returns
    -------
    catalogue : obspy.core.event Catalog object
        Catalogue of events contained in input .csv file.
        
        
    """
    GG_columns = ['RunAuth', 'Type', 'Ind', 'Ass', 'Year', 'Month', 'Day', 
                  'Hour', 'Minute', 'Second', 'STCode', 'STCorr', 'LongitudeE',
                  'LatitudeS', 'Depth', 'ZCode', 'Mx', 'Mval', 'TimUnc', 
                  'Smaj', 'Smin', 'SmajAz', 'DepUnc', 'Arrivals', 'SDResids',
                  'Sites', 'Gap', 'Gap2', 'NearestSG', 'Nearest', 'HDPlace',
                  'Origin (text)', 'YearF', 'Magnitude text']
    GG_dtype = [str, str, int, str, int, int, int, int, int, float, str, float, 
                float, float, float, str, str, float, float, float, float, 
                float, float, int, float, int, float, float, str, float, str, str, 
                str, str]
    
    if path == "":
        files = list()
    else:
        files = glob.glob(os.path.join(path, '*.csv'))
    #end if
    
    if len(files) == 0:
        return Catalog(events=None)
    else:
        df = pd.DataFrame(columns=GG_columns)
        
        for file in files:
            df_temp = pd.read_csv(file)
            df_temp.columns = GG_columns
            df = pd.concat([df, df_temp], axis=0)
        #end for
        
        df = df.replace(r'^\s*$', np.nan, regex=True)
        df.fillna('0', inplace=True)    
        
        df = df.astype({GG_columns[i]: GG_dtype[i] \
                        for i in range(len(GG_columns))})
        
        df['n_phase'] = [0 for i in range(len(df))]
        df['longitude'] = df['LongitudeE']
        df['latitude'] = -df['LatitudeS']
        df['event_id'] = [str('#GG' + str(i + 1)) for i in range(len(df))]
        
        df = df.drop(df[df['Second'] >= 60].index)
        df = df.drop(df[df['Minute'] >= 60].index)
        df = df.drop(df[df['Hour'] >= 24].index)
        df.reset_index(inplace=True, drop=True)
    #end if
    lst = list()
    for i in range(len(df)):
        data = df.loc[i]
        event_id = data['event_id']
        YYYY = data['Year']
        MM = data['Month']
        DD = data['Day']
        hh = data['Hour']
        mm = data['Minute']
        ss = data['Second']
        lon = data['longitude']
        lat = data['latitude']
        depth = data['Depth']
        n_phase = data['n_phase']
        azim_gap = data['Gap']
        if azim_gap == 0:
            azim_gap = 999.0
        #end if
        mag_type = data['Mx']
        mag = data['Mval']
        mb = ms = ml = mw = -999.0
        if mag_type == 'mb' or mag_type == 'MB' or mag_type == 'Mb' or \
            mag_type == 'mB':
            mb = mag
        elif mag_type == 'ms' or mag_type == 'MS' or mag_type == 'Ms' or \
            mag_type == 'mS':
            ms = mag
        elif mag_type == 'ml' or mag_type == 'ML' or mag_type == 'Ml' or \
            mag_type == 'mL':
            ml = mag
        elif mag_type == 'mw' or mag_type == 'MW' or mag_type == 'Mw' or \
            mag_type == 'mW':
            mw = mag
        #end if
        lst.append([event_id, YYYY, MM, DD, hh, mm, ss, lon, lat, depth, 
                    n_phase, mb, ms, ml, mw, azim_gap])
    #end for
    return list_to_catalogue(lst)
#end func
    
def list_to_catalogue(lst):
    """
    Converts a list of events into an obspy.core.event Catalog object. List 
    format:
        event_rows = [event_id, YYYY, MM, DD, hh, mm, ss, lon, lat, depth, 
                      n_phase, mb, ms, ml, mw, azim_gap]
        pick_rows = [event_id, sta, cha, loc, net, lon, lat, elev, phase, 
                     YYYY, MM, DD, hh, mm, ss, ang_dist]
        
        
    Parameters
    ----------
    lst : list
        List of event information to be converted.
        
    
    Returns
    -------
    catalogue : obspy.core.event Catalog object
        Catalogue of all events in lst.
        
        
    """
    
    catalogue = Catalog(events=None)
    
    event_list = [row for row in lst if row[0][0] == '#']
    
    for event in event_list:
        event_id, YYYY, MM, DD, hh, mm, ss, longitude, latitude, depth, \
            n_phase, mb, ms, ml, mw, azim_gap = event
        
        if np.isnan(azim_gap):
            azim_gap = 999.0
        #end if
        event_id = event_id[1:]
        origin_time = UTCDateTime(year=YYYY, month=MM, day=DD, hour=hh,
                                  minute=mm, second=int(ss), 
                                  microsecond=int(1e6*(ss - int(ss))))
        pick_list = [row for row in lst if row[0] == event_list[0][1:]]
        picks = list()
        arrivals = list()
        
        i = 0
        for pick in pick_list:
            i = i + 1
            event_id, sta, cha, loc, net, lon, lat, elev, phase, YYYY, MM, \
                DD, hh, mm, ss, dist = pick
            
            if lon is None or lat is None:
                azim = 0
                ecdist = 0
            else:
                azim = azimuth(longitude, latitude, lon, lat)
                backazim = azimuth(lon, lat, longitude, latitude)
                ecdist = ang_dist(longitude, latitude, lon, lat)
            #end if
            
            pick_obj = Pick(resource_id = str(event_id + '_p_' + str(i)),
                            time = UTCDateTime(year=YYYY, month=MM, day=DD, \
                                   hour=hh, minute=mm, second=int(ss), \
                                   microsecond=int(1e6*(ss - int(ss)))),
                            backazimuth = backazim,
                            creation_info = CreationInfo(agency_id = 'X',
                                                         author = 'X'),
                            waveform_id = WaveformStreamID(network_code=net,
                                                           station_code=sta,
                                                           channel_code=cha,
                                                           location_code=loc),
                            time_errors = QuantityError(),
                            filter_id = None,
                            method_id = None,
                            horizontal_slowness = None,
                            horizontal_slowness_errors = QuantityError(),
                            backazimuth_errors = QuantityError(),
                            slowness_method_id = None,
                            onset = None,
                            polarity = None,
                            evaluation_mode = None,
                            evaluation_status = None,
                            comments = list())
            picks.append(pick_obj)
            
            arrival_obj = Arrival(resource_id = str(event_id + '_a_' + str(i)),
                                  pick_id = str(event_id + '_p_' + str(i)),
                                  phase = phase,
                                  azimuth = azim,
                                  distance = ecdist,
                                  creation_info = CreationInfo(agency_id = 'X',
                                                               author = 'X'),
                                  time_correction = None,
                                  takeoff_angle = None,
                                  takeoff_angle_errors = None,
                                  time_residual = None,
                                  horizontal_slowness_residual = None,
                                  backazimuth_residual = None,
                                  time_weight = None,
                                  horizontal_slowness_weight = None,
                                  backazimuth_weight = None,
                                  earth_model_id = None,
                                  comments = list())
            arrivals.append(arrival_obj)
        #end for
        
        mb_obj = Magnitude(resource_id = str(event_id + '_mb'),
                           mag = mb,
                           mag_errors = QuantityError(),
                           magnitude_type = 'mb',
                           origin_id = str(event_id + '_origin'),
                           creation_info = CreationInfo(agency_id = 'X',
                                                        author = 'X'))
        ms_obj = Magnitude(resource_id = str(event_id + '_ms'),
                           mag = ms,
                           mag_errors = QuantityError(),
                           magnitude_type = 'ms',
                           origin_id = str(event_id + '_origin'),
                           creation_info = CreationInfo(agency_id = 'X',
                                                        author = 'X'))
        ml_obj = Magnitude(resource_id = str(event_id + '_ml'),
                           mag = mb,
                           mag_errors = QuantityError(),
                           magnitude_type = 'ml',
                           origin_id = str(event_id + '_origin'),
                           creation_info = CreationInfo(agency_id = 'X',
                                                        author = 'X'))
        mw_obj = Magnitude(resource_id = str(event_id + '_mw'),
                           mag = mw,
                           mag_errors = QuantityError(),
                           magnitude_type = 'mw',
                           origin_id = str(event_id + '_origin'),
                           creation_info = CreationInfo(agency_id = 'X',
                                                        author = 'X'))
        pref_mag = ml_obj.resource_id
        if mw > 0:
            pref_mag = mw_obj.resource_id
        if ms > 0:
            pref_mag = ms_obj.resource_id
        if ml > 0:
            pref_mag = ml_obj.resource_id
        if mb > 0:
            pref_mag = mb_obj.resource_id
        #end if
        
        origin_object = Origin(resource_id = str(event_id + '_origin'),
                               time = origin_time,
                               longitude = longitude,
                               latitude = latitude,
                               depth = depth,
                               creation_info = CreationInfo(agency_id = 'X',
                                                            author = 'X'),
                               arrivals = arrivals,
                               time_errors = QuantityError(),
                               longitude_errors = QuantityError(),
                               latitude_errors = QuantityError(),
                               depth_errors = QuantityError(),
                               depth_type = None,
                               time_fixed = None,
                               epicentre_fixed = None,
                               reference_system_id = None,
                               method_id = None,
                               earth_model_id = None,
                               quality = OriginQuality(used_phase_count=None,
                                                       standard_error=None,
                                                       minimum_distance=None,
                                                       maximum_distance=None,
                                                       azimuthal_gap=azim_gap),
                               origin_type = None,
                               origin_uncertainty = \
                                   OriginUncertainty( \
                                       min_horizontal_uncertainty=None,
                                       max_horizontal_uncertainty=None,
                                       azimuth_max_horizontal_uncertainty=None,
                                       preferred_description=None),
                               region = None,
                               evaluation_mode = None,
                               evaluation_status = None,
                               comments = list(),
                               composite_times = list())
                                   
        event_object = Event(resource_id = event_id,
                        	 preferred_origin_id = str(event_id + '_origin'),
                        	 picks = picks,
                        	 amplitudes = list(),
                        	 focal_mechanisms = list(),
                        	 origins = [origin_object],
                        	 magnitudes = [mb_obj, ml_obj, ms_obj, mw_obj],
                        	 station_magnitudes = list(),
                        	 creation_info = CreationInfo(agency_id='X', 
                                                          author='X'),
                        	 event_type = None,
                        	 event_type_certainty = 'known',
                        	 preferred_magnitude_id = pref_mag,
                        	 comments = list(),
                        	 event_descriptions = list())
                             
        catalogue.append(event_object)
    #end for
        
    return catalogue
#end func
    
def merge_catalogues(catalogues):
    """
    Function to merge a list of obspy.core.event Catalog objects into a single 
    catalogue.
    
    
    Parameters
    ----------
    catalogues : list of obspy.core.event Catalog objects
        Catalogues to be merged.
        
        
    Returns
    -------
    merged_catalogue : obspy.core.event Catalog object
        Merged catalogue.
        
        
    """
    merged_catalogue = Catalog(events=None)
    for cat in catalogues:
        merged_catalogue.extend(cat)
    #end for
    return merged_catalogue
#end func
    
def write_catalogue_to_file(catalogue, path):
    """
    Function to save a catalogue of events as .pkl files. 
    
    
    Parameters
    ----------
    catalogue : obspy.core.event Catalog object
        Catalogue of events to write to file.
        
    path : string
        Directory in which to save files.
    
    
    """
    t0 = time.time()
    file_ext = '.pkl'
    l = len(catalogue)
    for i in range(len(catalogue)):
        with open('out.txt', 'a') as temp:
            temp.write(str('Saving event ' + str(i+1) + ' of ' + str(l) + \
                           ', time = ' + str(time.time() - t0) + '\n'))
        #end with
        cat = catalogue[i:i+1]
        try:
            ot = str(str(cat[0].origins[0].time).split('.')[0] + 'Z') \
                .replace(':', '-')
        except:
            continue
        #end try
        savefile = os.path.join(path, ot, file_ext)
        with open(savefile, 'wb') as outp:
            pickle.dump(cat, outp, pickle.HIGHEST_PROTOCOL)
        #end with
    #end for
#end func
    
def filter_ISC(arr, params):
    """
    Filter ISC catalogue by the following criteria:
        1) Events in aus_bbox have azimuthal gap less than aus_azim_gap_thr.
        2) Events in io_bbox_1 and io_bbox_2 have magnitude larger than 
           io_mag_thr.
        3) Events elsewhere have magnitude larger than mag_thr.
                
    
    Parameters
    ----------
    arr : numpy.ndarray
        Structured array containing fields event_id, lon, lat, azim_gap, mag
        
    params : dictionary
        Dictionary containing parameters for filtering
        
    
    Returns
    arr : numpy.ndarray
        Structured array of filtered events
        
        
    """
    mag_thr = params['mag_thr']
    io_mag_thr = params['io_mag_thr']
    aus_azim_gap_thr = params['aus_azim_gap_thr']
    aus_bbox = params['aus_bbox']
    io_bbox_1 = params['io_bbox_1']
    io_bbox_2 = params['io_bbox_2']
    
    #Filter by bbox restrictions
    in_aus_bbox = np.all([arr['lon']>aus_bbox[0], arr['lat']>aus_bbox[1],
                          arr['lon']<aus_bbox[2], arr['lat']<aus_bbox[3]], 
                         axis=0)
    in_io_bbox_1 = np.all([arr['lon']>io_bbox_1[0], arr['lat']>io_bbox_1[1],
                           arr['lon']<io_bbox_1[2], arr['lat']<io_bbox_1[3]], 
                          axis=0)
    in_io_bbox_2 = np.all([arr['lon']>io_bbox_2[0], arr['lat']>io_bbox_2[1],
                           arr['lon']<io_bbox_2[2], arr['lat']<io_bbox_2[3]], 
                          axis=0)
    in_other = np.all([~in_aus_bbox, ~in_io_bbox_1, ~in_io_bbox_2], axis=0)
    
    filt = list()
    filt.append(np.all([arr['azim_gap']<aus_azim_gap_thr, in_aus_bbox], 
                       axis=0))
    filt.append(np.all([arr['mag']>io_mag_thr, in_io_bbox_1], axis=0))
    filt.append(np.all([arr['mag']>io_mag_thr, in_io_bbox_2], axis=0))
    filt.append(np.all([arr['mag']>mag_thr, in_other], axis=0))
    
    filt = np.any(filt, axis=0)
    
    return arr[filt]
#end func
    
def filter_USGS(arr, params):
    """
    Filter USGS catalogue by the following criteria:
        1) Events in aus_bbox have azimuthal gap less than aus_azim_gap_thr.
        2) Events in io_bbox_1 and io_bbox_2 have magnitude larger than 
           io_mag_thr.
        3) Events elsewhere have magnitude larger than mag_thr.
                
    
    Parameters
    ----------
    arr : numpy.ndarray
        Structured array containing fields event_id, lon, lat, azim_gap, mag
        
    params : dictionary
        Dictionary containing parameters for filtering
        
    
    Returns
    arr : numpy.ndarray
        Structured array of filtered events
        
        
    """
    mag_thr = params['mag_thr']
    io_mag_thr = params['io_mag_thr']
    aus_azim_gap_thr = params['aus_azim_gap_thr']
    aus_bbox = params['aus_bbox']
    io_bbox_1 = params['io_bbox_1']
    io_bbox_2 = params['io_bbox_2']
    
    #Filter by bbox restrictions
    in_aus_bbox = np.all([arr['lon']>aus_bbox[0], arr['lat']>aus_bbox[1],
                          arr['lon']<aus_bbox[2], arr['lat']<aus_bbox[3]], 
                         axis=0)
    in_io_bbox_1 = np.all([arr['lon']>io_bbox_1[0], arr['lat']>io_bbox_1[1],
                           arr['lon']<io_bbox_1[2], arr['lat']<io_bbox_1[3]], 
                          axis=0)
    in_io_bbox_2 = np.all([arr['lon']>io_bbox_2[0], arr['lat']>io_bbox_2[1],
                           arr['lon']<io_bbox_2[2], arr['lat']<io_bbox_2[3]], 
                          axis=0)
    in_other = np.all([~in_aus_bbox, ~in_io_bbox_1, ~in_io_bbox_2], axis=0)
    
    filt = list()
    filt.append(np.all([arr['azim_gap']<aus_azim_gap_thr, in_aus_bbox], 
                       axis=0))
    filt.append(np.all([arr['mag']>io_mag_thr, in_io_bbox_1], axis=0))
    filt.append(np.all([arr['mag']>io_mag_thr, in_io_bbox_2], axis=0))
    filt.append(np.all([arr['mag']>mag_thr, in_other], axis=0))
    
    filt = np.any(filt, axis=0)
    
    return arr[filt]
#end func
    
def filter_GA(arr, params):
    """
    Filter GA catalogue by the following criteria:
        1) Events have azimuthal gap less than aus_azim_gap_thr.
                
    
    Parameters
    ----------
    arr : numpy.ndarray
        Structured array containing fields event_id, lon, lat, azim_gap, mag
        
    params : dictionary
        Dictionary containing parameters for filtering
        
    
    Returns
    arr : numpy.ndarray
        Structured array of filtered events
        
        
    """
    azim_gap_thr = params['aus_azim_gap_thr']
    filt = list()
    filt.append(arr['azim_gap'] < azim_gap_thr)
    filt = np.any(filt, axis=0)
    
    return arr[filt]
#end func
    
def filter_EHB(arr, params):
    """
    Filter Engdahl catalogue by the following criteria:
        None
                
    
    Parameters
    ----------
    arr : numpy.ndarray
        Structured array containing fields event_id, lon, lat, azim_gap, mag
        
    params : dictionary
        Dictionary containing parameters for filtering
        
    
    Returns
    arr : numpy.ndarray
        Structured array of filtered events
        
        
    """
    return arr
#end func
    
def filter_other(arr, params):
    """
    Filter other/unspecified catalogue by the following criteria:
        1) Events in aus_bbox have azimuthal gap less than aus_azim_gap_thr.
        2) Events in io_bbox_1 and io_bbox_2 have magnitude larger than 
           io_mag_thr.
        3) Events elsewhere have magnitude larger than mag_thr.
                
    
    Parameters
    ----------
    arr : numpy.ndarray
        Structured array containing fields event_id, lon, lat, azim_gap, mag
        
    params : dictionary
        Dictionary containing parameters for filtering
        
    
    Returns
    arr : numpy.ndarray
        Structured array of filtered events
        
        
    """
    mag_thr = params['mag_thr']
    io_mag_thr = params['io_mag_thr']
    aus_azim_gap_thr = params['aus_azim_gap_thr']
    aus_bbox = params['aus_bbox']
    io_bbox_1 = params['io_bbox_1']
    io_bbox_2 = params['io_bbox_2']
    
    #Filter by bbox restrictions
    in_aus_bbox = np.all([arr['lon']>aus_bbox[0], arr['lat']>aus_bbox[1],
                          arr['lon']<aus_bbox[2], arr['lat']<aus_bbox[3]], 
                         axis=0)
    in_io_bbox_1 = np.all([arr['lon']>io_bbox_1[0], arr['lat']>io_bbox_1[1],
                           arr['lon']<io_bbox_1[2], arr['lat']<io_bbox_1[3]], 
                          axis=0)
    in_io_bbox_2 = np.all([arr['lon']>io_bbox_2[0], arr['lat']>io_bbox_2[1],
                           arr['lon']<io_bbox_2[2], arr['lat']<io_bbox_2[3]], 
                          axis=0)
    in_other = np.all([~in_aus_bbox, ~in_io_bbox_1, ~in_io_bbox_2], axis=0)
    
    filt = list()
    filt.append(np.all([arr['azim_gap']<aus_azim_gap_thr, in_aus_bbox], 
                       axis=0))
    filt.append(np.all([arr['mag']>io_mag_thr, in_io_bbox_1], axis=0))
    filt.append(np.all([arr['mag']>io_mag_thr, in_io_bbox_2], axis=0))
    filt.append(np.all([arr['mag']>mag_thr, in_other], axis=0))
    
    filt = np.any(filt, axis=0)
    
    return arr[filt]
#end func
    
def replace_GA_with_GG(GA, GG, event_dict, thr_time, thr_dist):
    """
    Replace GA hypocentres with GG hypocentres if a match is found. Criteria:
        1) Origin time difference less than thr_time seconds
        2) Epicentral distance less than thr_dist degrees
        
        
    Parameters
    ----------
    GA : numpy.ndarray
        Structured array containing fields event_id, lon, lat, time, source
    
    GG : numpy.ndarray
        Structured array containing fields event_id, lon, lat, time, source
    
    event_dict : dictionary
        Dictionary of obspy.core.event Event objects
    
    thr_time : float
        Threshold time (seconds)
    
    thr_dist : float
        Threshold distance (degrees)
        
        
    Returns
    -------
    GA : numpy.ndarray
        Structured array of events with replaced hypocentres.
    """
    timediff_GG = np.abs(GA['time'][:, np.newaxis].T - \
                         GG['time'][:, np.newaxis])
    angdist_GG = ang_dist(GA['lon'][:, np.newaxis].T, 
                          GA['lat'][:, np.newaxis].T,
                          GG['lon'][:, np.newaxis],
                          GG['lat'][:, np.newaxis])
    match_GG = np.logical_and(timediff_GG < thr_time, angdist_GG < thr_dist)
    
    matches = np.array(np.where(match_GG)).T
    with open('out.txt', 'a') as temp:
        temp.write(str(str(len(matches)) + '\n'))
    np.save('GG_matches.npy', matches)
    for i, j in matches:
        GA['lon'][j] = GG['lon'][i]
        GA['lat'][j] = GG['lat'][i]
        GA['time'][j] = GG['time'][i]
        event_id_GA = GA['event_id'][j]
        event_id_GG = GG['event_id'][i]
        GA['source'][j] = 'GG'
        try:
            event_dict[event_id_GA].preferred_origin().longitude = \
                event_dict[event_id_GG].preferred_origin().longitude
            event_dict[event_id_GA].preferred_origin().latitude = \
                event_dict[event_id_GG].preferred_origin().latitude
            event_dict[event_id_GA].preferred_origin().depth = \
                event_dict[event_id_GG].preferred_origin().depth
            event_dict[event_id_GA].preferred_origin().time = \
                event_dict[event_id_GG].preferred_origin().time
        except:
            event_dict[event_id_GA].origins[0].longitude = \
                event_dict[event_id_GG].origins[0].longitude
            event_dict[event_id_GA].origins[0].latitude = \
                event_dict[event_id_GG].origins[0].latitude
            event_dict[event_id_GA].origins[0].depth = \
                event_dict[event_id_GG].origins[0].depth
            event_dict[event_id_GA].origins[0].time = \
                event_dict[event_id_GG].origins[0].time
        #end try
    #end for
    with open('out.txt', 'a') as temp:
        temp.write(str(str(np.sum(GA['source'] == 'GG')) + '\n'))
    return GA, event_dict
#end func
    
def remove_duplicates(event_array, event_dict, thr_time, thr_dist):
    """
    Sort events by origin time and remove duplicates if a match is found. 
    Criteria:
        1) Origin time difference less than thr_time seconds
        2) Epicentral distance less than thr_dist degrees
    If duplicates are found, preference is to keep Engdahl > ISC > USGS > GA >
    other, or if a tie is found, keep earliest event.
        
        
    Parameters
    ----------
    event_array : numpy.ndarray
        Structured array containing fields event_id, lon, lat, time, source
    
    event_dict : dictionary
        Dictionary of obspy.core.event Event objects
    
    thr_time : float
        Threshold time (seconds)
    
    thr_dist : float
        Threshold distance (degrees)
        
        
    Returns
    -------
    event_array : numpy.ndarray
        Structured array of events with replaced hypocentres.
        
    event_dict : dictionary
        Dictionary of obspy.core.event Event objects with updated origins and
        picks.
        
        
    """
    priority = {'EHB': 1, 'ISC': 2, 'USGS':3, 'GA': 4, 'GG': 4, 'other': 5}
    
    event_array['timediff'][1:] = \
        event_array['time'][1:] - event_array['time'][:-1]
    event_array['angdist'][1:] = \
        ang_dist(event_array['lon'][1:], event_array['lat'][1:],
                 event_array['lon'][:-1], event_array['lat'][:-1])
        
    small_timediff_rows = (event_array['timediff'] < thr_time)
    small_ang_dist_rows = (event_array['angdist'] < thr_dist)
    duplicate_rows = np.logical_and(small_timediff_rows, small_ang_dist_rows)
    duplicate_ind = np.where(duplicate_rows)[0]
    
    drop_rows = np.zeros_like(duplicate_rows).astype(bool)
    
    for ind in duplicate_ind:
        i0 = ind - 1
        i1 = ind
        event_id0 = event_array['event_id'][i0]
        event_id1 = event_array['event_id'][i1]
        source0 = event_array['source'][i0]
        source1 = event_array['source'][i1]
        
        if priority[source0] > priority[source1]:
            #Keep 1
            pass
        elif priority[source1] > priority[source0]:
            #Replace 1 with 0
            event_dict[event_id1] = event_dict[event_id0]
            event_array['source'][i1] = source0
        else: #Identical, replace 1 with 0 (keep earliest)
            event_dict[event_id1] = event_dict[event_id0]
            event_array['source'][i1] = source0
        #end if
        drop_rows[i0] = True
    #end for
    
    return event_array[~drop_rows], event_dict
#end func   

def filter_catalogue(catalogue, thr_time=15, thr_dist=0.5, mag_thr=5.5, 
                     io_mag_thr=5.0, aus_azim_gap_thr=180, 
                     aus_bbox=[105, -40, 160, -13],
                     io_bbox_1=[15, -60, 90, 30],
                     io_bbox_2=[90, -60, 180, -45]):
    """
    Function to filter a catalogue of events and to remove duplicates. If 
    events occur within thr_time and thr_dist of each other, they are classed 
    as duplicates.
    
    
    Parameters
    ----------
    catalogue : obspy.core.event Catalog object
        Catalogue to be filtered.
        
    thr_time : float (optional)
        Time threshold to classify events as duplicates (seconds).
        
    thr_dist : float (optional)
        Distance threshold to classify events as duplicates (degrees).
        
    mag_thr : float (optional)
        Minimum magnitude to keep events for.
        
    io_mag_thr : float (optional)
        Minimum magnitude for which to keep events in the Indian Ocean bounding 
        boxes.
        
    aus_azim_gap_thr : float (optional)
        Maximum azimuthal gap threshold for which to keep events in the 
        Australian bounding box.
        
    aus_bbox : list of float (optional)
        Minimum longitude, minimum latitude, maximum longitude, and maximum 
        latitude for the Australian bounding box.
        
    ind_oc_bbox_1 : list of float (optional)
        Minimum longitude, minimum latitude, maximum longitude, and maximum 
        latitude for the 1st Indian Ocean bounding box.
        
    ind_oc_bbox_2 : list of float (optional)
        Minimum longitude, minimum latitude, maximum longitude, and maximum 
        latitude for the 2nd Indian Ocean bounding box.
        
    
    Returns
    -------
    catalogue_filt : obspy.core.event Catalog object
        Filtered catalogue.
        
        
    """
    
    params = {'mag_thr': mag_thr, 'io_mag_thr': io_mag_thr, 
              'aus_azim_gap_thr': aus_azim_gap_thr, 'aus_bbox': aus_bbox, 
              'io_bbox_1': io_bbox_1, 'io_bbox_2': io_bbox_2}
    event_dict = {event.resource_id.id: event for event in catalogue}
    
    GG = list()
    ISC = list()
    GA = list()
    EHB = list()
    USGS = list()
    other = list()
    
    for event in catalogue:
        event_id = event.resource_id.id
        try:
            origin_time = event.preferred_origin().time.timestamp
            lon = event.preferred_origin().longitude
            lat = event.preferred_origin().latitude
            magval = event.preferred_magnitude().mag
            azim_gap = event.preferred_origin().quality.azimuthal_gap
        except:
            magval = 0.0
            for magnitude in event.magnitudes:
                magval = np.max([magval, magnitude.mag])
            #end for
            if event.origins[0].quality is None:
                azim_gap = 999.0
            else:
                azim_gap = event.origins[0].quality.azimuthal_gap
            #end if
            origin_time = event.origins[0].time.timestamp
            lon = event.origins[0].longitude
            lat = event.origins[0].latitude
        #end try
        source = get_source(event.resource_id.id)
        
        if source == 'GG':
            GG.append((event_id, source, origin_time, lon, lat, magval, 
                       azim_gap, 999.0, 999.0))
        elif source == 'ISC':
            ISC.append((event_id, source, origin_time, lon, lat, magval, 
                        azim_gap, 999.0, 999.0))
        elif source == 'GA':
            GA.append((event_id, source, origin_time, lon, lat, magval, 
                       azim_gap, 999.0, 999.0))
        elif source == 'EHB':
            EHB.append((event_id, source, origin_time, lon, lat, magval, 
                        azim_gap, 999.0, 999.0))
        elif source == 'USGS':
            USGS.append((event_id, source, origin_time, lon, lat, magval, 
                         azim_gap, 999.0, 999.0))
        else:
            other.append((event_id, source, origin_time, lon, lat, magval, 
                          azim_gap, 999.0, 999.0))
        #end if
    #end for
    
    dtype=[('event_id', 'U50'), ('source', 'U50'), ('time', float), 
           ('lon', float), ('lat', float), ('mag', float), 
           ('azim_gap', float), ('timediff', float), ('angdist', float)]
    
    GG = np.array(GG, dtype=dtype)
    ISC = np.array(ISC, dtype=dtype)
    GA = np.array(GA, dtype=dtype)
    EHB = np.array(EHB, dtype=dtype)
    USGS = np.array(USGS, dtype=dtype)
    other = np.array(other, dtype=dtype)
    
    GA, event_dict = replace_GA_with_GG(np.sort(GA, order='time'), 
                                        np.sort(GG, order='time'), 
                                        event_dict, thr_time, thr_dist)
    
    ISC = filter_ISC(ISC, params)
    GA = filter_GA(GA, params)
    EHB = filter_EHB(EHB, params)
    USGS = filter_USGS(USGS, params)
    other = filter_other(other, params)
    
    event_array = np.sort(np.hstack([ISC, GA, EHB, USGS, other]), order='time')
    
    event_array_filt, event_dict_filt = \
        remove_duplicates(event_array, event_dict, thr_time, thr_dist)
        
    catalogue_filt = Catalog(events=None)
    for event in event_array_filt:
        event_id = event['event_id']
        event_dict_filt[event_id].creation_info.author = event['source']
        catalogue_filt.append(event_dict_filt[event_id])
    #end for
    
    return catalogue_filt
#end func
        
def get_source(event_id):
    """
    Function to determine the source an event came from based on its event ID.
    
    
    Parameters
    ----------
    event_id : string
        Name of event.
        
        
    Returns
    -------
    source : string
        Name of event source.
        
        
    """
    if event_id[:2] == 'GG':
        return 'GG'
    elif event_id[:10] == 'quakeml:us':
        return 'USGS'
    elif 'ISC' in event_id:
        return 'ISC'
    elif 'ga' in event_id:
        if 'engdahl' in event_id:
            return 'EHB'
        else:
            return 'GA'
        #end if
    else:
        return 'other'
    #end if
#end func    
    
def process():
    """
    Read in event catalogues, filter by various bounding box restrictions and
    other criteria, remove duplicates, and produce output event binary files.
    
    
    Criteria
    --------
    1) No filtering performed on Bob Engdahl catalogue.
    2) USGS and ISC events in bounding box 105 < lon < 160, -40 < lat < -13 
        have azimuthal gap less than 180 degrees.
    3) USGS and ISC events in bounding boxes 15 < lon < 90, -60 < lat < 30 and 
        90 < lon < 180, -60 < lat < -45 have magnitude above 5.0.
    4) USGS and ISC events outside the areas described above have magnitude 
        above 5.5.
    5) GA/GG events have azimuthal gap less than 180 degrees.
        
    
    Arguments
    ---------
    -GA, --GA_path : string
        Directory in which GA event .xml files are stored.
        
    -EHB, --EHB_path : string
        Directory in which Engdahl event .xml files are stored.
        
    -ISC, --ISC_path : string
        Directory in which ISC event .xml files are stored.
        
    -PP, --preprocessed_path : string
        The output directory for Read_XML_Multiprocessing.py.
        
    -GG, --GG_path : string
        Directory in which Gary Gibson event .csv files are stored.
        
    -USGS, --USGS_path : string
        Directory in which USGS event .csv files are stored.
        
    -O, --output_path : string
        Directory in which to save output files.
        
    --read_from_preprocessed : boolean
        True if xml files have been pre-processed using 
        Read_XML_Multiprocessing.py, in which case .pkl files are read from 
        preprocessed_path. If else, .xml files are read from GA_path and/or 
        ISC_path and/or EHB_path.
        
        
    Usage
    -----
    python3 Filter_Catalogues.py --preprocessed_path .../preprocessed_path/ 
    --GG_path .../GG_path/ --output_path .../output_path/ 
    --read_from_preprocessed True/False

    python3 Filter_Catalogues.py --GA_path .../GA_path/ --ISC-path 
    .../ISC_path/ --EHB_path .../EHB_path/ --GG_path .../EHB_path/ --USGS_path 
    .../USGS_path/ --output_path .../output_path/ 
    
    
    """
    #print('Beginning process')
    with open('out.txt', 'a') as temp:
        temp.write(str('Beginning process' + '\n'))
    parser = argparse.ArgumentParser(description='Catalogue merger')

    parser.add_argument("-EHB", "--EHB_path", type=str, default="")
    parser.add_argument("-ISC", "--ISC_path", type=str, default="")
    parser.add_argument("-GA", "--GA_path", type=str, default="")
    parser.add_argument("-GG", "--GG_path", type=str, default="")
    parser.add_argument("-USGS", "--USGS_path", type=str, default="")
    parser.add_argument("-PP", "--preprocessed_path", type=str, default="")
    parser.add_argument("-O", "--output_path", type=str, default="")
    parser.add_argument("--read_from_preprocessed", type=bool, default=False)
    args = parser.parse_args()
    #print('Parsed arguments')
    with open('out.txt', 'a') as temp:
        temp.write(str('Parsed arguments' + '\n'))
    ISC_path = args.ISC_path
    GA_path = args.GA_path
    EHB_path = args.EHB_path
    USGS_path = args.USGS_path
    GG_path = args.GG_path
    preprocessed_path = args.preprocessed_path
    output_path = args.output_path
    read_from_preprocessed = args.read_from_preprocessed
    
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    #end if
    
    t = [time.time()]
    
    catalogues = list()
    if read_from_preprocessed:
        #print('Reading preprocessed catalogues')
        with open('out.txt', 'a') as temp:
            temp.write(str('Reading preprocessed catalogues' + '\n'))
        preprocessed_cat = read_preprocessed(preprocessed_path)
        catalogues.append(preprocessed_cat)
        t.append(time.time())
        #print('Done reading preprocessed catalogues, time =', t[-1] - t[-2])
        with open('out.txt', 'a') as temp:
            temp.write(str('Done reading preprocessed catalogues, time = ' \
                           + str(t[-1] - t[-2]) + '\n'))
    else:
        #print('Reading ISC catalogue')
        with open('out.txt', 'a') as temp:
            temp.write(str('Reading ISC catalogue' + '\n'))
        ISC_cat = read_event_xml_files(ISC_path)
        catalogues.append(ISC_cat)
        t.append(time.time())
        #print('Done reading ISC catalogue, time =', t[-1] - t[-2])
        with open('out.txt', 'a') as temp:
            temp.write(str('Done reading pISC catalogues, time = ' \
                           + str(t[-1] - t[-2]) + '\n'))
        
        #print('Reading GA catalogue')
        with open('out.txt', 'a') as temp:
            temp.write(str('Reading GA catalogue' + '\n'))
        GA_cat = read_event_xml_files(GA_path)
        catalogues.append(GA_cat)
        t.append(time.time())
        #print('Done reading GA catalogue, time =', t[-1] - t[-2])
        with open('out.txt', 'a') as temp:
            temp.write(str('Done reading GA catalogues, time = ' \
                           + str(t[-1] - t[-2]) + '\n'))
        
        #print('Reading Engdahl catalogue')
        with open('out.txt', 'a') as temp:
            temp.write(str('Reading Engdahl catalogue' + '\n'))
        EHB_cat = read_event_xml_files(EHB_path)
        catalogues.append(EHB_cat)
        t.append(time.time())
        #print('Done reading Engdahl catalogue, time =', t[-1] - t[-2])
        with open('out.txt', 'a') as temp:
            temp.write(str('Done reading Engdahl catalogues, time = ' \
                           + str(t[-1] - t[-2]) + '\n'))
            
        #print('Reading USGS catalogue')
        with open('out.txt', 'a') as temp:
            temp.write(str('Reading USGS catalogue' + '\n'))
        USGS_cat = read_event_xml_files(USGS_path)
        catalogues.append(USGS_cat)
        t.append(time.time())
        #print('Done reading USGS catalogue, time =', t[-1] - t[-2])
        with open('out.txt', 'a') as temp:
            temp.write(str('Done reading USGS catalogues, time = ' \
                           + str(t[-1] - t[-2]) + '\n'))
    #end if
    
    #print('Reading Gary Gibson catalogue')
    with open('out.txt', 'a') as temp:
        temp.write(str('Reading Gary Gibson catalogue' + '\n'))
    GG_cat = read_GG_csv_files(GG_path)
    catalogues.append(GG_cat)
    t.append(time.time())
    #print('Done reading Gary Gibson catalogue, time =', t[-1] - t[-2])
    #print('Time spent reading files =', t[-1] - t[0])
    with open('out.txt', 'a') as temp:
        temp.write(str('Done reading Gary Gibson catalogues, time = ' \
                       + str(t[-1] - t[-2]) + '\n'))
        temp.write(str('Time spent reading files = ' \
                       + str(t[-1] - t[0]) + '\n'))
    
    catalogue = merge_catalogues(catalogues)
    t.append(time.time())
    #print('Time spent merging catalogues =', t[-1] - t[-2])
    with open('out.txt', 'a') as temp:
        temp.write(str('Time spent merging catalogues = ' \
                       + str(t[-1] - t[-2]) + '\n'))
    
    catalogue_filt = filter_catalogue(catalogue)
    t.append(time.time())
    #print('Time spent filtering catalogue =', t[-1] - t[-2])
    with open('out.txt', 'a') as temp:
        temp.write(str('Time spent filtering catalogue = ' \
                       + str(t[-1] - t[-2]) + '\n'))
    
    write_catalogue_to_file(catalogue_filt, output_path)
    t.append(time.time())
    #print('Time spent writing catalogue to file =', t[-1] - t[-2])
    with open('out.txt', 'a') as temp:
        temp.write(str('Time spent writing to file = ' \
                       + str(t[-1] - t[-2]) + '\n'))
    return
#end func
    
if __name__ == "__main__":
    process()
#end if