"""
Description
-----------
This script is used to convert the output of the catalogue compilation workflow
and the pick harvesting workflow into the .xml format required by SeisComp3.

Developer: Lachlan Adams 
Contact: lachlan.adams@ga.gov.au or lachlan.adams.1996@outlook.com

"""

from obspy.core.event import Catalog, Event, Origin, Pick, WaveformStreamID, \
                             Arrival, Magnitude, CreationInfo, QuantityError, \
                             OriginQuality, OriginUncertainty
from obspy import UTCDateTime
import numpy as np
import argparse, configparser, os
from mpi4py import MPI

def azimuth(lon1, colat1, lon2, colat2, units='degrees'):
    """
    Function to calculate the azimuth from (lon1, lat1) to (lon2, lat2).
    
    
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
    azim : float
        Azimuth from (lon1, colat1) to (lon2, colat2).
        
        
    """
    if units == 'degrees':
        degrad = np.pi/180.0
        a = lon1*degrad
        b = colat1*degrad
        x = lon2*degrad
        y = colat2*degrad
    else:
        a = lon1
        b = colat1
        x = lon2
        y = colat2
    #end if
    azim = np.arctan(np.sin(x - a)/(np.sin(b)/np.tan(y) - \
                                    np.cos(b)*np.cos(x - a)))
    
    offset1 = np.pi*(colat1 < colat2)
    offset2 = np.pi*(np.logical_and(colat1 == colat2, b > np.pi/2.0))
    
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
    
def convert_array_to_catalogue(arr, event_ids):
    """
    Convert a structured array of picks into an obspy.core.event Catalog.
    
    
    Parameters
    ----------
    arr : numpy.ndarray
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
        
    event_ids : list
        List of event IDs for which to create the catalogue. These will be
        split between processors before this function is called if 
        multiprocessing is used.
        
        
    Returns
    catalogue : obspy.core.event Catalog object
        Catalogue of events with the picks in 'arr' that have event IDs in 
        'event_ids'.
        
        
    """
    catalogue = Catalog(events=None)
    for event_id in event_ids:
        event_arr = arr[arr['event_id'] == event_id]
        if len(event_arr) == 0:
            continue
        #end if
        event = convert_array_to_event(event_arr)
        catalogue.append(event)
    #end for
    return catalogue
#end func
    
def convert_array_to_event(arr): 
    """
    Convert an array of picks into an obspy.core.event Event object.
    
    
    Parameters
    ----------
    arr : numpy.ndarray
        Structured array containing all picks for the event with the following 
        data type.
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
    
    
    Returns
    -------
    event_object : obspy.core.event Event object
        Object containing origin, picks, magnitude, and other data for the 
        event.
        
        
    """
    event_id = arr['event_id'][0]
    origin_id = str(event_id + '_origin')
    magnitude_id = str(event_id + '_magnitude')
    lon = arr['elon'][0]
    lat = 90.0 - arr['ecolat'][0]
    depth = arr['edepth'][0]
    origin_time = UTCDateTime(arr['origin_time'][0])
    mag = arr['mag'][0]
    
    origin_list = list()
    arrival_list = list()
    pick_list = list()
    mag_list = list()
    azim_list = list()
    ecdist_list = list()
    
    for row in arr:
        pick_id = row['pick_id']
        arrival_id = pick_id.replace('pick', 'arrival')
        sta = row['stat']
        net = row['net']
        cha = row['cha']
        slon = row['slon']
        slat = 90.0 - row['scolat']
        phase = row['phase']
        arrival_time = UTCDateTime(row['arrival_time'])
        residual = row['residual']
        azim = azimuth(lon, lat, slon, slat)
        backazim = azimuth(slon, slat, lon, lat)
        ecdist = ang_dist(slon, slat, lon, lat)
        
        azim_list.append(azim)
        ecdist_list.append(ecdist)
        
        pick_obj = Pick(resource_id = pick_id,
                        time = arrival_time,
                        backazimuth = backazim,
                        creation_info = CreationInfo(agency_id = 'GA',
                                                     author = 'GA'),
                        waveform_id = WaveformStreamID(network_code=net,
                                                       station_code=sta,
                                                       channel_code=cha,
                                                       location_code=''),
                        phase_hint = phase,
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
        pick_list.append(pick_obj)
        
        arrival_obj = Arrival(resource_id = arrival_id,
                              pick_id = pick_id,
                              phase = phase,
                              azimuth = azim,
                              distance = ecdist,
                              creation_info = CreationInfo(agency_id = 'GA',
                                                           author = 'GA'),
                              time_correction = None,
                              takeoff_angle = None,
                              takeoff_angle_errors = None,
                              time_residual = residual,
                              horizontal_slowness_residual = None,
                              backazimuth_residual = None,
                              time_weight = None,
                              horizontal_slowness_weight = None,
                              backazimuth_weight = None,
                              earth_model_id = None,
                              comments = list())
        arrival_list.append(arrival_obj)
    #end for
    
    azim_arr = np.sort(azim_list)
    azim_gaps = np.hstack([azim_arr[1:] - azim_arr[:-1], 
                           np.array([azim_arr[0] + 360 - azim_arr[-1]])])
    azim_gap = np.max(azim_gaps)
    
    min_dist = np.min(ecdist_list)
    max_dist = np.max(ecdist_list)
        
    mag_obj = Magnitude(resource_id = magnitude_id,
                        mag = mag,
                        mag_errors = QuantityError(),
                        magnitude_type = 'M',
                        origin_id = str(event_id + '_origin'),
                        creation_info = CreationInfo(agency_id = 'GA',
                                                        author = 'GA'))
    mag_list.append(mag_obj)
    
    origin_obj = Origin(resource_id = origin_id,
                        time = origin_time,
                        longitude = lon,
                        latitude = lat,
                        depth = depth,
                        creation_info = CreationInfo(agency_id = 'GA',
                                                     author = 'GA'),
                        arrivals = arrival_list,
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
                                                minimum_distance=min_dist,
                                                maximum_distance=max_dist,
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
    origin_list.append(origin_obj)
        
    event_object = Event(resource_id = event_id,
                         preferred_origin_id = origin_id,
                         picks = pick_list,
                    	 amplitudes = list(),
                    	 focal_mechanisms = list(),
                    	 origins = origin_list,
                    	 magnitudes = mag_list,
                    	 station_magnitudes = list(),
                    	 creation_info = CreationInfo(agency_id='GA', 
                                                      author='GA'),
                    	 event_type = None,
                    	 event_type_certainty = 'known',
                    	 comments = list(),
                    	 event_descriptions = list())
        
    return event_object
#end func
    
def write_events_to_file(catalogue, path, output_format):
    """
    Write a catalogue of obspy.core.event objects each into a separate .xml 
    file in the format 'output_format', which may be 'SC3ML', 'QuakeML', or
    many others.
    
    
    Parameters
    ----------
    catalogue : obspy.core.event Catalog object
        Catalogue of events to save.
        
    path : string
        Output path.
        
    output_format : string
        Format for the output .xml files. May be 'SC3ML', 'QuakeML', or many 
        others.
    
    
    """
    for i in range(len(catalogue)):
        dtstring = str(catalogue[i].origins[0].time)
        dtstring = dtstring.replace(':', '-').split('.')[0]
        filename = os.path.join(path, str(dtstring + 'Z.xml'))
        catalogue[i:i+1].write(filename, format=output_format)
    #end for
#end func
    
def split_list_sorted(lst, npartitions):
    """
    Split a sorted list evenly into 'npartitions' partitions.
    E.g. [1, 2, 3, 4, 5, 6, 7, 8] -> [[1, 4, 7], [2, 5, 8], [3, 6]]
    
    
    Parameters
    ----------
    lst : list
        List to be split.
    
    npartitions : integer
        Number of partitions of 'lst' to be created.
        
    
    Returns
    -------
    new_lst : list
        List of the partitions of 'lst'.
        
        
    """
    new_lst = [list() for i in range(npartitions)]
    for i in range(len(lst)):
        new_lst[i % npartitions].append(lst[i])
    #end for
    return new_lst
#end func

def partition_by_event(picks, nproc):
    events = np.unique(picks['event_id'])
    dct = {}
    for i in range(len(events)):
        dct[events[i]] = i
    #end for
    n_picks = np.zeros(len(events)).astype(int)
    for i in range(len(picks)):
        n_picks[dct[picks['event_id'][i]]] -= 1
    #end for
    
    events = [item for _, item in sorted(zip(n_picks, events))]
    events_split = split_list_sorted(events, nproc)
    dct = {}
    
    for i in range(nproc):
        for event in events_split[i]:
            dct[event] = i
        #end for
    #end for
    
    partitions = np.array([dct[picks[i]['event_id']] \
                           for i in range(len(picks))])
    
    picks_split = [picks[partitions == i] for i in range(nproc)]
    return picks_split, events_split
#end if
    
def process():
    """
    Read in an array of picks, split between processors, convert the arrays 
    into obspy.core.event Catalog objects, and save each event to an output
    .xml file.
    
    
    Arguments
    ---------
    event_file : string
        Name of file containing array of picks.
        
    output_path : string
        Output path.
        
    output_format : string
        Format to use for the output .xml file. May be 'SC3ML', 'QuakeML', or 
        many others.
        
        
    Usage
    -----
    mpirun -np $number_of_processors python convert_array_to_catalogue.py
    --event_file .../event_file.npy --output_path .../output_path/ 
    --output_format <output_format>
    
    
    """
    parser = argparse.ArgumentParser(description='Catalogue converter')

    parser.add_argument("--event_file", type=str, required=True)
    parser.add_argument("--config_file", type=str, required=True)
    parser.add_argument("--output_path", type=str, default=".")
    parser.add_argument("--output_format", type=str, default='SC3ML')
    args = parser.parse_args()
    
    config = configparser.ConfigParser()
    config.read(args.config_file)
    config = config['default']
    ignore_temp_networks = config['ignore_temp_networks'] == 'True'
    if ignore_temp_networks:
        temp_networks = config['temp_networks'].split(', ')
    else:
        temp_networks = list()
    #end if
    
    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()
    
    events_split = None
    if rank == 0:
        file = args.event_file
        picks = np.load(file)
        picks_split, events_split = partition_by_event(picks, nproc)
        for i in range(nproc):
            np.save(os.path.join(args.output_path, '%s.npy'%str(i).zfill(3)), 
                    picks_split[i])
        #end for
    #end if
    
    events_split = comm.scatter(events_split, root=0)
    picks_split = np.load(os.path.join(args.output_path, 
                                       '%s.npy'%str(rank).zfill(3)))
    ind = ~np.isin(picks_split['net'], np.array(temp_networks))
    picks_split = picks_split[ind]
    
    os.remove(os.path.join(args.output_path, '%s.npy'%str(rank).zfill(3)))
    
    catalogue = convert_array_to_catalogue(picks_split, events_split)
    
    write_events_to_file(catalogue, args.output_path, args.output_format)
#end func
    
if __name__ == '__main__':    
    process()
#end if
