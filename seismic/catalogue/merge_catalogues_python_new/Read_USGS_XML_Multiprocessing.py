"""
Description
-----------
This script is almost identical to Read_USGS_XML_Multiprocessing.py, however 
since it is impossible to filter data from the USGS ftp server, a filtering 
routine is included to remove events which don't fit certain bounding box 
restrictions (see Filter_Catalogue.py description).

Developer: Lachlan Adams 
Contact: lachlan.adams@ga.gov.au or lachlan.adams.1996@outlook.com

"""

from mpi4py import MPI
from obspy import Catalog, read_events
import numpy as np
import os, glob, pickle, argparse

def filter_USGS(arr, params):
    """
    Filter catalogue to only save events which fulfil criteria specified in
    params.
    
    
    Parameters
    ----------
    arr : numpy.ndarray
        Structured array containing event information.
        
    params: dictionary
        Dictionary containing filtering criteria.
        
        
    Returns
    -------
    arr : numpy.ndarray
        Structured array containing filtered event information.
        
        
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
    
def filter_catalogue(catalogue, thr_time=15, thr_dist=0.5, mag_thr=5.5, 
                     io_mag_thr=5.0, aus_azim_gap_thr=180, 
                     aus_bbox=[105, -40, 160, -13],
                     io_bbox_1=[15, -60, 90, 30],
                     io_bbox_2=[90, -60, 180, -45]):
    """
    Function to filter USGS catalogue so redundant/unusable events are 
    excluded. 
    
    
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
    
    USGS = list()
    
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
        
        USGS.append((event_id, origin_time, lon, lat, magval, azim_gap))
    #end for
    
    dtype=[('event_id', 'U50'), ('time', float), ('lon', float), 
           ('lat', float), ('mag', float), ('azim_gap', float)]
    
    USGS = np.array(USGS, dtype=dtype)
    USGS = filter_USGS(USGS, params)
    catalogue_filt = Catalog(events=None)
    for event in USGS:
        event_id = event['event_id']
        catalogue_filt.append(event_dict[event_id])
    #end for
    
    return catalogue_filt
#end func

def read_event_xml_files(files, output_path, split_smaller, rank):
    """
    Function to read event information from xml files. First tries to read in
    SC3ML format, and if this fails, attempts to read in QuakeML format.
    
    
    Parameters
    ----------
    files : list of string
        Files to read in.
        
    output_path : string
        Path in which to save output.
        
    split_smaller : boolean
        Flag to tell algorithm to make one output file for all inputs, or one
        output file for each input.
        
    rank : integer
        Rank of the processor, used to create output filename.
        
    
    """
    
    if len(files) == 0:
        return
    #end if
    
    if not split_smaller:
        catalogue = Catalog(events = None)
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
                if len(event.origins) == 0:
                    continue
                #end if
                if event.preferred_origin() is None:
                    event.preferred_origin_id = event.origins[0].resource_id
                #end if
                if event.preferred_magnitude() is None:
                    station_count_list = list()
                    for magnitude in event.magnitudes:
                        if magnitude.station_count is not None:
                            station_count_list.append(magnitude.station_count)
                        else:
                            station_count_list.append(0)
                        #end if
                    #end for
                    if len(station_count_list) == 0:
                        event.preferred_magnitude_id = None
                    else:
                        max_count = max(station_count_list)
                        max_ind = station_count_list.index(max_count)
                        event.preferred_magnitude_id = \
                            event.magnitudes[max_ind].resource_id
                    #end if
                #end if
                catalogue.append(event)
            #end for
        #end for
        catalogue = filter_catalogue(catalogue)
        with open(os.path.join(output_path, str(str(rank).zfill(3) + '.pkl')), 
                  'wb') as outp:
            pickle.dump(catalogue, outp, pickle.HIGHEST_PROTOCOL)
        #end with
    else:
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
            catalogue = Catalog(events = None)
            for event in newcat:
                if len(event.origins) == 0:
                    continue
                #end if
                if event.preferred_origin() is None:
                    event.preferred_origin_id = event.origins[0].resource_id
                #end if
                if event.preferred_magnitude() is None:
                    station_count_list = list()
                    for magnitude in event.magnitudes:
                        if magnitude.station_count is not None:
                            station_count_list.append(magnitude.station_count)
                        else:
                            station_count_list.append(0)
                        #end if
                    #end for
                    if len(station_count_list) == 0:
                        event.preferred_magnitude_id = None
                    else:
                        max_count = max(station_count_list)
                        max_ind = station_count_list.index(max_count)
                        event.preferred_magnitude_id = \
                            event.magnitudes[max_ind].resource_id
                    #end if
                #end if
                catalogue.append(event)
            #end for
            catalogue = filter_catalogue(catalogue)
            if len(catalogue) == 0:
                continue
            #end if
            savefile = str(file.split('/')[-1].split('.')[0] + '.pkl')
            with open(os.path.join(output_path, savefile), 'wb') as outp:
                pickle.dump(catalogue, outp, pickle.HIGHEST_PROTOCOL)
            #end with
        #end for
    return
#end func
    
def split_list(lst, npartitions):
    """
    Partition a list.
    
    
    Parameters
    ----------
    lst : list
        List to be partitioned.
        
    npartitions : integer
        Number of partitions.
        
        
    Returns
    -------
    lst_split : list
        List of lists, containing the partitions of lst.
        
        
    """
    lst_split = np.array_split(lst, npartitions)
    lst_split = [list(item) for item in lst_split]
    return lst_split
#end func    
    
def process():
    """
    Read event xml files from USGS input directories and pickle the 
    catalogues in binary files.
    
    
    Usage
    -----
    mpirun -np $number_of_processors python Read_XML_Multiprocessing.py
    --paths path_1 path_2 ... --output_path .../output_path/ 
    --split_smaller True/False
    
    
    """
    parser = argparse.ArgumentParser(description='XML file reader')
    
    parser.add_argument("--paths", nargs='+', default=[])
    parser.add_argument("--already_read_path", type=str, default=None)
    parser.add_argument("--output_path", type=str, default="")
    parser.add_argument("--split_smaller", type=bool, default=False)
    args = parser.parse_args()
    
    output_path = args.output_path
    paths = args.paths
    already_read_path = args.already_read_path
    split_smaller = args.split_smaller
    
    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()
    
    files = files_split = None
    
    if rank == 0:
        files = list()
        for path in paths:
            files.extend(glob.glob(os.path.join(path, '**/*.xml'), 
                                   recursive=True))
        #end for
        if already_read_path is not None:
            already_read_files = glob.glob(os.path.join(already_read_path, 
                                                        '*.pkl'))
            already_read_files = [file.split('/')[-1].split('.')[0] \
                                  for file in already_read_files]
            files = [file for file in files \
                     if file.split('/')[-1].split('.')[0] \
                     not in already_read_files]
        #end if
        files_split = split_list(files, nproc)
    #end if
    
    files = comm.scatter(files_split, root=0)
    read_event_xml_files(files, output_path, split_smaller, rank)
#end func
    
if __name__ == '__main__':
    process()
#end func