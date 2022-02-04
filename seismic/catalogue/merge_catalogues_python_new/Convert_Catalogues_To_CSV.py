"""
Description
-----------
This script is used to read the output of Filter_Catalogues and convert the 
catalogues to .csv format on multiple processors. 

Developer: Lachlan Adams 
Contact: lachlan.adams@ga.gov.au or lachlan.adams.1996@outlook.com

"""

from obspy.core.event import Catalog
from obspy import UTCDateTime
import copy, csv, glob, os, time, pickle, argparse
import numpy as np
import xml.etree.ElementTree as ET
from mpi4py import MPI

def read_station_xml(path):
    """
    Function to read station xml files to obtain station coordinates.
    
    
    Parameters
    ----------
    path : string
        Directory where station xml files are stored.
        
    
    Returns
    -------
    dct : dictionary
        Dictionary with keys network_code.station_code and entries [lontitude,
        latitude, elevation, recording_start_time, recording_end_time].
        
        
    """
    dct = {}
    
    if path == "":
        return dct
    else:
        files = glob.glob(os.path.join(path, '*.xml'))
        
        for file in files:
            filename = os.path.join(path, file)
            namespaces = dict([node for _, node \
                               in ET.iterparse(file, events=['start-ns'])])
            tree = ET.parse(filename)
            root = tree.getroot()
            ns = {'prefix':namespaces['']}
            networks = root.findall('prefix:Network', ns)
            
            for network in networks:
                network_code = network.attrib['code']
                startDate = UTCDateTime(network.attrib['startDate'])
                if 'endDate' in network.attrib.keys():
                    endDate = UTCDateTime(network.attrib['endDate'])
                else:
                    endDate = UTCDateTime(3000, 1, 1)
                #end if
                stations = network.findall('prefix:Station', ns)
                for station in stations:
                    station_code = station.attrib['code']
                    lon = float(station.find('prefix:Longitude', ns).text)
                    lat = float(station.find('prefix:Latitude', ns).text)
                    elev = float(station.find('prefix:Elevation', ns).text)
                    dct['%s.%s'%(network_code, station_code)] = \
                        [lon, lat, elev, startDate, endDate]
                #end for
            #end for
        #end for
        return dct
    #end if
#end func

def read_station_kml(path):
    """
    Function to read station kml files to obtain station coordinates.
    
    
    Parameters
    ----------
    path : string
        Directory where station kml files are stored.
        
    
    Returns
    -------
    dct : dictionary
        Dictionary with keys network_code.station_code and entries [lontitude,
        latitude, elevation, recording_start_time, recording_end_time].
        
        
    """
    dct = {}
    if path == "":
        return dct
    else:
        files = glob.glob(os.path.join(path, '*.kml'))
        
        for file in files:
            filename = os.path.join(path, file)
            namespaces = dict([node for _, node \
                               in ET.iterparse(file, events=['start-ns'])])
            tree = ET.parse(filename)
            root = tree.getroot()
            ns = {'prefix':namespaces['']}
            networks = root.findall('prefix:Document', ns)
            
            for network in networks:
                network_code = 'IR'
                stations = network.findall('prefix:Placemark', ns)
                for station in stations:
                    station_code = station.find('prefix:name', ns).text
                    description = station.find('prefix:description', ns)
                    attrib = description.text.split('<br>')
                    attrib = {item.split(':')[0]: item.split(':')[1].strip() \
                              for item in attrib if ': ' in item}
                    keys = attrib.keys()
                    if 'Lon' not in keys or 'Lat' not in keys:
                        continue
                    #end if
                    lon = float(attrib['Lon'])
                    lat = float(attrib['Lat'])
                    if 'Elev' in keys:
                        elev = float(attrib['Elev'])
                    else:
                        elev = 0.0
                    #end if
                    if 'Opened' in keys:
                        startDate = UTCDateTime(attrib['Opened'])
                    else:
                        startDate = UTCDateTime(0)
                    #end if
                    if 'Closed' in keys:
                        endDate = UTCDateTime(attrib['Closed'])
                    else:
                        endDate = UTCDateTime(3000, 1, 1)
                    #end if
                    dct['%s.%s'%(network_code, station_code)] = \
                        [lon, lat, elev, startDate, endDate]
                #end for
            #end for
        #end for
    #end if
    return dct
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

def split_list(lst, npartitions):
    """
    Partition a list.
    Preserves ordering of list making sums of indices in each approximately 
    equal, e.g. 
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] -> [[1, 4, 7, 10], [2, 5, 8], [3, 6, 9]].
    This was designed to distribute files between processors such that the sum
    of file sizes sent to each processor is approximately equal.
    
    
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
    temp = copy.deepcopy(lst)
    split_lst = [[] for n in range(npartitions)]
    for i in range(len(temp)):
        split_lst[i % npartitions].append(temp.pop(0))
    #end for
    return split_lst
#end func
    
def read_preprocessed(path):
    """
    Function to read event information already ingested from xml files. 
    Information is assumed to be stored in .pkl files.
    
    
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
        for file in files:
            with open(file, 'rb') as inp:
                newcat = pickle.load(inp)
            catalogue.extend(list(newcat))
        #end for
    #end if
    
    return catalogue
#end func
    
def catalogue_to_list(catalogue, station_dict, IR_station_dict):
    """
    Function to convert an obspy.core.event Catalog object into a list. 
        
    
    Parameters
    ----------
    catalogue : obspy.core.event Catalog object
        Catalogue to be converted into list.
        
    station_dict : dictionary
        Dictionary of station information with keys network_code.station_code 
        and entries [lontitude, latitude, elevation, recording_start_time, 
        recording_end_time].
        
    IR_station_dict : dictionary
        Dictionary of station information for the International Registry, with 
        keys network_code.station_code and entries [lontitude, latitude, 
        elevation, recording_start_time, recording_end_time]. Used if station
        is not found in station_dict.
        
        
    Returns
    -------
    lst : list
        List of events from catalogue. Format is:
            event_rows = [#source, YYYY, MM, DD, hh, mm, ss, lon, lat, depth, 
                          n_phase, mb, ms, ml, mw, event_id, azim_gap]
            pick_rows = [event_id, sta, cha, loc, net, lon, lat, elev, phase, 
                         YYYY, MM, DD, hh, mm, ss, ang_dist]
        
        
    """
    
    lst = list()
    
    for event in catalogue:
        author = event.creation_info.author
        origin = event.preferred_origin()
        if origin is None:
            origin = event.origins[0]
        #end if
        if origin.depth is None:
            origin.depth = 0.0
        #end if
        picks = event.picks
        if hasattr(origin.quality, 'azimuthal_gap'):
            azim_gap = origin.quality.azimuthal_gap
            if azim_gap is None:
                azim_gap = -999.0
            #end if
        else:
            azim_gap = -999.0
        #end if
        mb = ms = ml = mw = -999.0
        for mag in event.magnitudes:
            if mag.magnitude_type == 'mb' or mag.magnitude_type == 'MB' or \
                mag.magnitude_type == 'Mb' or mag.magnitude_type == 'mB':
                mb = mag.mag
            elif mag.magnitude_type == 'ms' or mag.magnitude_type == 'MS' or \
                mag.magnitude_type == 'Ms' or mag.magnitude_type == 'mS':
                ms = mag.mag
            elif mag.magnitude_type == 'ml' or mag.magnitude_type == 'ML' or \
                mag.magnitude_type == 'Ml' or mag.magnitude_type == 'mL':
                ml = mag.mag
            elif mag.magnitude_type == 'mw' or mag.magnitude_type == 'MW' or \
                mag.magnitude_type == 'Mw' or mag.magnitude_type == 'mW':
                mw = mag.mag
            #end if
        #end for
        lst.append([str('#' + author), origin.time.year, 
                origin.time.month, origin.time.day, origin.time.hour, 
                origin.time.minute, 
                origin.time.second + origin.time.microsecond/1e6, 
                origin.longitude, origin.latitude, origin.depth, 
                len(picks), mb, ms, ml, mw, str(event.resource_id), azim_gap])
        for pick in picks:
            net_code = pick.waveform_id.network_code
            stat_code = pick.waveform_id.station_code
            lon, lat, elev, net_code, dist = \
                find_station_information(stat_code, net_code, origin.longitude, 
                                         origin.latitude, origin.time, 
                                         station_dict, IR_station_dict)
            lst.append([str(event.resource_id), stat_code, 
                        pick.waveform_id.channel_code, 
                        pick.waveform_id.location_code, net_code, lon, lat, 
                        elev, pick.phase_hint, pick.time.year, pick.time.month,
                        pick.time.day, pick.time.hour, pick.time.minute, 
                        pick.time.second + pick.time.microsecond/1e6, dist])
        #end for
    #end for
    return lst
#end for
    
def find_station_information(stat, net, lon, lat, ot, station_dict, 
                             IR_station_dict):
    """
    Function to find the coordinates of a recording station. If a match 
    isn't found, all stations with matching station code are searched. If more 
    than one is found, the station which has a recording time overlapping with 
    the origin time of the event is used. If multiple matches are found, the 
    station which is closest to the event is used.
    
    
    Parameters
    ----------
    stat : string
        Station code.
        
    net : string
        Network code.
        
    lon : float
        Event longitude.
        
    lat : float
        Event latitude.
    
    ot : float
        Origin time.
        
    station_dict : dictionary
        Dictionary of station information with keys network_code.station_code 
        and entries [lontitude, latitude, elevation, recording_start_time, 
        recording_end_time].
        
    IR_station_dict : dictionary
        Dictionary of station information for the International Registry, with 
        keys network_code.station_code and entries [lontitude, latitude, 
        elevation, recording_start_time, recording_end_time]. Used if station
        is not found in station_dict.
        
    
    Returns
    -------
    slon : float
        Station longitude.
        
    slat : float
        Station latitude.
        
    selev : float
        Station elevation.
        
    net : string
        Network code. Needed in case it had to be updated to find a match in 
        the dictionary.
        
    dist : float
        Epicentral distance between source and receiver.
        
        
    """
    key = str(net + '.' + stat)
    if key in station_dict.keys():
        slon, slat, selev, _, _ = station_dict[key]
        dist = ang_dist(lon, lat, slon, slat)
        return slon, slat, selev, net, dist
    else:
        key = str('IR.' + stat)
        if key in IR_station_dict.keys():
            slon, slat, selev, _, _ = IR_station_dict[key]
            dist = ang_dist(lon, lat, slon, slat)
            keys = [key for key in station_dict.keys() \
                if key.split('.')[1] == stat]
            if len(keys) == 0:
                return slon, slat, selev, 'IR', dist
            else:
                rows = [station_dict[key] for key in keys]
                slons = np.array([row[0] for row in rows])
                slats = np.array([row[1] for row in rows])
                dists = ang_dist(slon, slat, slons, slats)
                if np.nanmin(dists) < 0.01:
                    ind = np.argwhere(dists == np.nanmin(dists))[0][0]
                    slon, slat, selev, _, _ = rows[ind]
                    net = keys[ind].split('.')[0]
                    dist = ang_dist(lon, lat, slon, slat)
                    return slon, slat, selev, net, dist
                else:
                    return slon, slat, selev, 'IR', dist
                #end if
            #end if
        else:
            slon = slat = selev = None
            dist = 999.0
            return slon, slat, selev, net, dist
        #end if
    #end if
    """
    else:
        keys = [key for key in station_dict.keys() \
                if key.split('.')[1] == stat]
        if len(keys) > 1: #More than one matching record
            rows = [station_dict[key] + [key] for key in keys]
            rows_tw = [row for row in rows \
                       if row[3].timestamp < ot and row[4].timestamp > ot]
            if len(rows_tw) == 0: #None with correct time window, use closest
                dists = [ang_dist(lon, lat, row[0], row[1]) for row in rows]
                row = [row for _, row in sorted(zip(dists, rows))][0]
                slon, slat, selev, _, _, key = row
                net = key.split('.')[0]
                dist = np.min(dists)
                return slon, slat, selev, net, dist
            elif len(rows_tw) == 1: #One match
                row = rows_tw[0]
                slon, slat, selev, _, _, key = row
                net = key.split('.')[0]
                dist = ang_dist(lon, lat, slon, slat)
                return slon, slat, selev, net, dist
            else: #More than one match, use closest
                dists = [ang_dist(lon, lat, row[0], row[1]) for row in rows_tw]
                row = [row for _, row in sorted(zip(dists, rows_tw))][0]
                slon, slat, selev, _, _, key = row
                net = key.split('.')[0]
                dist = np.min(dists)
                return slon, slat, selev, net, dist
            #end if
        elif len(keys) == 1:
            key = keys[0]
            net = key.split('.')[0]
            slon, slat, selev, _, _ = station_dict[key]
            dist = ang_dist(lon, lat, slon, slat)
            return slon, slat, selev, net, dist
        else:
            key = str('IR.' + stat)
            if key in IR_station_dict.keys():
                slon, slat, selev, _, _ = IR_station_dict[key]
                dist = ang_dist(lon, lat, slon, slat)
                return slon, slat, selev, 'IR', dist
            else:
                slon = slat = selev = None
                dist = 999.0
                return slon, slat, selev, net, dist
            #end if
        #end if
    #end if
    """
#end func
    
def write_list_to_csv(lst, path, rank, nproc):
    """
    Function to save a list of events as a .csv file. The prameters 'rank' and 
    'nproc' are used to create a filename.
    
    
    Parameters
    ----------
    lst : list
        List of events to write to file.
        
    path : string
        Directory in which to save file.
        
    rank : integer
        Rank of the processor.
        
    nproc : integer
        Number of processors in use.
    
    
    """
    i = 0
    newlst = list()
    for row in lst:
        if row[0][0] != '#':
            row.pop(0)
            newlst.append(row)
        else:
            if row[10] == 0:
                continue
            else:
                i += 1
            #end if
            newlst.append(row + [i])
        #end if
    #end for
    ndigits = int(np.ceil(np.log(nproc)/np.log(10)))
    savefile = os.path.join(path, str('events_' + str(rank).zfill(ndigits) + \
                                      '.csv'))
    if os.path.exists(savefile):
        os.remove(savefile)
    #end if
    with open(savefile, 'w+', newline ='') as file:
        write = csv.writer(file)
        write.writerows(newlst)
    #end with
#end func
    
def process():
    """
    Read obspy.core.event catalogues from .pkl binary files, convert the
    catalogues into a list format, add station information for each pick, and
    save the catalogue in .csv format, on multiple processors.
    Format for the output is:
        event_rows = [#source, YYYY, MM, DD, hh, mm, ss, lon, lat, depth, 
                          n_phase, mb, ms, ml, mw, event_id, azim_gap]
        pick_rows = [event_id, sta, cha, loc, net, lon, lat, elev, phase, 
                         YYYY, MM, DD, hh, mm, ss, ang_dist]
    
    
    Arguments
    ---------
    preprocessed_path : string
        Directory containing .pkl files to be converted.
        
    station_path : string
        Directory containing station .xml and/or .kml files, containing station
        coordinates.
        
    output_path : string
        Directory in which to save output .csv files.
        
    
    Usage 
    -----
    mpirun -np $number_of_processors python Convert_Catalogues_To_CSV.py
    --preprocessed_path .../preprocessed_path --station_path .../station_path
    --output_path .../output_path
    
    
    """
    parser = argparse.ArgumentParser(description='Catalogue merger')

    parser.add_argument("-PP", "--preprocessed_path", type=str, default=".")
    parser.add_argument("-ST", "--station_path", type=str, default=".")
    parser.add_argument("-O", "--output_path", type=str, default=".")
    args = parser.parse_args()
    
    preprocessed_path = args.preprocessed_path
    station_path = args.station_path
    output_path = args.output_path
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nproc = comm.Get_size()
    
    files_split = station_dict = IR_station_dict = None
    t = [time.time()]
    
    if rank == 0:
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        #end if
        
        #print('Reading station information')
        with open('out2.txt', 'a') as temp:
            temp.write(str('Reading station information' + '\n'))
        #end with
        station_dict = read_station_xml(station_path)
        IR_station_dict = read_station_kml(station_path)
        t.append(time.time())
        
        #print('Reading preprocessed catalogues')
        with open('out2.txt', 'a') as temp:
            temp.write(str('Reading station information' + '\n'))
        #end with
        files = glob.glob(os.path.join(preprocessed_path, '**/*.pkl'), 
                          recursive=True)
        filesizes = [os.path.getsize(file) for file in files]
        files_sorted = [x for _, x in sorted(zip(filesizes, files))]
        files_split = split_list(files_sorted, nproc)
    #end if
    
    comm.barrier()
    files_split = comm.scatter(files_split, root=0)
    station_dict = comm.bcast(station_dict, root=0)
    IR_station_dict = comm.bcast(IR_station_dict, root=0)
    #catalogue = Catalog(events=None)
    lst = list()
        
    if len(files_split) != 0:
        i = 0
        for file in files_split:
            i = i + 1
            if i % 10 == 0:
                with open('out2.txt', 'a') as temp:
                    temp.write(str('Reading file ' + str(i) + ' of ' + \
                                   str(len(files_split)) + ' on rank ' + \
                                   str(rank) + '\n'))
                #end with
            #end if
            with open(file, 'rb') as inp:
                newcat = pickle.load(inp)
            #end with
            lst.extend(catalogue_to_list(newcat, station_dict, 
                                         IR_station_dict))
        #end for
    #end if
    
    write_list_to_csv(lst, output_path, rank, nproc)
    t.append(time.time())
    if rank == 0:
        #print('Time spent writing list to csv =', t[-1] - t[-2])
        with open('out2.txt', 'a') as temp:
            temp.write(str('Time spent writing list to csv on rank ' + \
                           str(rank) + ' = ' + str(t[-1] - t[-2]) + '\n'))
        #end with
    #end if
    return
#end func
    
if __name__ == "__main__":
    process()
#end if
