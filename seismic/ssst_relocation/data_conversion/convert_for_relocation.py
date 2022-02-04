"""
Description
-----------
This script is used to convert the output of the catalogue compilation workflow
and the pick harvesting workflow into the format required by the event 
relocation and phase redefinition algorithm.

Developer: Lachlan Adams 
Contact: lachlan.adams@ga.gov.au or lachlan.adams.1996@outlook.com

"""

import argparse, glob, os, time
import numpy as np
from obspy import UTCDateTime

def import_from_csv(file, delimiter=',', newline=''):
    """
    Import a list of lists from a .csv file.
    
    
    Parameters
    ----------
    file : string
        Name of file
        
    delimiter : string
        Separator used between items
        
    newline : string
        New line character
        
        
    Returns
    -------
    lst : list
        List of lists of elements from file
        
        
    """
    import csv
    lst = list()
    with open(file, newline=newline) as csvfile:
        reader = csv.reader(csvfile, delimiter=delimiter)
        for row in reader:
            lst.append([item.strip() for item in row])
        #end for
    #end with
    return lst
    """
    with open(file, 'r', newline=newline) as csvfile:
        lines = csvfile.readlines()
        for line in lines:
            line = line.split(delimiter)
        #end for
    #end with
    return lines
    """
#end func
    
def convert_event_list_to_array(event_list, pick_array):
    """
    Convert the list of events used by pick harvester routine into the format
    required by the source specific station term relocation method.
    Stacks picks read from event_list under those already contained within
    pick_array.
    
    
    Parameters
    ----------
    event_list : list
        List with rows formatted as the following.
        event_rows = [#source, YYYY, MM, DD, hh, mm, ss, lon, lat, depth, 
                      n_phase, mb, ms, ml, mw, event_id, azim_gap, index]
        pick_rows = [sta, cha, loc, net, lon, lat, elev, phase, YYYY, MM, DD, 
                     hh, mm, ss, ang_dist]
    
    pick_array : numpy.ndarray
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
    
    
    Returns
    -------
    pick_array : numpy.ndarray
        Structured array with data type described above.
        
        
    """
    new_row_list = list()
    nph_for_events = {}
    
    snr = -1.0
    qualityMeasureCWT = -1.0
    domFreq = -1.0
    qualityMeasureSlope = -1.0
    bandIndex = 0
    nSigma = 0
    
    for row in event_list:
        if row[0][0] == '#':
            source, YYYY, MM, DD, hh, mm, ss, elon, elat, edepthkm, nph, \
                mb, ms, ml, mw, ev_id, azim_gap, int_id = row
            event_id = str('smi:local/' + str(int(int_id)))
            ecolat = 90.0 - float(elat)
            edepth = 1e3*float(edepthkm)
            origin_time = UTCDateTime(int(YYYY), int(MM), int(DD), int(hh), 
                                      int(mm), float(ss)).timestamp
            if edepth == '':
                edepth = 0
            #end if
            if float(mb) > 0:
                mag = mb
            elif float(ms) > 0:
                mag = ms
            elif float(ml) > 0:
                mag = ml
            elif float(mw) > 0:
                mag = mw
            else:
                mag = 0
            #end if
            i = 0
            nph_for_events[event_id] = i
        else:
            i = i + 1
            sta, cha, loc, net, slon, slat, selev, phase, YYYY, MM, DD, \
                hh, mm, ss, ang_dist = row
            if slon == '' or slat == '':
                continue
            #end if
            if selev == '':
                selev = 0
            #end if
            scolat = 90.0 - float(slat)
            ptt = 0
            tcor = 0
            residual = 0
            arrival_time = UTCDateTime(int(YYYY), int(MM), int(DD), int(hh), 
                                       int(mm), float(ss)).timestamp
            pick_id = str('smi:local/' + str(int_id) + '_pick_' + str(i))
            new_row = (event_id, pick_id, sta, net, cha, elon, ecolat, edepth,
                       origin_time, mag, slon, scolat, selev, phase, 
                       arrival_time, ptt, tcor, residual, snr, 
                       qualityMeasureCWT, domFreq, qualityMeasureSlope, 
                       bandIndex, nSigma)
            new_row_list.append(new_row)
            nph_for_events[event_id] = i
        #end if
    #end for
    new_pick_array = np.array(new_row_list, dtype=pick_array.dtype)
    pick_array = np.hstack([pick_array, new_pick_array])
    return pick_array, nph_for_events
#end func
    
def convert_pick_list_to_array(pick_list, pick_array, nph_for_events, phase,
                               output_path):
    """
    Converts the list of picks produced by the pick harvester routine into the
    format required by the source specific station term relocation method, and
    adds it to the list of already existing picks.
    
    
    Parameters
    ----------
    pick_list : list
        List of picks with the following format.
        pick_rows = [event_id, origin_time, mag, elon, elat, edepth (km), net,
                     sta, cha, arrival_time, slon, slat, az, baz, dist, 
                     residual, snr, qualityMeasureCWT, domFreq, 
                     qualityMeasureSlope, bandIndex, nsigma]
    
    pick_array : numpy.ndarray
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
    
    nph_for_events : dictionary
        Dictionary with keyword-value pair event_id: npicks for each event.
        Used to track the number of picks pertaining to each event so that a
        pick ID may be assigned (event_id + '_pick_' + (npick+1)).
        
    phase : string
        'P' or 'S'.
        
    
    Returns
    -------
    pick_array : numpy.ndarray
        Structured array with data type described above.
    
        
    """
    new_row_list = list()
    for row in pick_list:
        if row[0][0] == '#': continue
        ev_id, origin_time, mag, elon, elat, edepthkm, net, sta, cha, \
            arrival_time, slon, slat, az, baz, dist, residual, snr, \
            qualityMeasureCWT, domFreq, qualityMeasureSlope, bandIndex, \
            nSigma = row
        # event ID is stored as 'xxxxx.0' which requires conversion into a 
        # float before into an integer.
        event_id = str('smi:local/' + str(int(float(ev_id))))
        ecolat = 90.0 - float(elat)
        scolat = 90.0 - float(slat)
        ptt = 0
        tcor = 0
        selev = 0
        edepth = 1e3*float(edepthkm)
        try:
            nph_for_events[event_id] = nph_for_events[event_id] + 1
            pick_id = str(event_id + '_pick_' + str(nph_for_events[event_id]))
        except:
            with open(os.path.join(output_path, 'out.txt'), 'a') as file:
                file.write(str('Event ' + event_id + ' not found!\n'))
            continue
        #end try
        new_row = (event_id, pick_id, sta, net, cha, elon, ecolat, edepth,
                   origin_time, mag, slon, scolat, selev, phase, arrival_time, 
                   ptt, tcor, residual, snr, qualityMeasureCWT, domFreq, 
                   qualityMeasureSlope, bandIndex, nSigma)
        new_row_list.append(new_row)
    new_pick_array = np.array(new_row_list, dtype=pick_array.dtype)
    pick_array = np.hstack([pick_array, new_pick_array])
    return pick_array, nph_for_events
#end func
    
def filter_duplicates(pick_array):
    """
    Filters list of picks to remove duplicate items, where station and network
    codes match and arrival time differs by less than 1 second.
    
    
    Parameters
    ----------
    pick_array : numpy.ndarray
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
        
    
    Returns
    -------
    pick_array : numpy.ndarray
        Input pick array with duplicates removed.
        
        
    """
    pick_array_sorted = np.sort(pick_array, order='arrival_time')
    tdiff = np.hstack([np.array(0), 
                       pick_array_sorted['arrival_time'][1:] - \
                       pick_array_sorted['arrival_time'][:-1]])
    small_tdiff = tdiff < 1.0
    duplicates = np.zeros(len(pick_array_sorted)).astype(bool)
    inds = np.where(small_tdiff)[0]
    
    for i in inds:
        pick0 = pick_array_sorted[i-1]
        pick1 = pick_array_sorted[i]
        t1 = pick0['net'] == pick1['net']
        t2 = pick0['stat'] == pick1['stat']
        if t1 and t2:
            duplicates[i] = True
        #end if
    #end for
    pick_array_filt = pick_array_sorted[~duplicates]
    return pick_array_filt
#end func
    
def process():
    """
    Reads event list output by catalogue compilation workflow, and pick lists
    output by pick harvester. Converts these lists into the format required by
    the source specific station term method, and filters to remove duplicate
    picks if station codes and network codes match, where arrival times differ
    by less than one second.
    
    
    Arguments
    ---------
    event_file : string
        Name of file output by catalogue compilation workflow.
        
    p_combined : string
        Directory containing list of extra P wave picks.
        
    s_combined : string
        Directory containing list of extra S wave picks.
        
    output_path : string
        Directory in which to place output file 'events.npy'.
        
        
    Usage
    -----
    python convert_for_relocation.py --event_file .../event_file.csv 
    --p_combined .../p_combined.txt --s_combined .../s_combined.txt 
    --output_path .../output_path/
    
    
    """
    parser = argparse.ArgumentParser(description='Catalogue converter')

    parser.add_argument("--event_file", type=str, required=True)
    parser.add_argument("--p_combined", type=str, default="")
    parser.add_argument("--s_combined", type=str, default="")
    parser.add_argument("--output_path", type=str, default=".")
    args = parser.parse_args()
    
    event_file = args.event_file
    p_combined_path = args.p_combined
    s_combined_path = args.s_combined
    output_path = args.output_path
    
    t0 = time.time()
    
    with open(os.path.join(output_path, 'out.txt'), 'w') as file:
        file.write(str('Reading event .csv file, start time =' + \
                       str(time.time()-t0) + '\n'))

    event_list = import_from_csv(event_file, delimiter=',')
    with open(os.path.join(output_path, 'out.txt'), 'a') as file:
        for row in event_list[:3]:
            file.write(str(str(row) + '\n'))
    
    p_combined_list = list()
    s_combined_list = list()
    
    with open(os.path.join(output_path, 'out.txt'), 'a') as file:
        file.write(str('Reading extra "P" picks, start time =' + \
                       str(time.time()-t0) + '\n'))
        
    if p_combined_path != "":
        for file in glob.glob(os.path.join(p_combined_path, '*.txt')):
            p_combined_list.extend(import_from_csv(file, delimiter=' '))
        #end for
    #end if
    
    with open(os.path.join(output_path, 'out.txt'), 'a') as file:
        for row in p_combined_list[:3]:
            file.write(str(str(row) + '\n'))
    
    with open(os.path.join(output_path, 'out.txt'), 'a') as file:
        file.write(str('Reading extra "S" picks, start time =' + \
                       str(time.time()-t0) + '\n'))
    
    if s_combined_path != "":
        for file in glob.glob(os.path.join(s_combined_path, '*.txt')):
            s_combined_list.extend(import_from_csv(file, delimiter=' '))
        #end for
    #end if
    
    with open(os.path.join(output_path, 'out.txt'), 'a') as file:
        for row in s_combined_list[:3]:
            file.write(str(str(row) + '\n'))
    
    with open(os.path.join(output_path, 'out.txt'), 'a') as file:
        file.write(str('Finished reading files, time =' + \
                       str(time.time()-t0) + '\n'))
    
    pick_array = np.empty(0, dtype = [('event_id', 'U20'), #Many characters
                                      ('pick_id', 'U30'), #Many characters
                                      ('stat', 'U10'), #Usually 3-5 characters
                                      ('net', 'U5'), #Usually 1-3 characters
                                      ('cha', 'U10'), #Usually 1-6 characters
                                      ('elon', 'single'), #Float between -180 and 180, high precision
                                      ('ecolat', 'single'), #Float between 0 and 180, high precision
                                      ('edepth', 'single'), #Float between 0 and 6.4*10^6, low precision
                                      ('origin_time', 'double'), #Float up to 10^10, high precision
                                      ('mag', 'half'), #Float between 0 and 10, low precision
                                      ('slon', 'single'), #Float between -180 and 180, high precision
                                      ('scolat', 'single'), #Float between 0 and 180, high precision
                                      ('selev', 'half'), #Float between 0 and 8000, low precision
                                      ('phase', 'U8'), #Up to 8 charaters
                                      ('arrival_time', 'double'), #Float up to 10^10, high precision
                                      ('ptt', 'single'), #Float between 0 and 10000, high precision
                                      ('tcor', 'half'), #Float between -10 and 10, high precision
                                      ('residual', 'half'), #Float between -10 and 10, high precision
                                      ('snr', 'half'), #Float with unknown range?
                                      ('qualityMeasureCWT', 'half'), #Float with unknown range?
                                      ('domFreq', 'half'), #Float with unknown range?
                                      ('qualityMeasureSlope', 'half'), #Float with unknown range?
                                      ('bandIndex', 'uint8'), #Integer between 0 and 10?
                                      ('nSigma', 'uint8')]) #Integer between 0 and 10?
    
    with open(os.path.join(output_path, 'out.txt'), 'a') as file:
        file.write(str('Converting event .csv file to array, start time =' + \
                       str(time.time()-t0) + '\n'))
    
    pick_array, nph_for_events = \
        convert_event_list_to_array(event_list, pick_array)
    
    with open(os.path.join(output_path, 'out.txt'), 'a') as file:
        file.write(str('Converting "P" picks to array, start time =' + \
                       str(time.time()-t0) + '\n'))
    
    pick_array, nph_for_events = \
        convert_pick_list_to_array(p_combined_list, pick_array, nph_for_events, 
                                   'P', output_path)
    
    with open(os.path.join(output_path, 'out.txt'), 'a') as file:
        file.write(str('Converting "S" picks to array, start time =' + \
                       str(time.time()-t0) + '\n'))
    
    pick_array, nph_for_events = \
        convert_pick_list_to_array(s_combined_list, pick_array, nph_for_events, 
                                   'S', output_path)
    
    with open(os.path.join(output_path, 'out.txt'), 'a') as file:
        file.write(str('Filtering duplicates, start time =' + \
                       str(time.time()-t0) + '\n'))
        
    pick_array = filter_duplicates(pick_array)
    
    with open(os.path.join(output_path, 'out.txt'), 'a') as file:
        file.write(str('Writing output file, start time =' + \
                       str(time.time()-t0) + '\n'))
    
    filename = os.path.join(output_path, 'events.npy')
    np.save(filename, pick_array)
    return
#end func
    
if __name__ == '__main__':
    process()
#end if