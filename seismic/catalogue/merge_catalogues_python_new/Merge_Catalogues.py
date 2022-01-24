"""
Description
-----------
This script is used to merge the .csv output of 
Convert_Catalogues_To_CSV.py, and re-order the events by date.

Developer: Lachlan Adams 
Contact: lachlan.adams@ga.gov.au or lachlan.adams.1996@outlook.com

"""

from obspy import UTCDateTime
import numpy as np
import argparse, csv, glob, os

def reorder_list_by_date(lst):
    """
    Takes a list of events and picks, extracts the event rows, sorts them by 
    origin time, and returns the sorted list with attached picks. Also converts
    depth to kilometres so it can be used by the pick harvester routine.
    
    
    Parameters
    ----------
    lst : list
        List of event information.
        event_rows = [#source, YYYY, MM, DD, hh, mm, ss, lon, lat, depth, 
                      n_phase, mb, ms, ml, mw, event_id, azim_gap, index]
        pick_rows = [sta, cha, loc, net, lon, lat, elev, phase, YYYY, MM, DD, 
                     hh, mm, ss, ang_dist]
        
    Returns
    -------
    lst : list
        Sorted list.
        
        
    """
    print('Reordering list')
    events = {}
    first = True
    event = time = None
    for row in lst:
        if row[0][0] == '#':
            if first:
                first = False
            else:
                flag = 0
                while flag == 0:
                    if str(time) in events.keys():
                        time = time + 1
                    else:
                        flag = 1
                    #end if
                #end while                         
                events[str(time)] = event
            #end if
            event = list()
            time = int(UTCDateTime(int(row[1]), int(row[2]), 
                                   int(row[3]), int(row[4]), 
                                   int(row[5]), float(row[6])).timestamp)
        #end if
        event.append(row)
    #end for
    
    keys = np.array([int(key) for key in events.keys()])
    sorted_keys = np.sort(keys)
    sorted_list = list()
    for key in sorted_keys:
        sorted_list.extend(events[str(key)])
    #end for
    
    i = 0
    for row in sorted_list:
        if row[0][0] == '#':
            i = i + 1
            row[-1] = i
            row[-3] = str(i) #str('smi:local/' + str(i))
            row[9] = str(float(row[9])/1e3)
        #end if
    #end for
            
    return sorted_list
#end func
    
def write_to_csv(filename, lst):
    """
    Writes the event list to a .csv file.
    
    
    Parameters
    ----------
    filename : string
        File name.
        
    lst : list
        List of event information.
        event_rows = [#source, YYYY, MM, DD, hh, mm, ss, lon, lat, depth, 
                      n_phase, mb, ms, ml, mw, event_id, azim_gap, index]
        pick_rows = [sta, cha, loc, net, lon, lat, elev, phase, YYYY, MM, DD, 
                     hh, mm, ss, ang_dist]
        
        
    """
    print('Writing to file')
    with open(filename, 'w+', newline='') as file:
        write = csv.writer(file)
        write.writerows(lst)
    #end with
#end func
    
def read_from_csv(path):
    """
    Reads event information from .csv files.
    
    
    Parameters
    ----------
    path : string
        Directory to read files from.
        
        
    Returns
    -------
    lst : list
        List of event information.
        event_rows = [#source, YYYY, MM, DD, hh, mm, ss, lon, lat, depth, 
                      n_phase, mb, ms, ml, mw, event_id, azim_gap, index]
        pick_rows = [sta, cha, loc, net, lon, lat, elev, phase, YYYY, MM, DD, 
                     hh, mm, ss, ang_dist]
        
        
    """
    lst = []
    files = glob.glob(os.path.join(path, '*.csv'))
    for filename in files:
        print('Reading from file', filename)
        with open(filename, newline='') as file:
            reader = csv.reader(file, delimiter=',')
            for row in reader:
                lst.append(row)
            #end for
        #end with
    return lst
#end func
    
def process():
    """
    Read input event files in list format, merge them, and produce a single 
    output. Format for the output is:
        event_rows = [#source, YYYY, MM, DD, hh, mm, ss, lon, lat, depth, 
                          n_phase, mb, ms, ml, mw, event_id, azim_gap]
        pick_rows = [event_id, sta, cha, loc, net, lon, lat, elev, phase, 
                         YYYY, MM, DD, hh, mm, ss, ang_dist]
        
    
    Arguments
    ---------
    input_path : string
        Input path.
        
    output_path : string
        Output path.
        
        
    Usage 
    -----
    python Merge_Catalogues.py --input_path .../input_path --output_path 
    .../output_path
    
    
    """
    parser = argparse.ArgumentParser(description='Catalogue merger')
    
    parser.add_argument("-I", "--input_path", type=str, default="")
    parser.add_argument("-O", "--output_file", type=str, default="")
    args = parser.parse_args()
    
    input_path = args.input_path
    output_filename = args.output_file
    
    print('Reading input files')
    lst = read_from_csv(input_path)
    
    print('Sorting event list')
    lst = reorder_list_by_date(lst)
    
    print('Writing output file')
    write_to_csv(output_filename, lst)
#end func
    
if __name__ == '__main__':
    process()
#end if