"""
Description
-----------
This script is used to shorten the time required to read in catalogues from 
.xml files. On multiple processors, catalogues are read in and merged, then 
output as a .pkl file. The output is to be read by Filter_Catalogues.py.

Developer: Lachlan Adams 
Contact: lachlan.adams@ga.gov.au or lachlan.adams.1996@outlook.com

"""

from mpi4py import MPI
from obspy import Catalog, read_events
import numpy as np
import os, glob, pickle, argparse

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
    Read event xml files from various input directories and pickle the 
    catalogues in binary files.
    
    
    Arguments
    ---------
    paths : list of string
        List of directories to read event .xml files from.
        
    output_path : string
        Directory in which to save output binary files.
        
    split_smaller : boolean
        If true, one output file is created for each input file, or if not, one
        output file is created per processor.
        
    already_read_path : string
        If the process has already been run previously with the flag 
        'split_smaller', then the names of output files will be the same as the 
        input files (less file extension). If already_read_path is provided, 
        the algorithm will search the directory for these output files, so the 
        corresponding inputs are ignored.
    
    
    Usage
    -----
    mpirun -np $number_of_processors python Read_XML_Multiprocessing.py
    --paths path_1 path_2 ... --output_path .../output_path/ 
    --split_smaller True/False --already_read_path .../already_read_path/
    
    
    """
    parser = argparse.ArgumentParser(description='XML file reader')
    
    parser.add_argument("--paths", nargs='+', default=['.'])
    parser.add_argument("--already_read_path", type=str, default=None)
    parser.add_argument("--output_path", type=str, default=".")
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