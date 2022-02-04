"""
Description
-----------
This script is used to perform phase redefinition for picks attributed to 
seismic events, and to perform event relocation using static or source-specific 
station terms.

Developer: Lachlan Adams 
Contact: lachlan.adams@ga.gov.au or lachlan.adams.1996@outlook.com

"""

import argparse, configparser, os, time, warnings
import numpy as np
from mpi4py import MPI
from Travel_Times import process_tt_tables, read_ellipcorr_table, \
                         predict_travel_times
from Station_Corrections import calculate_station_corrections
                       
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

def partition_by_station(picks, npartitions):
    """
    Partition pick array such that all picks for a station are placed within 
    the same partition.
    
    
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
        
    npartitions : integer
        Number of partitions.
        
    
    Returns
    -------
    picks_split : list of numpy.ndarray
        List containing the partitions of "picks" array.
        
    statnames_split : list
        List of the names of each station in a partition.
        
    
    """
    
    
    statnames = np.unique(picks['stat'])
    statnames_split = list(np.array_split(statnames, npartitions))
    dct = {}
    
    for i in range(npartitions):
        for stat in statnames_split[i]:
            dct[stat] = i
        #end for
    #end for
    
    partitions = np.array([dct[picks[i]['stat']] for i in range(len(picks))])
    
    picks_split = [picks[partitions == i] for i in range(npartitions)]
    return picks_split, statnames_split
#end func
    
def partition_by_event(picks, npartitions):
    """
    Partition pick array such that all picks for an event are placed within the 
    same partition. The list is partitioned so that each partition has 
    approximately the same cardinality. 
    
    
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
        
    npartitions : integer
        Number of partitions.
        
    
    Returns
    -------
    picks_split : list of numpy.ndarray
        List containing the partitions of "picks" array.
        
    events_split : list
        List of the names of each event in a partition.
        
        
    """
    
    
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
    events_split = split_list_sorted(events, npartitions)
    dct = {}
    
    for i in range(npartitions):
        for event in events_split[i]:
            dct[event] = i
        #end for
    #end for
    
    partitions = np.array([dct[picks[i]['event_id']] \
                           for i in range(len(picks))])
    
    picks_split = [picks[partitions == i] for i in range(npartitions)]
    return picks_split, events_split
#end if
    
def write(output_file, *args):
    string = ''.join((str(item) for item in args)) + '\n'
    with open(output_file, 'a') as file:
        file.write(string)
    print(string)
#end func

def process():
    """
    Process to perform source specific station term relocation of a set of 
    seismic events.
    

    Process
    -------
    1. Read and interpolate travel time tables, and ellipticity correction 
        coefficient tables.
    2. Read event/pick information file.
    3. Partition pick list evenly between multiple processors and compute 
        predicted travel times, and travel time residuals, then bring list back 
        to a single processor.
    4. Partition pick list by station names and compute source specific station 
        term corrections on multiple processors, then bring list back to a 
        single processor.
    5. Partition pick list by event names and relocate hypocentres, then bring 
        list back to a single processor.
    6. Find events which have an unstable hypocentre and write list of their 
        names to disk. These events are no longer used to compute the travel 
        time corrections.
        
    
    Arguments
    ---------
    config_file : string
        Name of configuration file to use.
    
    tt_table_path : string
        Directory containing pre-computed travel time tables.
        
    input_path : string
        Directory containing pick array save file and (if iteration > 0) list
        of events with unstable hypocentres.
        
    elcordir : string
        Directory containing ellipticity correction coefficient tables.
        
    iteration : integer
        Iteration number.
        
    output_path : string
        Directory to output updated pick array and list of events with unstable 
        hypocentres
        
    relocation_algorithm : string
        If relocation_algorithm is 'iloc' then the iLoc relocation algorithm is 
        used, or if not, a heavily modified version of the XCorLoc algorithm is 
        used.
        
        
    """
    
    
    # Parse arguments and config
    warnings.simplefilter("ignore", category=RuntimeWarning)
    
    parser = argparse.ArgumentParser(description='SSST relocator')

    parser.add_argument("--config_file", type=str, required=True)
    parser.add_argument("--tt_table_path", type=str, required=True)
    parser.add_argument("--input_path", type=str, required=True)
    parser.add_argument("--elcordir", type=str, required=True)
    parser.add_argument("--iteration", type=int, required=True)
    parser.add_argument("--output_path", type=str, default='.')
    parser.add_argument("--relocation_algorithm", type=str, default='x')
    
    """
    #Testing
    import sys
    sys.argv = ['/g/data/ha3/la8536/SSST/relocation/main.py', "--config_file", 
                '/g/data/ha3/la8536/SSST/relocation/relocation_config.txt',
                "--tt_table_path", '/g/data/ha3/la8536/SSST/TT_tables/',
                "--input_path", '/g/data/ha3/la8536/SSST/old/input_events/',
                "--elcordir", '/g/data/ha3/la8536/SSST/TT_tables/',
                "--iteration", '0',
                "--output_path", '/g/data/ha3/la8536/SSST/test/']
    sys.argv = ['/home/centos/ladams/ssst_relocation/relocation/main.py', 
                "--config_file", 
                '/home/centos/ladams/ssst_relocation/relocation/relocation_config.txt',
                "--tt_table_path", 
                '/home/centos/ladams/ssst_relocation/TT_tables/',
                "--input_path", 
                '/home/centos/ladams/ssst_relocation/output_events/ssst_iloc/',
                "--elcordir", '/home/centos/ladams/ssst_relocation/TT_tables/',
                "--iteration", '0', "--output_path", 
                '/home/centos/ladams/ssst_relocation/output_events/ssst_iloc/',
                "--relocation_algorithm", 'iloc']
    """
    
    args = parser.parse_args()
    config_file = args.config_file
    tt_table_path = args.tt_table_path
    input_path = args.input_path
    elcordir = args.elcordir
    iteration = args.iteration
    output_path = args.output_path
    relocation_algorithm = args.relocation_algorithm
    outfile = os.path.join(output_path, 'out.txt')

    config = configparser.ConfigParser()
    config.read(config_file)
    
    if relocation_algorithm == 'iloc':
        config = config['iLoc']
    else:
        config = config['default']
    #end if
    
    no_relocation = config['no_relocation'] == 'True'
    if iteration == 0:
        config['method'] = 'None'
        config['reloc_dlat'] = '1.0'
        config['reloc_nit'] = '10'
    #end if
    
    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()
    
    # Process travel time tables
    if rank == 0:
        t0 = time.time()
        write(outfile, 'Processing travel time tables, time = ', 
              time.time() - t0)
    #end if
    
    TT_dict, phase_list = process_tt_tables(tt_table_path)
    ellipcorr_dict, _ = read_ellipcorr_table(os.path.join(elcordir, 
                                                          'elcordir.tbl'))
    
    # Read pick information file
    unstable_events = None
    if rank == 0:
        write(outfile, 'Reading event information from file, time = ', 
              time.time() - t0)
        if iteration == 0:
            filename = os.path.join(input_path, 'picks.npy')
            picks = np.load(filename)
            unstable_events = list()
        elif iteration == 1:
            filename = os.path.join(input_path, 'picks.npy')
            picks = np.load(filename)
            filename = os.path.join(input_path, 'unstable_events.txt')
            with open(filename, 'r') as file:
                unstable_events = file.readlines()
            #end with
        else:
            filename = os.path.join(input_path, 
                                    str('picks' + str(iteration - 1) + '.npy'))
            picks = np.load(filename)
            filename = os.path.join(input_path, 'unstable_events.txt')
            with open(filename, 'r') as file:
                unstable_events = file.readlines()
            #end with
        #end if
        picks_split = list(np.array_split(picks, nproc))
        for i in range(nproc):
            np.save(os.path.join(output_path, '%s.npy'%str(i).zfill(3)), 
                    picks_split[i])
        #end for
        write(outfile, 'Predicting travel times, time = ', time.time() - t0)
    #end if
    
    # Calculate travel time residuals  
    unstable_events = comm.bcast(unstable_events, root=0)
    picks_split = np.load(os.path.join(output_path, 
                                       '%s.npy'%str(rank).zfill(3)))
    
    picks_split = predict_travel_times(picks_split, phase_list, TT_dict, 
                                       ellipcorr_dict, config)
    
    for i in range(nproc):
        if rank == i:
            np.save(os.path.join(output_path, '%s.npy'%str(rank).zfill(3)), 
                    picks_split)
            comm.barrier()
        #end if
    #end for
    
    # Gather data
    statnames_split = None
    if rank == 0:
        write(outfile, 'Completed predicting travel times, time = ', 
              time.time() - t0)
        picks_split = list()
        
        for i in range(nproc):
            picks_split.append(np.load(os.path.join(output_path, 
                                                    '%s.npy'%str(i).zfill(3))))
            os.remove(os.path.join(output_path, '%s.npy'%str(i).zfill(3)))
        #end for
        
        picks = np.hstack(picks_split)
        
        # Save picks before relocation
        np.save(os.path.join(output_path, 
                             str('picks' + str(iteration) + '_temp.npy')),
                picks)
        
        picks_split, statnames_split = partition_by_station(picks, nproc)
        
        for i in range(nproc):
            np.save(os.path.join(output_path, '%s.npy'%str(i).zfill(3)), 
                    picks_split[i])
        #end for
        write(outfile, 'Calculating travel time corrections, time = ', 
              time.time() - t0)
    #end if
    
    # Calculate and apply travel time corrections
    statnames_split = comm.scatter(statnames_split, root=0)
    picks_split = np.load(os.path.join(output_path, 
                                       '%s.npy'%str(rank).zfill(3)))
    
    picks_split = \
        calculate_station_corrections(statnames_split, picks_split, rank, 
                                      config, unstable_events=unstable_events)
    for i in range(nproc):
        if rank == i:
            np.save(os.path.join(output_path, '%s.npy'%str(rank).zfill(3)), 
                    picks_split)
            comm.barrier()
        #end if
    #end for
    
    # Gather data
    events_split = None
    if rank == 0:
        write(outfile, 'Finished calculating travel time corrections, time = ',
              time.time() - t0)
        picks_split = list()
        
        for i in range(nproc):
            picks_split.append(np.load(os.path.join(output_path, 
                                                    '%s.npy'%str(i).zfill(3))))
            os.remove(os.path.join(output_path, '%s.npy'%str(i).zfill(3)))
        #end for
        
        picks = np.hstack(picks_split)
        
        if no_relocation == True:
            write(outfile, 'Writing event information to file, time = ', 
                  time.time() - t0)
            
            filename = os.path.join(output_path, 
                                    str('picks' + str(iteration) + '.npy'))
            np.save(filename, picks)
        else:        
            picks_split, events_split = partition_by_event(picks, nproc)
            for i in range(nproc):
                np.save(os.path.join(output_path, '%s.npy'%str(i).zfill(3)), 
                        picks_split[i])
            #end for
            if relocation_algorithm == 'iloc':
                """
                If relocation algorithm is iloc:
                    1. Partition pick array by event ID between processors.
                    2. Push time corrections to database.
                    3. Save pick arrays to file.
                """
                from Relocation import push_time_corrections_to_database
                push_time_corrections_to_database(picks)
                write(outfile, 'Relocating events using iLoc, time = ', 
                      time.time() - t0)
            else:
                """
                If relocation algorithm is not iloc:
                    1. Partition pick array by event ID between processors.
                    2. Save pick arrays to file.
                """
                write(outfile, 'Relocating events, time = ', time.time() - t0)
            #end if
        #end if
    #end if 
    
    # Relocate events and calculate new residuals    
    if no_relocation == True:
        comm.barrier()
    else:
        events_split = comm.scatter(events_split, root=0)
        
        if relocation_algorithm == 'iloc':
            """
            If relocation algorithm is iloc:
                1. Compute new hypocentres using iloc.
                2. Extract hypocentres from database.
                3. Load pick save file.
                4. Update picks wiht new hypocentres.
            """
            from Relocation import compute_new_hypocentre_iloc, \
                                   update_hypocentres_from_database
            compute_new_hypocentre_iloc(events_split, output_path, config, 
                                        rank)
            comm.barrier()
            
            hypo_dict = None
            if rank == 0:
                from Relocation import extract_hypocentres_from_database
                hypo_dict = extract_hypocentres_from_database()
            #end if
            
            hypo_dict = comm.bcast(hypo_dict, root=0)
            picks_split = np.load(os.path.join(output_path, 
                                               '%s.npy'%str(rank).zfill(3)))
            
            picks_split, unstable_events = \
                update_hypocentres_from_database(events_split, picks_split, 
                                                 hypo_dict, config, 
                                                 unstable_events=\
                                                 unstable_events)
        else:
            """
            If relocation algorithm is not iloc:
                1. Load pick array save files.
                2. Compute new hypocentres.
                3. Store pick array in memory for calculation of residuals with
                    respect to new hypocentre.
            """
            from Relocation import compute_new_hypocentre
            picks_split = \
                np.load(os.path.join(output_path, '%s.npy'%str(rank).zfill(3)))
                
            picks_split, unstable_events = \
                compute_new_hypocentre(events_split, picks_split, TT_dict, 
                                       ellipcorr_dict, output_path, rank, 
                                       config, unstable_events=unstable_events)
        #end if
        
        picks_split = predict_travel_times(picks_split, phase_list, TT_dict, 
                                           ellipcorr_dict, config)
                
        for i in range(nproc):
            if rank == i:
                np.save(os.path.join(output_path, '%s.npy'%str(rank).zfill(3)), 
                        picks_split)
                comm.barrier()
            #end if
        #end for
        
        unstable_events = comm.gather(unstable_events, root=0)  
    #end if
    
    # Gather data
    if rank == 0 and no_relocation == False:
        write(outfile, 'Completed relocating events, time = ', 
              time.time() - t0)
        
        unstable_events = list(np.unique([item for lst in unstable_events \
                                              for item in lst]))
        picks_split = list()
        for i in range(nproc):
            picks_split.append(np.load( \
                os.path.join(output_path, '%s.npy'%str(i).zfill(3))))
            os.remove(os.path.join(output_path, '%s.npy'%str(i).zfill(3)))
            os.remove(os.path.join(output_path, 'out%s.txt'%str(i).zfill(3)))
        #end for
        picks = np.hstack(picks_split)
        
        write(outfile, 'Writing event information to file, time = ', 
              time.time() - t0)
        
        filename = os.path.join(output_path, 
                                str('picks' + str(iteration) + '.npy'))
        np.save(filename, picks)
        
        filename = os.path.join(output_path, 'unstable_events.txt')
        with open(filename, 'w') as file:
            for event in unstable_events:
                file.write(str(str(event) + '\n'))
            #end for
        #end with
    #end if
#end func

if __name__ == '__main__':
    process()
#end func