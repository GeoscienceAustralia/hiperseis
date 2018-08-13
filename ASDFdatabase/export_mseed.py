#!/bin/env python
"""
Description:
    Small utility for exporting mseed files from an asdf file in parallel.
References:

CreationDate:   06/08/18
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     04/12/17   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

from mpi4py import MPI
import os, sys

from os.path import join, exists
from collections import defaultdict
import numpy as np
from obspy import Stream, Trace, UTCDateTime
from ASDFdatabase.seisds import SeisDB
from xcorqc.correlator import Dataset
from multiprocessing import Pool, TimeoutError



ASDF_FILENAME = '/g/data/ha3/Passive/_ANU/7G(2013-2015)/ASDF/7G(2013-2015).h5'
OUTPUT_PATH = '/g/data/ha3/rakib/cleansedData/7G/MSEED'

def process(ds, sn_list):
    startDate = UTCDateTime("2013-01-01T00:00:00").timestamp
    endDate   = UTCDateTime("2016-01-01T00:00:00").timestamp
    step = 3600*24

    for sn in sn_list:
        logf = open(os.path.join(OUTPUT_PATH, '%s.log.txt'%(sn)), "w+")

        if(not os.path.exists(os.path.join(OUTPUT_PATH, sn))):
            os.system('mkdir %s'%os.path.join(OUTPUT_PATH, sn))

        currentTime = startDate

        trCounts = np.zeros(3)
        logf.write('Exporting mseed files for station: %s\n'%(sn))
        while(currentTime < endDate):
            zst = ds.ds_jason_db.fetchDataByTime(ds.ds, sn, ds.zchannel,
                                                 currentTime,
                                                 currentTime+step)        

            nst = ds.ds_jason_db.fetchDataByTime(ds.ds, sn, ds.nchannel,
                                                 currentTime,
                                                 currentTime+step)        

            est = ds.ds_jason_db.fetchDataByTime(ds.ds, sn, ds.echannel,
                                                 currentTime,
                                                 currentTime+step)        

            for i, st in enumerate([zst, nst, est]):
                if(st is None): continue

                fsName = '%s.%s.MSEED'%(st.traces[0].id, UTCDateTime(currentTime).strftime("%y-%m-%d.T%H:%M:%S"))

                try:
                    st.merge(method=1, fill_value=0)                
                    st.write(os.path.join(os.path.join(OUTPUT_PATH, sn), fsName), 
                             format="MSEED")
                    trCounts[i] += 1

                    #st.plot()
                except:
                    logf.write('Failed to write stream: %s\n'%(fsName))
                    logf.flush()
            # end for
            #break
            currentTime += step
        # wend
        logf.write('\t Exported (%d)Z, (%d)N and (%d)E traces.\n' % (trCounts[0], trCounts[1], trCounts[2]))
        logf.flush()
        logf.close()
    # end for
# end func

if (__name__ == '__main__'):
    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()
    proc_stations = defaultdict(list)
    ds = Dataset(ASDF_FILENAME)

    if(rank == 0):
        # split work over stations in ds1 for the time being

        count = 0
        for iproc in np.arange(nproc):
            for istation in np.arange(np.divide(len(ds.stations), nproc)):
                proc_stations[iproc].append(ds.stations[count])
                count += 1
        # end for
        for iproc in np.arange(np.mod(len(ds.stations), nproc)):
            proc_stations[iproc].append(ds.stations[count])
            count += 1
        # end for
    # end if

    # broadcast workload to all procs
    proc_stations = comm.bcast(proc_stations, root=0)
    
    print proc_stations
    process(ds, proc_stations[rank])
# end if
