from mpi4py import MPI
import glob, os, sys
from os.path import join, exists
from collections import defaultdict

from math import radians, cos, sin, asin, sqrt
import numpy as np
import scipy
from scipy.spatial import cKDTree

import xcorqc
from obspy import Stream, Trace, UTCDateTime
import pyasdf
import json

import click

from ASDFdatabase.seisds import SeisDB
from xcorqc import IntervalStackXCorr, xcorr2

# define utility functions
def rtp2xyz(r, theta, phi):
    xout = np.zeros((r.shape[0], 3))
    rst = r * np.sin(theta);
    xout[:, 0] = rst * np.cos(phi)
    xout[:, 1] = rst * np.sin(phi)
    xout[:, 2] = r * np.cos(theta)
    return xout
# end func

def xyz2rtp(x, y, z):
    rout = np.zeros((x.shape[0], 3))
    tmp1 = x * x + y * y
    tmp2 = tmp1 + z * z
    rout[0] = np.sqrt(tmp2)
    rout[1] = np.arctan2(sqrt(tmp1), z)
    rout[2] = np.arctan2(y, x)
    return rout
# end func

class Dataset:
    def __init__(self, asdf_file_name, station_names = '*'):

        self._data_path = asdf_file_name
        self.ds = None
        self.has_jason_db = False
        self.ds_jason_db = None
        self._earth_radius = 6371 #km

        try:
            self.ds = pyasdf.ASDFDataSet(self._data_path, mode='r')
        except:
            raise NameError('Error reading file : %s'%(self._data_path))
        # end try

        # look for json db
        files = glob.glob(os.path.dirname(self._data_path) + '/*.json')
        for f in files:
            if(os.path.splitext(os.path.basename(self._data_path))[0] in f):
                try:
                    self.ds_jason_db = SeisDB(f)
                    self.has_jason_db = True
                except:
                    raise RuntimeError('Failed to load json file:%s' % (f))
                #end try
                break
            # end if
        # end for

        # Gather station metadata
        station_subset = set(station_names.split(' ')) if station_names != '*' else station_names
        self.stations = []
        self.stations_metadata = defaultdict(list)
        zchannels = set()
        rtps = []
        for station in self.ds.waveforms:
            sn = station._station_name.split('.')[1]

            if(station_subset != '*'):
                if sn not in station_subset: continue

            self.stations.append(sn)
            self.stations_metadata[sn] = station

            rtps.append([self._earth_radius,
                         np.radians(90 - station.coordinates['latitude']),
                         np.radians(station.coordinates['longitude'])])

            for nw in station.StationXML:
                for st in nw:
                    for ch in st:
                        if 'Z' in ch.code: zchannels.add(ch.code)
        # end for

        rtps = np.array(rtps)
        xyzs = rtp2xyz(rtps[:, 0], rtps[:, 1], rtps[:, 2])

        self._tree = cKDTree(xyzs)
        self._cart_location = defaultdict(list)
        for i, s in enumerate(self.stations):
            self._cart_location[s] = xyzs[i, :]
        # end for

        assert len(zchannels)==1, 'Multiple z-channels found %s'%(zchannels)
        self.zchannel = zchannels.pop()
    # end func

    def get_closest_stations(self, station_name, other_dataset, nn=1):
        assert isinstance(station_name, str) or isinstance(station_name, unicode), 'station_name must be a string'
        assert isinstance(other_dataset, Dataset), 'other_dataset must be an instance of Dataset'
        station_name = station_name.upper()

        assert station_name in self.stations, 'station %s not found'%(station_name)

        d, l = other_dataset._tree.query(self._cart_location[station_name], nn)

        l = l[l<len(other_dataset.stations)]
        assert len(l), 'No stations found..'

        return list(np.array(other_dataset.stations)[l])
    # end func
# end class

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('data-source1',
                type=click.Path('r'))
@click.argument('data-source2',
                type=click.Path('r'))
@click.argument('output-path', required=True,
                type=click.Path(exists=True))
@click.argument('interval-seconds', required=True,
                type=int)
@click.argument('window-seconds', required=True,
                type=int)

@click.option('--ds1-dec-factor', default=1, help="Decimation factor for data-source-1")
@click.option('--ds2-dec-factor', default=1, help="Decimation factor for data-source-2")
@click.option('--nearest-neighbours', default=-1, help="Number of nearest neighbouring stations in data-source-2"
                                                      " to correlate against a given station in data-source-1. If"
                                                      " set to -1, correlations for a cross-product of all stations"
                                                      " in both data-sets are produced -- note, this is computationally"
                                                      " expensive.")
@click.option('--fmin', default=0.3, help="Lowest frequency for bandpass filter")
@click.option('--fmax', default=1., help="Highest frequency for bandpass filter")
@click.option('--station-names1', default='*', type=str,
              help="Station name(s) (space-delimited) to process in data-source-1; default is '*', which processes all available stations.")
@click.option('--station-names2', default='*', type=str,
              help="Station name(s) (space-delimited) to process in data-source-2; default is '*', which processes all available stations.")
@click.option('--start-time', default='1970-01-01T00:00:00',
              type=str,
              help="Date and time (in UTC format) to start from; default is year 1900.")
@click.option('--end-time', default='2100-01-01T00:00:00',
              type=str,
              help="Date and time (in UTC format) to stop at; default is year 2100.")
@click.option('--read-buffer-size', default=10,
              type=int,
              help="Data read buffer size; default is 10 x 'interval_seconds'. This parameter allows fetching data in bulk,"
                   " which can improve efficiency, but has no effect on the results produced")
def process(data_source1, data_source2, output_path, interval_seconds, window_seconds, ds1_dec_factor,
            ds2_dec_factor, nearest_neighbours, fmin, fmax, station_names1, station_names2, start_time,
            end_time, read_buffer_size):
    """
    DATA_SOURCE1: Path to ASDF file \n
    DATA_SOURCE2: Path to ASDF file \n
    OUTPUT_PATH: Output folder \n
    INTERVAL_SECONDS: Length of time window (s) over which to compute cross-correlations; e.g. 86400 for 1 day \n
    WINDOW_SECONDS: Length of stacking window (s); e.g 3600 for an hour. INTERVAL_SECONDS must be a multiple of
                    WINDOW_SECONDS; no stacking is performed if they are of the same size.
    """

    read_buffer_size *= interval_seconds
    station_names1 = str(station_names1)
    station_names2 = str(station_names2)

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()
    proc_stations = defaultdict(list)

    ds1 = Dataset(data_source1, station_names1)
    ds2 = Dataset(data_source2, station_names2)

    if(rank == 0):
        # split work over stations in ds1 for the time being

        count = 0
        for iproc in np.arange(nproc):
            for istation in np.arange(np.divide(len(ds1.stations), nproc)):
                proc_stations[iproc].append(ds1.stations[count])
                count += 1
        # end for
        for iproc in np.arange(np.mod(len(ds1.stations), nproc)):
            proc_stations[iproc].append(ds1.stations[count])
            count += 1
        # end for
    # end if

    # broadcast workload to all procs
    proc_stations = comm.bcast(proc_stations, root=0)

    for st1 in proc_stations[rank]:
        st2list = None
        if(nearest_neighbours != -1):
            if (data_source1==data_source2):
                st2list = set(ds1.get_closest_stations(st1, ds2, nn=nearest_neighbours+1))
                if (st1 in  st2list): st2list.remove(st1)
                st2list = list(st2list)
            else:
                st2list = ds1.get_closest_stations(st1, ds2, nn=nearest_neighbours)
        else:
            st2list = ds2.stations

        for st2 in st2list:
            startTime = UTCDateTime(start_time)
            endTime = UTCDateTime(end_time)

            x, xCorrResDict, wcResDict = IntervalStackXCorr(ds1.ds, ds1.ds_jason_db, ds2.ds, ds2.ds_jason_db, startTime,
                                                            endTime, [st1], [st2], ds1.zchannel, ds2.zchannel, ds1_dec_factor,
                                                            ds2_dec_factor, read_buffer_size, interval_seconds,
                                                            window_seconds, fmin, fmax, output_path, 2)

        # end for
    # end for
# end func

if (__name__ == '__main__'):
    process()
# end if

