"""
Description:
    Wrapper Class for providing fast access to data contained within a set of ASDF files
References:

CreationDate:   12/12/18
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     12/12/18   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

from mpi4py import MPI
import os
import glob
import atexit
import logging
import pickle
import numpy as np

from obspy.core import Stream, UTCDateTime
from obspy import read, Trace
import pyasdf
import ujson as json
from scipy.spatial import cKDTree
from collections import defaultdict
from seismic.ASDFdatabase.utils import rtp2xyz
from FederatedASDFDataSetDBVariant import FederatedASDFDataSetDBVariant

HAS_RTREE = True
try:
    from FederatedASDFDataSetMemVariant import FederatedASDFDataSetMemVariant
except:
    HAS_RTREE = False
# end if

class FederatedASDFDataSet():
    def __init__(self, asdf_source, variant='db', use_json_db=False, logger=None):
        """
        :param asdf_source: path to a text file containing a list of ASDF files:
               Entries can be commented out with '#'
        :param variant: can be either 'mem' or 'db'
        :param use_json_db: whether to use json db if available. Has no effect if variant is 'db'
        :param logger: logger instance
        """
        self.variant = variant
        self.use_json_db = use_json_db
        self.logger = logger
        self.asdf_source = asdf_source
        self._unique_coordinates = None
        self._earth_radius = 6371 #km

        if(self.variant == 'mem'):
            if(HAS_RTREE):
                self.fds = FederatedASDFDataSetMemVariant(asdf_source, use_json_db, logger)
            else:
                raise Exception("'mem' variant cannot be used, because package rtree is missing")
            # end if
        elif(self.variant == 'db'):
            self.fds = FederatedASDFDataSetDBVariant(asdf_source, logger)
        else:
            raise Exception("Invalid variant: must be 'mem' or 'db'")
        # end if

        # Populate coordinates
        self._unique_coordinates = defaultdict(list)

        rtps_dict = defaultdict()
        for ds_dict in self.fds.asdf_station_coordinates:
            for key in ds_dict.keys():
                self._unique_coordinates[key] = [ds_dict[key][0], ds_dict[key][1]]

                rtps_dict[key] = [self._earth_radius,
                                  np.radians(90 - ds_dict[key][1]),
                                  np.radians(ds_dict[key][0])]
            # end for
        # end for

        rtps_list = []
        for k in rtps_dict.keys():
            rtps_list.append(rtps_dict[k])
        # end for
        rtps = np.array(rtps_list)
        xyzs = rtp2xyz(rtps[:, 0], rtps[:, 1], rtps[:, 2])

        self._tree = cKDTree(xyzs)
        self._key_list = np.array(rtps_dict.keys())
    # end func

    @property
    def unique_coordinates(self):
        """

        :return: dictionary containing [lon, lat] coordinates indexed by 'net.sta'
        """
        return self._unique_coordinates
    # end func

    def get_closest_stations(self, lon, lat, nn=1):
        """

        :param lon: longitude (degree)
        :param lat: latitude (degrees)
        :param nn: number of closest stations to fetch
        :return: A tuple containing a list of closest 'network.station' names and a list of distances
                 (in ascending order) in kms
        """
        xyz = rtp2xyz(np.array([self._earth_radius]),
                      np.array([np.radians(90 - lat)]),
                      np.array([np.radians(lon)]))
        d, l = self._tree.query(xyz, nn)

        if isinstance(l, int):
            l = [l]

        if (len(d.shape)==1):
            d = np.expand_dims(d, axis=0)

        l = l[l<len(self.unique_coordinates)]

        if isinstance(l, int):
            l = [l]

        return (list(self._key_list[l]), d[0, :len(l)])
    # end func

    def get_global_time_range(self, network, station, location=None, channel=None):
        """
        :param network: network code
        :param station: station code
        :param location: location code (optional)
        :param channel: channel code (optional)
        :return: tuple containing min and max times as UTCDateTime objects. If no matching records are found
                 min is set to 2100-01-01T00:00:00.000000Z and max is set to 1900-01-01T00:00:00.000000Z
        """

        return self.fds.get_global_time_range(network, station, location=location, channel=channel)
    # end func

    def get_stations(self, starttime, endtime, network=None, station=None, location=None, channel=None):
        """
        :param starttime: start time string in UTCDateTime format; can also be an instance of Obspy UTCDateTime
        :param endtime: end time string in UTCDateTime format; can also be an instance of Obspy UTCDateTime
        :param network: network code (optional)
        :param station: station code (optional)
        :param location: location code (optional)
        :param channel: channel code (optional)

        :return: a list containing [net, sta, loc, cha, lon, lat] in each row
        """
        results = self.fds.get_stations(starttime, endtime, network, station, location, channel)
        return results
    # end func

    def get_waveform_count(self, network, station, location, channel, starttime,
                           endtime, trace_count_threshold=200):
        """
        Count the number of traces within the given parameters of network, station, etc..
        and date range. This is a fast method of determing whether any trace data exists
        in a given time period, if you don't actually need the waveform data itself.

        :param network: network code
        :param station: station code
        :param location: location code
        :param channel: channel code
        :param starttime: start time string in UTCDateTime format; can also be an instance of Obspy UTCDateTime
        :param endtime: end time string in UTCDateTime format; can also be an instance of Obspy UTCDateTime
        :param trace_count_threshold: returns 0 if the number of traces within the time-rage provided
                                      exceeds the threshold (default 200). This is particularly useful for filtering
                                      out data from bad stations, e.g. those from the AU.Schools network
        :return: The number of streams containing waveform data over the time-range provided
        """
        return self.fds.get_waveform_count(network, station, location, channel,
                                           starttime, endtime, trace_count_threshold)

    def get_waveforms(self, network, station, location, channel, starttime,
                      endtime, automerge=False, trace_count_threshold=200):
        """
        :param network: network code
        :param station: station code
        :param location: location code
        :param channel: channel code
        :param starttime: start time string in UTCDateTime format; can also be an instance of Obspy UTCDateTime
        :param endtime: end time string in UTCDateTime format; can also be an instance of Obspy UTCDateTime
        :param automerge: merge traces (default False)
        :param trace_count_threshold: returns an empty Stream if the number of traces within the time-rage provided
                                      exceeds the threshold (default 200). This is particularly useful for filtering
                                      out data from bad stations, e.g. those from the AU.Schools network
        :return: an Obspy Stream containing waveform data over the time-rage provided
        """
        s = self.fds.get_waveforms(network, station, location, channel, starttime,
                              endtime, automerge, trace_count_threshold)
        return s
    # end func

    def local_net_sta_list(self):
        """
        This function provides an iterator over the entire data volume contained in all the ASDF files listed in the
        text file during instantiation. When FederatedASDFDataSet is instantiated in an MPI-parallel environment,
        meta-data for the entire data volume are equally partitioned over all processors -- in such instances, this
        function provides an iterator over the data allocated to a given processor. This functionality underpins
        parallel operations, e.g. picking arrivals.

        :return: tuples containing [net, sta, start_time, end_time]; start- and end-times are instances of Obspy
                 UTCDateTime
        """
        for item in self.fds.local_net_sta_list():
            yield item
        # end for
    # end func
# end class
