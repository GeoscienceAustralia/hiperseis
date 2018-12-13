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
from collections import defaultdict
from rtree import index

from FederatedASDFDataSetMemVariant import FederatedASDFDataSetMemVariant
from FederatedASDFDataSetDBVariant import FederatedASDFDataSetDBVariant

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

        if(self.variant == 'mem'):
            self.fds = FederatedASDFDataSetMemVariant(asdf_source, use_json_db, logger)
        elif(self.variant == 'db'):
            self.fds = FederatedASDFDataSetDBVariant(asdf_source, logger)
        else:
            raise Exception("Invalid variant: must be 'mem' or 'db'")
        # end if
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
