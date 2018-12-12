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
        results = self.fds.get_stations(starttime, endtime, network, station, location, channel)
        return results
    # end func

    def get_waveforms(self, network, station, location, channel, starttime,
                      endtime, automerge=False, trace_count_threshold=200):

        s = self.fds.get_waveforms(network, station, location, channel, starttime,
                              endtime, automerge, trace_count_threshold)
        return s
    # end func

    def local_net_sta_list(self):
        for item in self.fds.local_net_sta_list():
            yield item
        # end for
    # end func
# end class
