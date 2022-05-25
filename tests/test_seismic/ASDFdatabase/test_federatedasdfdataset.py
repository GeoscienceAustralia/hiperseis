#!/bin/env python
"""
Description:
    Tests various aspects of the FederatedASDFDataSet class

References:

CreationDate:   5/20/19
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     5/20/19   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
import os
import pytest
from ordered_set import OrderedSet as set
import numpy as np
import tempfile
import sqlite3
from obspy.core import UTCDateTime
import logging

path = os.path.dirname(os.path.abspath(__file__))

# Unzip expected results
tempdir = tempfile.mkdtemp()
cmd = 'cp %s/data/asdf_test_data.h5.gz %s'%(path, tempdir)
os.system(cmd)
cmd = 'gzip -d %s'%('%s/asdf_test_data.h5.gz '%tempdir)
os.system(cmd)

# Initialize input data
asdf_file_list = os.path.join(tempdir, 'asdf_file_list.txt')

f = open(asdf_file_list, 'w+')
f.write('%s/asdf_test_data.h5'%(tempdir))
f.close()

@pytest.fixture(params=[0.0025, 0.01, 0.04])
def buffer_mb(request):
    # the first value above triggers _FederatedASDFDataSetImpl to
    # fetch data in several attempts, each time expanding the
    # data-buffer size.
    return request.param
# end func

@pytest.fixture(params=[1, 2, 4, 8, 10])
def num_neighbours(request):
    return request.param
# end func

def test_db_integrity():
    fds = FederatedASDFDataSet(asdf_file_list)

    # get number of waveforms from the db directly
    conn = sqlite3.connect(fds.fds.db_fn)
    query = 'select count(*) from wdb;'
    db_waveform_count = conn.execute(query).fetchall()[0][0]

    # fetch waveform counts for each unique combination of net, sta, loc, cha
    waveform_count = 0
    rows = fds.get_stations('1900-01-01T00:00:00', '2100-01-01T00:00:00')
    for row in rows:
        n, s, l, c, _, _, _ = row

        waveform_count += fds.get_waveform_count(n, s, l, c, '1900:01:01T00:00:00', '2100:01:01T00:00:00')
    # end for

    assert waveform_count == db_waveform_count
# end func

def test_get_stations():
    fds = FederatedASDFDataSet(asdf_file_list)

    rows = np.array(fds.get_stations('1900-01-01T00:00:00', '2100-01-01T00:00:00'))

    station_set = set()
    for n, s in rows[:, 0:2]: station_set.add((n, s))

    # There are eight stations in the h5 file
    assert len(station_set) == 8
# end func

def test_get_coordinates():
    fds = FederatedASDFDataSet(asdf_file_list)

    rows = np.array(fds.get_stations('1900-01-01T00:00:00', '2100-01-01T00:00:00'))

    station_set = set()
    for n, s in rows[:, 0:2]: station_set.add((n, s))

    # we should have coordinates for each station
    assert len(fds.unique_coordinates) == len(station_set)
# end func

def test_get_global_time_range():
    fds = FederatedASDFDataSet(asdf_file_list)

    rows = np.array(fds.get_stations('1900-01-01T00:00:00', '2100-01-01T00:00:00'))

    station_set = set()
    for n, s in rows[:, 0:2]: station_set.add((n, s))

    minlist =[]
    maxlist = []
    for (n, s) in station_set:
        min, max = fds.get_global_time_range(n, s)
        minlist.append(min)
        maxlist.append(max)
    # end for

    min = UTCDateTime(np.array(minlist).min())
    max = UTCDateTime(np.array(maxlist).max())

    # Ensure aggregate min/max to corresponding values in the db
    assert min == UTCDateTime('2000-01-01T00:00:00.000000Z')
    assert max == UTCDateTime('2002-01-01T00:00:00.000000Z')
# end func


def test_get_stations():
    fds = FederatedASDFDataSet(asdf_file_list)

    local_netsta_list = list(fds.stations_iterator())
    rows = np.array(fds.get_stations('1900-01-01T00:00:00', '2100-01-01T00:00:00'))

    # Get a list of unique stations
    stations = set()
    for n, s in rows[:,0:2]:
        stations.add((n, s))
    # end for

    # On serial runs, all stations should be allocated to rank 0
    assert len(local_netsta_list) == len(stations)
# end func

def test_get_closest_stations(num_neighbours):
    fds = FederatedASDFDataSet(asdf_file_list)

    netsta, dist = fds.get_closest_stations(0, 0, num_neighbours)

    # There are a total of 8 stations in the data set.
    assert len(netsta) > 0 and len(netsta) <= 8 and len(netsta) <= num_neighbours
# end func

def test_get_waveform(buffer_mb):
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    handler = logging.FileHandler('%s/a_%d.txt'%(tempdir, int(buffer_mb)), mode='w')
    handler.setFormatter(formatter)
    logger = logging.getLogger('test')
    logger.setLevel(logging.DEBUG)
    logger.addHandler(handler)

    fds = FederatedASDFDataSet(asdf_file_list, logger=logger,
                               single_item_read_limit_in_mb=buffer_mb)

    rows = np.array(fds.get_stations('1900-01-01T00:00:00', '2100-01-01T00:00:00'))

    for n, s, l, c in rows[:, 0:4]:
        wc = fds.get_waveform_count(n, s, l, c, '1900-01-01T00:00:00', '2100-01-01T00:00:00')
        stream = fds.get_waveforms(n, s, l, c, '1900-01-01T00:00:00', '2100-01-01T00:00:00',
                                   trace_count_threshold=1e4)

        assert wc == len(stream)
        logger.info('%s.%s: %d traces fetched'%(n, s, len(stream)))
    # end for
# end func

