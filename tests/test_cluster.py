import os
import csv
import numpy
from obspy import read_events
from seismic.cluster.cluster import process_event, _read_stations

TESTS = os.path.dirname(__file__)
stations_file = os.path.join(TESTS, 'mocks', 'inventory', 'stations.csv')
saved_out = os.path.join(TESTS, 'mocks', 'events', 'ga2017qxlpiu.csv')


def test_single_event_output(xml, random_filename):
    output_file = random_filename(ext='.csv')
    with open(stations_file, 'r') as sta_f:
        with open(output_file, 'w') as csvfile:
            writer = csv.writer(csvfile)
            process_event(read_events(xml)[0],
                          stations=_read_stations(sta_f),
                          writer=writer,
                          nx=1440, ny=720, dz=25.0)
    inputs = numpy.genfromtxt(saved_out, delimiter=',')
    outputs = numpy.genfromtxt(output_file, delimiter=',')
    numpy.testing.assert_array_almost_equal_nulp(inputs, outputs)
