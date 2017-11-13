import os
import csv
import numpy
from obspy import read_events
from seismic.cluster.cluster import process_event, _read_stations

TESTS = os.path.dirname(__file__)
xmls = os.path.join(TESTS, 'mocks', 'events', '*.xml')
stations_file = os.path.join(TESTS, 'mocks', 'inventory', 'stations.csv')
saved_out = os.path.join(TESTS, 'mocks', 'events', 'ga2017qxlpiu.csv')


def test_single_event_output(xml, random_filename):
    p_file = random_filename(ext='_p.csv')
    s_file = random_filename(ext='_s.csv')
    event = read_events(xml).events[0]
    origin = event.preferred_origin()
    with open(stations_file, 'r') as sta_f:
        with open(p_file, 'w') as p_writer:
            with open(s_file, 'w') as s_writer:
                p_writer = csv.writer(p_writer)
                process_event(read_events(xml)[0],
                              stations=_read_stations(sta_f),
                              p_writer=p_writer,
                              s_writer=s_writer,
                              nx=1440, ny=720, dz=25.0,
                              wave_type='P S')
    inputs = numpy.genfromtxt(saved_out, delimiter=',')
    outputs = numpy.genfromtxt(p_file, delimiter=',')
    numpy.testing.assert_array_almost_equal(inputs, outputs)

    # s_file is created
    assert os.path.exists(s_file)

    # make sure number of arrivals match that of output lines
    assert len(origin.arrivals) == len(outputs)


# def test_multiple_event_output(random_filename):
#     p_file = random_filename(ext='_p.csv')
#     s_file = random_filename(ext='_s.csv')