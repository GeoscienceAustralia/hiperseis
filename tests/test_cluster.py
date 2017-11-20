import os
import csv
import glob
import numpy as np
import pytest
import pandas as pd
from subprocess import check_call
from obspy import read_events
from obspy.geodetics import locations2degrees
from seismic.cluster.cluster import (process_event,
                                     read_stations,
                                     process_many_events)

TESTS = os.path.dirname(__file__)
PASSIVE = os.path.dirname(TESTS)
EVENTS = os.path.join(TESTS, 'mocks', 'events')
xmls = glob.glob(os.path.join(EVENTS, '*.xml'))
engdhal_xmls = glob.glob(os.path.join(EVENTS, 'engdahl_sample', '*.xml'))
stations_file = os.path.join(PASSIVE, 'inventory', 'stations.csv')
saved_out = os.path.join(TESTS, 'mocks', 'events', 'ga2017qxlpiu.csv')


@pytest.fixture(params=xmls + engdhal_xmls)
def event_xml(request):
    return request.param


@pytest.fixture(params=['P S', 'p s', 'Pn Sn', 'Pg Sg'])
def arr_type(request):
    return request.param


@pytest.mark.filterwarnings("ignore")
def test_single_event_output(xml, random_filename):
    p_file = random_filename(ext='_p.csv')
    s_file = random_filename(ext='_s.csv')
    event = read_events(xml).events[0]
    origin = event.preferred_origin()

    with open(p_file, 'w') as p_wrtr:
        with open(s_file, 'w') as s_wrtr:
            p_writer = csv.writer(p_wrtr)
            s_writer = csv.writer(s_wrtr)
            process_event(read_events(xml)[0],
                          stations=read_stations(stations_file),
                          p_writer=p_writer,
                          s_writer=s_writer,
                          nx=1440, ny=720, dz=25.0,
                          wave_type='P S')

    inputs = np.genfromtxt(saved_out, delimiter=',')
    outputs = np.genfromtxt(p_file, delimiter=',')

    np.testing.assert_array_almost_equal(inputs, outputs)

    # s_file is created
    assert os.path.exists(s_file)

    # make sure number of arrivals match that of output lines
    # no s arrivals for this event
    assert len(origin.arrivals) == outputs.shape[0]


@pytest.mark.filterwarnings("ignore")
def test_single_event_arrivals(event_xml, random_filename, arr_type):
    p_file = random_filename(ext='_p.csv')
    s_file = random_filename(ext='_s.csv')
    event = read_events(event_xml).events[0]
    origin = event.preferred_origin()

    with open(p_file, 'w') as p_wrt:
        with open(s_file, 'w') as s_wrt:
            p_writer = csv.writer(p_wrt)
            s_writer = csv.writer(s_wrt)
            process_event(read_events(event_xml)[0],
                          stations=read_stations(stations_file),
                          p_writer=p_writer,
                          s_writer=s_writer,
                          nx=1440, ny=720, dz=25.0,
                          wave_type=arr_type)

    outputs_p = np.genfromtxt(p_file, delimiter=',')
    outputs_s = np.genfromtxt(s_file, delimiter=',')

    stations = read_stations(stations_file)

    p_arrivals = []
    s_arrivals = []
    p, s = arr_type.split()
    for arr in origin.arrivals:
        sta_code = arr.pick_id.get_referred_object(
            ).waveform_id.station_code
        if sta_code in stations:
            degrees_to_source = locations2degrees(
                origin.latitude, origin.longitude,
                float(stations[sta_code].latitude),
                float(stations[sta_code].longitude))
            if degrees_to_source < 90.0:
                if arr.phase == p:
                    p_arrivals.append(arr)
                if arr.phase == s:
                    s_arrivals.append(arr)

    if len(outputs_p.shape) == 1 and outputs_p.shape[0]:
        out_shape_p = 1
    else:
        out_shape_p = outputs_p.shape[0]

    if len(outputs_s.shape) == 1 and outputs_s.shape[0]:
        out_shape_s = 1
    else:
        out_shape_s = outputs_s.shape[0]

    # make sure number of arrivals match that of output lines
    assert len(p_arrivals) == out_shape_p
    assert len(s_arrivals) == out_shape_s

    # test that location2degress is never more than 90 degrees
    # test last columns, i.e., wave type

    if out_shape_p > 1:
        assert max(abs(outputs_p[:, -2])) < 90.001
        np.testing.assert_array_equal(outputs_p[:, -1], 1)
    elif out_shape_p == 1:
        assert abs(outputs_p[-2]) < 90.001
        assert outputs_p[-1] == 1

    if out_shape_s > 1:
        assert max(abs(outputs_s[:, -2])) < 90.001
        np.testing.assert_array_equal(outputs_s[:, -1], 2)
    elif out_shape_s == 1:
        assert abs(outputs_s[-2]) < 90.001
        assert outputs_s[-1] == 2


def test_sorted(random_filename):
    """
    check cluster sort and filter operation
    """
    outfile = random_filename()

    # gather
    gather = ['cluster', 'gather', os.path.join(EVENTS, 'engdahl_sample'),
              '-o', outfile]
    check_call(gather)

    # sort
    sorted_p = outfile + '_sorted_' + 'P' + '.csv'
    sorted_s = outfile + '_sorted_' + 'S' + '.csv'
    sort_p = ['cluster', 'sort', outfile + '_' + 'P' + '.csv',
              '-s', sorted_p]
    sort_s = ['cluster', 'sort', outfile + '_' + 'S' + '.csv',
              '-s', sorted_s]

    check_call(sort_p)
    check_call(sort_s)
    assert os.path.exists(outfile + '_sorted_' + 'P' + '.csv')
    assert os.path.exists(outfile + '_sorted_' + 'S' + '.csv')
    p_df = pd.read_csv(sorted_p)
    s_df = pd.read_csv(sorted_s)


def test_filtered():
    pass


@pytest.mark.filterwarnings("ignore")
def test_multiple_event_output(random_filename):

    events = read_events(os.path.join(EVENTS, '*.xml')).events
    outfile = random_filename()

    with open(stations_file, 'r') as sta_f:
        process_many_events(events,
                            stations=read_stations(sta_f),
                            nx=1440, ny=720, dz=25.0,
                            wave_type='P S',
                            output_file=outfile)

    # check files created
    assert os.path.exists(outfile + '_' + 'P' + '.csv')
    assert os.path.exists(outfile + '_' + 'S' + '.csv')

    # check all arrivals are present
    arrivals = []
    for e in events:
        origin = e.preferred_origin()
        arrivals += origin.arrivals

    p_arr = np.genfromtxt(outfile + '_' + 'P' + '.csv', delimiter=',')
    s_arr = np.genfromtxt(outfile + '_' + 'S' + '.csv', delimiter=',')

    assert len(arrivals) == len(p_arr) + len(s_arr)


def test_matched_files():
    pass
