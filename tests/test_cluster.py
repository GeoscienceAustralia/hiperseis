import os
import glob
import numpy as np
import pytest
import pandas as pd
from subprocess import check_call
from collections import Counter
from obspy import read_events
from obspy.geodetics import locations2degrees
from seismic.cluster.cluster import (process_event,
                                     _read_all_stations,
                                     read_stations,
                                     process_many_events,
                                     Grid)

TESTS = os.path.dirname(__file__)
PASSIVE = os.path.dirname(TESTS)
EVENTS = os.path.join(TESTS, 'mocks', 'events')
xmls = glob.glob(os.path.join(EVENTS, '*.xml'))
engdhal_xmls = glob.glob(os.path.join(EVENTS, 'engdahl_sample', '*.xml'))
stations_file = os.path.join(PASSIVE, 'inventory', 'stations.csv')
stations = _read_all_stations()
saved_out = os.path.join(TESTS, 'mocks', 'events', 'ga2017qxlpiu.csv')


@pytest.fixture(params=xmls + engdhal_xmls)
def event_xml(request):
    return request.param


@pytest.fixture(params=['P S', 'p s', 'Pn Sn', 'Pg Sg'])
def arr_type(request):
    return request.param


@pytest.fixture(params=['P S', 'Pn Sn'])
def pair_type(request):
    return request.param


@pytest.mark.filterwarnings("ignore")
def test_single_event_output(xml):
    event = read_events(xml).events[0]
    origin = event.preferred_origin()
    grid = Grid(nx=1440, ny=720, dz=25.0)
    p_arr, s_arr = process_event(read_events(xml)[0],
                                 stations=read_stations(stations_file),
                                 grid=grid,
                                 wave_type='P S')
    inputs = np.genfromtxt(saved_out, delimiter=',')
    np.testing.assert_array_almost_equal(inputs, np.array(p_arr, dtype=float),
                                         decimal=2)

    # make sure number of arrivals match that of output lines
    # no s arrivals for this event
    assert len(origin.arrivals) == len(p_arr)


@pytest.mark.filterwarnings("ignore")
def test_single_event_arrivals(event_xml, arr_type):
    event = read_events(event_xml).events[0]
    origin = event.preferred_origin()

    p_type, s_type = arr_type.split()

    grid = Grid(nx=1440, ny=720, dz=25.0)

    p_arr, s_arr = process_event(read_events(event_xml)[0],
                                 stations=stations,
                                 grid=grid,
                                 wave_type=arr_type)

    outputs_p = np.array(p_arr, dtype=float)
    outputs_s = np.array(s_arr, dtype=float)

    p_arrivals = []
    s_arrivals = []

    for arr in origin.arrivals:
        sta_code = arr.pick_id.get_referred_object(
            ).waveform_id.station_code
        if sta_code in stations:
            degrees_to_source = locations2degrees(
                origin.latitude, origin.longitude,
                float(stations[sta_code].latitude),
                float(stations[sta_code].longitude))
            if degrees_to_source < 90.0:
                if arr.phase == p_type:
                    p_arrivals.append(arr)
                if arr.phase == s_type:
                    s_arrivals.append(arr)

    # make sure number of arrivals match that of output lines
    assert len(p_arrivals) == len(outputs_p) and \
           len(s_arrivals) == len(outputs_s)

    # test that location2degress is never more than 90 degrees
    # test last columns, i.e., wave type

    if len(outputs_p) > 1:
        assert max(abs(outputs_p[:, -2])) < 90.001
        np.testing.assert_array_equal(outputs_p[:, -1], 1)
    elif len(outputs_p) == 1 and len(outputs_p[0]):
        assert abs(outputs_p[0][-2]) < 90.001
        assert outputs_p[0][-1] == 1

    if len(outputs_s) > 1:
        assert max(abs(outputs_s[:, -2])) < 90.001
        np.testing.assert_array_equal(outputs_s[:, -1], 2)
    elif len(outputs_s) == 1 and len(outputs_s[0]):
        assert abs(outputs_s[0][-2]) < 90.001
        assert outputs_s[0][-1] == 2


def test_sorted_filtered_matched(pair_type, random_filename):
    """
    check cluster sort and filter operation
    """
    outfile = random_filename()

    # gather
    gather = ['cluster', 'gather', os.path.join(EVENTS, 'engdahl_sample'),
              '-o', outfile, '-w', pair_type]
    check_call(gather)

    for wave_type in pair_type.split():  # check for both P and S
        _test_sort_and_filtered(outfile, wave_type)

    _test_matched(outfile, pair_type)


def _test_matched(outfile, wave_type):
    p, s = wave_type.split()
    sorted_p = outfile + '_sorted_' + p + '.csv'
    sorted_s = outfile + '_sorted_' + s + '.csv'

    matched_p = outfile + '_matched_' + p + '.csv'
    matched_s = outfile + '_matched_' + s + '.csv'

    match = ['cluster', 'match', sorted_p, sorted_s,
             '-p', matched_p, '-s', matched_s]
    check_call(match)

    # files created
    assert os.path.exists(matched_p) and os.path.exists(matched_s)

    pdf = pd.read_csv(matched_p, header=None)
    sdf = pd.read_csv(matched_s, header=None)
    outdf = pd.merge(pdf[[0, 1]],
                     sdf[[0, 1]],
                     how='inner',
                     on=[0, 1])

    # after inner join each df should have same number of rows
    assert outdf.shape[0] == pdf.shape[0] == sdf.shape[0]

    # make sure the arrays themselves match
    np.testing.assert_array_equal(outdf[0].values, pdf[0].values)
    np.testing.assert_array_equal(outdf[0].values, sdf[0].values)
    np.testing.assert_array_equal(outdf[1].values, pdf[1].values)
    np.testing.assert_array_equal(outdf[1].values, sdf[1].values)


def _test_sort_and_filtered(outfile, wave_type):
    sorted_p_or_s = outfile + '_sorted_' + wave_type + '.csv'
    residual = 5.0 if wave_type == 'P' or 'Pn' else 10.0
    sort_p = ['cluster', 'sort', outfile + '_' + wave_type + '.csv',
              residual, '-s', sorted_p_or_s]
    check_call(sort_p)
    assert os.path.exists(sorted_p_or_s)
    p_df = pd.read_csv(sorted_p_or_s)

    # tests for filter

    # after sorting and filtering, every group should have one row
    for _, group in p_df.groupby(by=['source_block', 'station_block']):
        # one extra due to pandas created extra index
        assert group.shape == (1, 13)

    # essentially the same thing as before
    assert len(p_df.groupby(by=['source_block', 'station_block'])) == \
        p_df.shape[0]

    # tests for sort
    # sum of source_block + station_block should be strictly increasing
    block_sum = p_df['source_block'].values + p_df['station_block'].values
    assert all(np.diff(block_sum) >= 1)


@pytest.mark.filterwarnings("ignore")
def test_multiple_event_output(random_filename):

    events = read_events(os.path.join(EVENTS, '*.xml')).events
    outfile = random_filename()

    grid = Grid(nx=1440, ny=720, dz=25.0)
    process_many_events(glob.glob(os.path.join(EVENTS, '*.xml')),
                        stations=stations,
                        grid=grid,
                        wave_type='P S',
                        output_file=outfile)

    # check files created
    assert os.path.exists(outfile + '_' + 'P_0' + '.csv')
    assert os.path.exists(outfile + '_' + 'S_0' + '.csv')

    # check all arrivals are present
    arrivals = []
    for e in events:
        origin = e.preferred_origin()
        arrivals += origin.arrivals

    p_arr = np.genfromtxt(outfile + '_' + 'P_0' + '.csv', delimiter=',')
    s_arr = np.genfromtxt(outfile + '_' + 'S_0' + '.csv', delimiter=',')

    assert len(arrivals) == len(p_arr) + len(s_arr)


@pytest.mark.filterwarning("ignore")
def test_parallel_gather(pair_type, random_filename):
    outfile_s = random_filename()
    outfile_p = random_filename()

    # gather single process
    gather_s = ['cluster', 'gather',
                os.path.join(EVENTS, 'engdahl_sample'),
                '-o', outfile_s, '-w', pair_type]
    check_call(gather_s)

    # gather multiple process
    gather_p = ['mpirun', '--allow-run-as-root', '-n', '4',
                'cluster', 'gather',
                os.path.join(EVENTS, 'engdahl_sample'),
                '-o', outfile_p, '-w', pair_type]
    check_call(gather_p)

    p, s = pair_type.split()
    assert os.path.exists(outfile_s + '_' + p + '.csv')
    assert os.path.exists(outfile_s + '_' + s + '.csv')
    assert os.path.exists(outfile_p + '_' + p + '.csv')
    assert os.path.exists(outfile_p + '_' + s + '.csv')
    sdf_p = pd.read_csv(outfile_s + '_' + p + '.csv', header=None)
    sdf_s = pd.read_csv(outfile_s + '_' + s + '.csv', header=None)
    pdf_p = pd.read_csv(outfile_p + '_' + p + '.csv', header=None)
    pdf_s = pd.read_csv(outfile_p + '_' + s + '.csv', header=None)

    # assert all columns of both df's contain the same number of elements
    for c in sdf_p.columns:
        assert Counter(sdf_p[c].values) == Counter(pdf_p[c].values)
        assert Counter(sdf_s[c].values) == Counter(pdf_s[c].values)
