import os
from os.path import exists
import shutil
import random
import string
import glob
import numpy as np
import pytest
from pytest import approx
import pandas as pd
from subprocess import check_call
from collections import Counter
from obspy import read_events
from obspy.geodetics import locations2degrees
from seismic.cluster.cluster import (process_event,
                                     _read_all_stations,
                                     read_stations,
                                     process_many_events,
                                     Grid,
                                     column_names,
                                     Region)

TESTS = os.path.dirname(__file__)
PASSIVE = os.path.dirname(TESTS)
EVENTS = os.path.join(TESTS, 'mocks', 'events')
xmls = glob.glob(os.path.join(EVENTS, '*.xml'))
engdhal_xmls = glob.glob(os.path.join(EVENTS, 'engdahl_sample', '*.xml'))
stations_file = os.path.join(PASSIVE, 'inventory', 'stations.csv')
stations = _read_all_stations()
saved_out = os.path.join(TESTS, 'mocks', 'events', 'ga2017qxlpiu.csv')
INV_PARAM_FILE = os.path.join(PASSIVE, 'raytracer', 'params', 'param2x2')


@pytest.fixture(params=xmls + engdhal_xmls, name='event_xml')
def ev_xml(request):
    return request.param


@pytest.fixture(params=['P S', 'p s', 'Pn Sn', 'Pg Sg'], name='arr_type')
def support_arr_type(request):
    return request.param


@pytest.fixture(params=['P S', 'Pn Sn'], name='pair_type', scope='module')
def arrival_type(request):
    return request.param


@pytest.mark.filterwarnings("ignore")
def test_single_event_output(xml):
    event = read_events(xml).events[0]
    origin = event.preferred_origin()
    grid = Grid(nx=1440, ny=720, dz=25.0)
    p_arr, s_arr, miss_sta, participating_sta = process_event(
        read_events(xml)[0], stations=read_stations(stations_file),
        grid=grid, wave_type='P S')

    inputs = np.genfromtxt(saved_out, delimiter=',')
    assert inputs == approx(np.array(p_arr, dtype=float), rel=1e-2)
    # make sure number of arrivals match that of output lines
    # no s arrivals for this event
    assert len(origin.arrivals) == len(p_arr)
    assert len(miss_sta) == 0
    assert len(participating_sta) == 10


@pytest.mark.filterwarnings("ignore")
def test_single_event_arrivals(event_xml, arr_type):
    event = read_events(event_xml).events[0]
    origin = event.preferred_origin()

    p_type, s_type = arr_type.split()

    grid = Grid(nx=1440, ny=720, dz=25.0)

    p_arr, s_arr, miss_sta, arr_sta = process_event(
        read_events(event_xml)[0],
        stations=stations,
        grid=grid,
        wave_type=arr_type)

    outputs_p = np.array(p_arr, dtype=float)
    outputs_s = np.array(s_arr, dtype=float)

    p_arrivals = []
    s_arrivals = []

    station_not_found = []
    arrival_stations = []

    for arr in origin.arrivals:
        sta_code = arr.pick_id.get_referred_object(
            ).waveform_id.station_code

        if sta_code in stations:
            arrival_stations.append(sta_code)
            degrees_to_source = locations2degrees(
                origin.latitude, origin.longitude,
                float(stations[sta_code].latitude),
                float(stations[sta_code].longitude))
            if degrees_to_source < 90.0:
                if arr.phase == p_type:
                    p_arrivals.append(arr)
                if arr.phase == s_type:
                    s_arrivals.append(arr)
        else:
            station_not_found.append(sta_code)

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

    assert len(station_not_found) <= 1
    if len(station_not_found):
        # only stations not in dict of the test events is 'ISG
        assert 'ISG' in station_not_found
    assert set(arrival_stations) == set(arr_sta)


# a very large residual allowed, imply we are not really using the filter
@pytest.fixture(params=[1e6, 1], name='residual_bool')
def res_bool(request):
    if request.param == 1e6:
        request.applymarker(pytest.mark.xfail)
    return request.param


@pytest.fixture(scope='module')
def cluster_outfiles(request, tmpdir_factory):
    print('In Setup')
    dir = str(tmpdir_factory.mktemp('seismic').realpath())
    outfiles = []
    pair_types = ['P S', 'Pn Sn']

    # gather
    for p in pair_types:
        fname = ''.join(random.choice(string.ascii_lowercase)
                        for _ in range(10))
        outfile_p = os.path.join(dir, fname)
        gather = [
                'mpirun', '--allow-run-as-root', '-n', '4',
                'cluster', 'gather', os.path.join(EVENTS, 'engdahl_sample'),
                '-o', outfile_p, '-w', p
            ]
        outfiles.append(outfile_p)
        pp, ss = p.split()
        check_call(gather)
        p_file = outfile_p + '_' + pp + '.csv'
        s_file = outfile_p + '_' + ss + '.csv'
        miss_st_file = outfile_p + '_missing_stations.csv'
        st_file = outfile_p + '_participating_stations.csv'

        assert exists(p_file) and exists(s_file) and exists(miss_st_file) and\
            exists(st_file)

    def tear_down():
        print('In tear down, removing dir: ', dir)
        shutil.rmtree(dir)
    request.addfinalizer(tear_down)
    return outfiles


def test_sorted_filtered_matched_zoned(residual_bool, pair_type,
                                       cluster_outfiles):
    """
    check cluster sort and filter operation
    """
    pair_types = ['P S', 'Pn Sn']
    outfile = cluster_outfiles[pair_types.index(pair_type)]
    for wave_type in pair_type.split():  # check for both P and S
        _test_sort_and_filtered(outfile, wave_type, residual_bool)

    _test_matched(outfile, pair_type)

    _test_zones(outfile, pair_type)


def _test_zones(outfile, pair_type):
    p, s = pair_type.split()

    matched_p = outfile + '_matched_' + p + '.csv'
    matched_s = outfile + '_matched_' + s + '.csv'

    region = '0 -50.0 100 160'

    region_p = outfile + 'region_{}.csv'.format(p)
    global_p = outfile + 'global_{}.csv'.format(p)
    region_s = outfile + 'region_{}.csv'.format(s)
    global_s = outfile + 'global_{}.csv'.format(s)

    zone_p = ['cluster', 'zone', '-z', region, matched_p,
              '-r', region_p,
              '-g', global_p]
    zone_s = ['cluster', 'zone', '-z', region, matched_s,
              '-r', region_s,
              '-g', global_s]

    check_call(zone_p)
    check_call(zone_s)

    for w in pair_type.split():
        _test_zone(outfile, region, w)


def _add_dicts(x, y):
    return {k: x.get(k, 0) + y.get(k, 0) for k in set(x) | set(y)}


def _test_zone(outfile, region, wave_type):
    region_p = outfile + 'region_{}.csv'.format(wave_type)
    global_p = outfile + 'global_{}.csv'.format(wave_type)
    matched_p = outfile + '_matched_' + wave_type + '.csv'

    prdf = pd.read_csv(region_p, header=None, names=column_names, sep=' ')
    pgdf = pd.read_csv(global_p, header=None, names=column_names, sep=' ')
    # ensure there are no overlaps between the regional and global df's
    null_df = prdf.merge(pgdf, how='inner', on=['source_block',
                                                'station_block'])
    assert null_df.shape[0] == 0

    # reconstruct the original matched files combining prdf an pgdf
    m_pdf = pd.read_csv(matched_p, header=None, names=column_names, sep=' ')

    # shapes match
    assert m_pdf.shape[0] == prdf.shape[0] + pgdf.shape[0]

    # assert elements match
    prdf_m_pdf = m_pdf.merge(prdf, how='inner',
                             on=['source_block', 'station_block'])
    assert prdf_m_pdf.shape[0] == prdf.shape[0]
    for c in ['source_block', 'station_block', 'event_number', 'P_or_S']:
        set(prdf[c].values).issubset(set(m_pdf[c].values))
        set(pgdf[c].values).issubset(set(m_pdf[c].values))
        c_prdf = Counter(prdf[c].values)
        c_pgdf = Counter(pgdf[c].values)
        c_mdf = Counter(m_pdf[c].values)

        assert c_mdf == _add_dicts(c_pgdf, c_prdf)

    region = [float(s) for s in region.split()]
    region = Region(*region)

    # check region bounds are satisfied
    assert all(np.logical_or(
        (prdf['source_latitude'].values < region.upperlat),
        (prdf['station_latitude'].values < region.upperlat)))

    assert all(np.logical_or(
        (prdf['source_latitude'].values > region.bottomlat),
        (prdf['station_latitude'].values > region.bottomlat)))

    assert all(np.logical_or(
        (prdf['source_longitude'].values < region.rightlon),
        (prdf['station_longitude'].values < region.rightlon)))

    assert all(np.logical_or(
        (prdf['source_longitude'].values > region.leftlon),
        (prdf['station_longitude'].values > region.leftlon)))

    _check_range(prdf)
    _check_range(pgdf)


def _test_output_stations_check(df):
    """
    This test checks that the stations in the output files are the same as
    that of the arrivals recorded in the event xmls and can be found in the
    `stations` dict.

    Ref: https://github.com/GeoscienceAustralia/passive-seismic/issues/61

    :param df: dataframe which will be checked for stations at any stage of
    the processing
    """



def _check_range(prdf):

    # check co-longitude range
    assert all(prdf['source_longitude'].values >= 0)
    assert all(prdf['station_longitude'].values >= 0)
    assert all(prdf['source_longitude'].values <= 360)
    assert all(prdf['station_longitude'].values <= 360)


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

    pdf = pd.read_csv(matched_p, header=None, sep=' ')
    sdf = pd.read_csv(matched_s, header=None, sep=' ')
    outdf = pd.merge(pdf[[0, 1]],
                     sdf[[0, 1]],
                     how='inner',
                     on=[0, 1])

    # after inner join each df should have same number of rows
    assert outdf.shape[0] == pdf.shape[0] == sdf.shape[0]

    # make sure the block numbers for both sources and stations match
    np.testing.assert_array_equal(outdf[0].values, pdf[0].values)
    np.testing.assert_array_equal(outdf[0].values, sdf[0].values)
    np.testing.assert_array_equal(outdf[1].values, pdf[1].values)
    np.testing.assert_array_equal(outdf[1].values, sdf[1].values)


def _test_sort_and_filtered(outfile, wave_type, residual_bool):
    sorted_p_or_s = outfile + '_sorted_' + wave_type + '.csv'
    p_s_res = 5.0 if wave_type == 'P' or 'Pn' else 10.0
    residual = residual_bool*p_s_res
    sort_p = ['cluster', 'sort', outfile + '_' + wave_type + '.csv',
              str(residual), '-s', sorted_p_or_s]
    check_call(sort_p)
    assert os.path.exists(sorted_p_or_s)
    p_df = pd.read_csv(sorted_p_or_s, header=None, names=column_names, sep=' ')

    # tests for residual filter
    assert all(abs(p_df['residual'].values) < p_s_res)

    # tests for median filter
    # after sorting and filtering, every group should have one row
    for _, group in p_df.groupby(by=['source_block', 'station_block']):
        assert group.shape == (1, 12)

    # essentially the same thing as before
    assert len(p_df.groupby(by=['source_block', 'station_block'])) == \
        p_df.shape[0]

    # tests for sort
    # source_block should be increasing
    assert all(np.diff(p_df['source_block'].values) >= 0)


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
    assert os.path.exists(outfile + '_P_0' + '.csv')
    assert os.path.exists(outfile + '_S_0' + '.csv')
    assert os.path.exists(outfile + '_missing_stations_0' + '.csv')
    assert os.path.exists(outfile + '_participating_stations_0' + '.csv')

    # check all arrivals are present
    arrivals = []
    arriving_stations = []
    for e in events:
        origin = e.preferred_origin()
        arrivals += origin.arrivals
        for a in arrivals:
            arriving_stations.append(a.pick_id.get_referred_object(
                ).waveform_id.station_code)

    p_arr = np.genfromtxt(outfile + '_' + 'P_0' + '.csv', delimiter=',')
    s_arr = np.genfromtxt(outfile + '_' + 'S_0' + '.csv', delimiter=',')

    assert len(arrivals) == len(p_arr) + len(s_arr)

    assert len(arrivals) == len(arriving_stations)

    output_arr_stations = set(pd.read_csv(
        outfile + '_participating_stations_0' + '.csv', header=None)[0])

    assert set(output_arr_stations) == set(arriving_stations)


@pytest.mark.filterwarning("ignore")
def test_parallel_gather(pair_type, random_filename):
    outfile_s = random_filename()
    outfile_p = random_filename()

    # gather single process
    gather_s = ['cluster', 'gather',
                os.path.join(EVENTS, 'engdahl_sample'),
                '-o', outfile_s, '-w', pair_type]
    check_call(gather_s)

    # gather multip`le process
    gather_p = ['mpirun', '--allow-run-as-root', '-n', '4',
                'cluster', 'gather',
                os.path.join(EVENTS, 'engdahl_sample'),
                '-o', outfile_p, '-w', pair_type]
    check_call(gather_p)

    p, s = pair_type.split()

    sdf_p = pd.read_csv(outfile_s + '_' + p + '.csv', header=None)
    sdf_s = pd.read_csv(outfile_s + '_' + s + '.csv', header=None)
    pdf_p = pd.read_csv(outfile_p + '_' + p + '.csv', header=None)
    pdf_s = pd.read_csv(outfile_p + '_' + s + '.csv', header=None)

    # assert all columns of both df's contain the same number of elements
    for c in sdf_p.columns:
        assert Counter(sdf_p[c].values) == Counter(pdf_p[c].values)
        assert Counter(sdf_s[c].values) == Counter(pdf_s[c].values)
