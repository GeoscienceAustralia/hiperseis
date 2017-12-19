from __future__ import print_function, absolute_import
import os
from os.path import exists
import shutil
import random
import string
import glob
import csv
import numpy as np
import pytest
from pytest import approx
import pandas as pd
from subprocess import check_call
from collections import Counter
from obspy import read_events
from obspy.core.event import Catalog
from obspy.geodetics import locations2degrees
from seismic.mpiops import rank
from seismic.cluster.cluster import (process_event,
                                     _read_all_stations,
                                     read_stations,
                                     process_many_events,
                                     Grid,
                                     column_names,
                                     Region,
                                     recursive_glob,
                                     STATION_LATITUDE,
                                     STATION_LONGITUDE,
                                     STATION_CODE,
                                     FREQUENCY,
                                     _in_region)

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

    # clear the station_code
    p_arr = [p[:11] + [p[12]] for p in p_arr]

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

    # clear the station_code
    p_arr = [p[:11] + [p[12]] for p in p_arr]
    s_arr = [s[:11] + [s[12]] for s in s_arr]

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
            degrees_to_source = locations2degrees(
                    origin.latitude, origin.longitude,
                    stations[sta_code].latitude,
                    stations[sta_code].longitude)
            if degrees_to_source < 90.0:
                if arr.phase == p_type:
                    p_arrivals.append(arr)
                if arr.phase == s_type:
                    s_arrivals.append(arr)
                if arr.phase in [p_type, s_type]:
                    arrival_stations.append(sta_code)
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


@pytest.fixture(params=[True, False], name='reject_stations')
def rej_stas(request):
    return request.param


@pytest.fixture(name='reject_stations_file')
def rej_sta_f(request, random_filename):
    reject_stations_file = random_filename(ext='.csv')
    stas = ['LRMC', 'MAT', 'MBWA', 'MDJ', 'MIR', 'MOO', 'MUN', 'NWAO']

    with open(reject_stations_file, 'w') as csvfile:
        writer = csv.writer(csvfile)
        for s in stas:
            writer.writerow([s])
    return reject_stations_file


def test_sorted_filtered_matched_zoned(residual_bool, pair_type,
                                       cluster_outfiles, reject_stations,
                                       reject_stations_file):
    """
    check cluster sort and filter operation
    """
    pair_types = ['P S', 'Pn Sn']
    outfile = cluster_outfiles[pair_types.index(pair_type)]
    for wave_type in pair_type.split():  # check for both P and S
        _test_sort_and_filtered(outfile, wave_type, residual_bool)

    _test_matched(outfile, pair_type)

    for w in pair_type.split():
        _test_zone(outfile, '0 -50.0 100 160', w, reject_stations,
                   reject_stations_file)


def _add_dicts(x, y):
    return {k: x.get(k, 0) + y.get(k, 0) for k in set(x) | set(y)}


def _test_zone(outfile, region, wave_type, reject_stations,
               reject_stations_file):

    region_p = outfile + 'region_{}.csv'.format(wave_type)
    global_p = outfile + 'global_{}.csv'.format(wave_type)
    matched_p = outfile + '_matched_' + wave_type + '.csv'

    zone_p = ['cluster', 'zone', '-z', region, matched_p,
              '-r', region_p,
              '-g', global_p]

    if reject_stations:
        zone_p += ['-j', reject_stations_file]

    check_call(zone_p)

    col_names = list(column_names)
    col_names.remove(STATION_CODE)

    prdf = pd.read_csv(region_p, header=None, names=col_names, sep=' ')
    pgdf = pd.read_csv(global_p, header=None, names=col_names, sep=' ')

    # ensure there are no overlaps between the regional and global df's
    null_df = prdf.merge(pgdf, how='inner', on=['source_block',
                                                'station_block'])
    assert null_df.shape[0] == 0

    # reconstruct the original matched files combining prdf an pgdf
    m_pdf = pd.read_csv(matched_p, header=None, names=column_names, sep=' ')

    if reject_stations:
        reg = Region(*[float(s) for s in region.split()])
        df_region, global_df, x_region_df = _in_region(reg,
                                                       m_pdf,
                                                       grid_size=0.0)
        rej_sta = pd.read_csv(reject_stations_file, header=None,
                              names=[STATION_CODE])
        r_region_df = df_region.merge(rej_sta, how='inner',
                                      on=STATION_CODE)
        r_global_df = global_df.merge(rej_sta, how='inner',
                                      on=STATION_CODE)
    else:
        r_region_df = pd.DataFrame()
        r_global_df = pd.DataFrame()

    # shapes match
    assert m_pdf.shape[0] == prdf.shape[0] + pgdf.shape[0] + \
                             r_global_df.shape[0] + r_region_df.shape[0]

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
        if reject_stations:
            c_r_prdf = Counter(r_region_df[c].values)
            c_r_pgdf = Counter(r_global_df[c].values)
        else:
            c_r_pgdf = Counter()
            c_r_prdf = Counter()

        assert c_mdf == _add_dicts(_add_dicts(_add_dicts(c_pgdf, c_prdf),
                                              c_r_prdf), c_r_pgdf)

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
    _test_output_stations_check(prdf, outfile)
    _test_output_stations_check(pgdf, outfile)

    # assert stats files created
    for f in [matched_p, region_p, global_p]:
        stats_file = os.path.splitext(f)[0] + '_stats.csv'
        assert exists(stats_file)

    # read arrivals
    reg = pd.read_csv(os.path.splitext(region_p)[0] + '_stats.csv')
    glo = pd.read_csv(os.path.splitext(global_p)[0] + '_stats.csv')
    orig = pd.read_csv(os.path.splitext(matched_p)[0] + '_stats.csv')

    # make sure frequency of the original match that of global and regional
    reg_d = {s: n for s, n in zip(reg[STATION_CODE], reg[FREQUENCY])}
    glo_d = {s: n for s, n in zip(glo[STATION_CODE], glo[FREQUENCY])}
    orig_d = {s: n for s, n in zip(orig[STATION_CODE], orig[FREQUENCY])}

    if reject_stations:
        for k in rej_sta[STATION_CODE]:
            orig_d.pop(k, None)

    assert _add_dicts(reg_d, glo_d) == orig_d


def _test_output_stations_check(df, outfile, co_longitude=True):
    """
    This test checks that the stations in the output files are the same as
    that of the arrivals recorded in the event xmls and can be found in the
    `stations` dict.

    Ref: https://github.com/GeoscienceAustralia/passive-seismic/issues/61

    :param df: dataframe which will be checked for stations at any stage of
    the processing
    :param outfile: base output filename
    """
    lat_lon = [(lat, lon) for lat, lon
               in zip(df[STATION_LATITUDE], df[STATION_LONGITUDE])]

    st_df = pd.DataFrame.from_dict(stations, orient='index')

    gathered_stations_file =\
        outfile + '_participating_stations_{}'.format(rank) + '.csv'

    if not os.path.exists(gathered_stations_file):
        gathered_stations_file = outfile + '_participating_stations.csv'

    arrival_station_codes = set(pd.read_csv(
        gathered_stations_file, header=None)[0])

    compared = 0
    for lat, lon in lat_lon:
        compared += 1
        # find row index of matching lat/lon
        row1 = set(st_df.index[abs(st_df['latitude'] - lat) < 1e-3].tolist())
        longs = st_df['longitude'] % 360 if co_longitude else \
            st_df['longitude']
        row2 = set(st_df.index[abs(longs - lon) < 1e-3].tolist())
        # make sure that the rows match for at least one station
        sta = row1.intersection(row2)
        assert len(sta) >= 1  # at least one match
        # make sure station code can be found in arrivals stations list output
        #  during gather operation
        assert any([s in arrival_station_codes for s in sta])

    # make sure  all lat/lon/entire df was checked
    assert compared == len(lat_lon) == df.shape[0]


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
    assert exists(matched_p) and exists(matched_s)

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
    assert exists(sorted_p_or_s)
    p_df = pd.read_csv(sorted_p_or_s, header=None, names=column_names, sep=' ')

    # tests for residual filter
    assert all(abs(p_df['residual'].values) < p_s_res)

    # tests for median filter
    # after sorting and filtering, every group should have one row
    for _, group in p_df.groupby(by=['source_block', 'station_block']):
        assert group.shape == (1, 13)

    # essentially the same thing as before
    assert len(p_df.groupby(by=['source_block', 'station_block'])) == \
        p_df.shape[0]

    # tests for sort
    # source_block should be increasing
    assert all(np.diff(p_df['source_block'].values) >= 0)


def recursive_read_events(xmls):
    cat = Catalog()
    for x in xmls:
        cat.events += read_events(x).events
    return cat.events


@pytest.mark.filterwarnings("ignore")
def test_multiple_events_gather(random_filename):
    xmls = recursive_glob(EVENTS, '*.xml')
    events = recursive_read_events(xmls)
    outfile = random_filename()
    grid = Grid(nx=1440, ny=720, dz=25.0)
    process_many_events(xmls,
                        stations=stations,
                        grid=grid,
                        wave_type='P S',
                        output_file=outfile)

    # check files created
    assert exists(outfile + '_P_{}'.format(rank) + '.csv')
    assert exists(outfile + '_S_{}'.format(rank) + '.csv')
    assert exists(outfile + '_missing_stations_{}'.format(rank) + '.csv')
    assert exists(outfile + '_participating_stations_{}'.format(rank) + '.csv')

    # check all arrivals are present
    wave_types = 'P S'.split()
    arrivals = []
    arriving_stations = []
    for e in events:
        origin = e.preferred_origin() or e.origins[0]
        this_event_arrivals = origin.arrivals
        if origin.latitude is None or origin.longitude is None or \
                        origin.depth is None:
            continue

        for a in this_event_arrivals:
            sta_code = a.pick_id.get_referred_object().waveform_id.station_code
            if sta_code not in stations:
                continue
            sta = stations[sta_code]
            degrees_to_source = locations2degrees(origin.latitude,
                                                  origin.longitude,
                                                  float(sta.latitude),
                                                  float(sta.longitude))

            # ignore stations more than 90 degrees from source
            if degrees_to_source > 90.0:
                continue

            if a.phase in wave_types:
                arrivals.append(a)
                arriving_stations.append(sta_code)

    p_arr = np.genfromtxt(outfile + '_' + 'P_{}'.format(rank) + '.csv',
                          delimiter=',')
    s_arr = np.genfromtxt(outfile + '_' + 'S_{}'.format(rank) + '.csv',
                          delimiter=',')

    assert len(arrivals) == len(p_arr) + len(s_arr) == len(arriving_stations)

    output_arr_stations = set(pd.read_csv(
        outfile + '_participating_stations_{}'.format(rank) + '.csv',
        header=None)[0])

    assert set(output_arr_stations) == set(arriving_stations)
    df = pd.DataFrame(p_arr)
    df.columns = column_names
    _test_output_stations_check(df, outfile, co_longitude=False)


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
