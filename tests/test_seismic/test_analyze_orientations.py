#!/usr/bin/env python

import os
import copy
import pytest

import numpy as np
import matplotlib.pyplot as plt
from rf import RFStream

from seismic.analyze_station_orientations import (analyze_station_orientations,
                                                  process_event_file)
from seismic.analyze_station_orientations import (DEFAULT_CONFIG_FILTERING,
                                                  DEFAULT_CONFIG_PROCESSING)
from seismic.receiver_fn.generate_rf import transform_stream_to_rf
from seismic.stream_processing import correct_back_azimuth
from seismic.receiver_fn.rf_plot_utils import plot_rf_stack


SYNTH_CURATION_OPTS = {
    "min_snr": 2.0,
    "rms_amplitude_bounds": {"R/Z": 1.0, "T/Z": 1.0}
}

EXPECTED_SYNTH_KEY = 'SY.OAA'


def compute_ned_stacked_rf(ned):
    # Compute RF from NetworkEventDataset and return R component
    rf_all = RFStream()
    for _sta, evid, stream in ned:
        rf_3ch = transform_stream_to_rf(evid, RFStream(stream),
                                        DEFAULT_CONFIG_FILTERING,
                                        DEFAULT_CONFIG_PROCESSING)
        if rf_3ch is None:
            continue
        rf_all += rf_3ch.select(component='R')
    # end for
    rf_stacked = rf_all.stack()
    return rf_stacked, rf_all
# end func


@pytest.fixture(scope='module')
def expected_receiver_fn(_master_event_dataset, tmpdir_factory):
    ned = copy.deepcopy(_master_event_dataset)
    baseline_plot_filename = 'test_analyze_orientations_baseline.png'
    baseline_outfile = os.path.join(tmpdir_factory.mktemp('expected_receiver_fn').strpath,
                                    baseline_plot_filename)
    baseline_rf, details = compute_ned_stacked_rf(ned)
    plot_rf_stack(details, save_file=baseline_outfile)
    print('Reference RF plot saved to tmp file {}'.format(baseline_outfile))
    plt.close()
    return baseline_rf
# end func


def test_unmodified_data(ned_original):
    # Test unmodified test data which should have zero orientation error
    result = analyze_station_orientations(ned_original, SYNTH_CURATION_OPTS,
                                          DEFAULT_CONFIG_FILTERING,
                                          DEFAULT_CONFIG_PROCESSING)
    print('Original synth data', result)
    assert EXPECTED_SYNTH_KEY in result
    assert result[EXPECTED_SYNTH_KEY] == 0.0
# end func


def test_ne_swapped(ned_channel_swapped):
    # Test swapping of N and E channels
    result = analyze_station_orientations(ned_channel_swapped, SYNTH_CURATION_OPTS,
                                          DEFAULT_CONFIG_FILTERING,
                                          DEFAULT_CONFIG_PROCESSING)
    print('N-E channels swapped', result)
    assert EXPECTED_SYNTH_KEY in result
    assert result[EXPECTED_SYNTH_KEY] != 0.0
# end func


def test_channel_negated(ned_channel_negated):
    # Test negation of one or both transverse channels
    result = analyze_station_orientations(ned_channel_negated, SYNTH_CURATION_OPTS,
                                          DEFAULT_CONFIG_FILTERING,
                                          DEFAULT_CONFIG_PROCESSING)
    print(ned_channel_negated.param, result)
    assert EXPECTED_SYNTH_KEY in result
    assert result[EXPECTED_SYNTH_KEY] != 0.0
# end func


def test_rotation_error(ned_rotation_error, expected_receiver_fn, tmpdir_factory):
    # Test simple rotation error
    ned_duplicate = copy.deepcopy(ned_rotation_error)
    result = analyze_station_orientations(ned_rotation_error, SYNTH_CURATION_OPTS,
                                          DEFAULT_CONFIG_FILTERING,
                                          DEFAULT_CONFIG_PROCESSING)
    print(ned_rotation_error.param, result)
    assert EXPECTED_SYNTH_KEY in result
    assert result[EXPECTED_SYNTH_KEY] == -ned_rotation_error.param

    # Apply correction and check that we can recover proper receiver function
    ned_duplicate.apply(lambda stream: correct_back_azimuth(None, stream,
                        baz_correction=result[EXPECTED_SYNTH_KEY]))
    stacked_rf_R, rfs = compute_ned_stacked_rf(ned_duplicate)
    actual_plot_filename = 'test_analyze_orientations_actual_{}.png'.format(ned_rotation_error.param)
    actual_outfile = os.path.join(tmpdir_factory.mktemp('expected_receiver_fn').strpath,
                                  actual_plot_filename)
    plot_rf_stack(rfs, save_file=actual_outfile)
    print('Actual RF plot saved to tmp file {}'.format(actual_outfile))
    expected_data = expected_receiver_fn[0].data
    actual_data = stacked_rf_R[0].data
    assert np.allclose(actual_data, expected_data, rtol=1.0e-3, atol=1.0e-5)
# end func


def test_analyze_file(synth_event_file, tmpdir):
    # Test file interface for orientation analysis
    output_file = os.path.join(tmpdir.strpath, 'test_analyze_file.json')
    process_event_file(synth_event_file, curation_opts=SYNTH_CURATION_OPTS,
                       dest_file=output_file)
    assert os.path.isfile(output_file)
    with open(output_file, 'r') as f:
        assert f.read() == '{{\n    "{}": 0.0\n}}'.format(EXPECTED_SYNTH_KEY)
# end func
