#!/usr/bin/env python

import os

from seismic.analyze_station_orientations import (analyze_station_orientations,
                                                  process_event_file)
from seismic.analyze_station_orientations import (DEFAULT_CONFIG_FILTERING,
                                                  DEFAULT_CONFIG_PROCESSING)


SYNTH_CURATION_OPTS = {
    "min_snr": 2.0,
    "rms_amplitude_bounds": {"R/Z": 1.0, "T/Z": 1.0}
}

EXPECTED_SYNTH_KEY = 'SY.OAA'


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


def test_rotation_error(ned_rotation_error):
    # Test simple rotation error
    result = analyze_station_orientations(ned_rotation_error, SYNTH_CURATION_OPTS,
                                          DEFAULT_CONFIG_FILTERING,
                                          DEFAULT_CONFIG_PROCESSING)
    print(ned_rotation_error.param, result)
    assert EXPECTED_SYNTH_KEY in result
    assert result[EXPECTED_SYNTH_KEY] == -ned_rotation_error.param
# end func


def test_analyze_file(synth_event_file, tmpdir):
    # Test file interface for orientation analysis
    output_file = os.path.join(tmpdir, 'test_analyze_file.json')
    process_event_file(synth_event_file, curation_opts=SYNTH_CURATION_OPTS,
                       dest_file=output_file)
    assert os.path.isfile(output_file)
    with open(output_file, 'r') as f:
        assert f.read() == '{{\n    "{}": 0.0\n}}'.format(EXPECTED_SYNTH_KEY)
# end func
