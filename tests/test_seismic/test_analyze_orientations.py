#!/usr/bin/env python

from seismic.analyze_station_orientations import (analyze_station_orientations,
                                                  process_event_file)
from seismic.analyze_station_orientations import (DEFAULT_CURATION_OPTS,
                                                  DEFAULT_CONFIG_FILTERING,
                                                  DEFAULT_CONFIG_PROCESSING)


SYNTH_CURATION_OPTS = {
    "min_snr": 2.0,
    "rms_amplitude_bounds": {"R/Z": 1.0, "T/Z": 1.0}
}

EXPECTED_SYNTH_KEY = 'SY.OAA'


def test_unmodified_data(ned_original):
    result = analyze_station_orientations(ned_original, SYNTH_CURATION_OPTS,
                                          DEFAULT_CONFIG_FILTERING,
                                          DEFAULT_CONFIG_PROCESSING)
    print('Original synth data', result)
    assert EXPECTED_SYNTH_KEY in result
    assert result[EXPECTED_SYNTH_KEY] == 0.0
# end func


def test_ne_swapped(ned_channel_swapped):
    result = analyze_station_orientations(ned_channel_swapped, SYNTH_CURATION_OPTS,
                                          DEFAULT_CONFIG_FILTERING,
                                          DEFAULT_CONFIG_PROCESSING)
    print('N-E channels swapped', result)
    assert EXPECTED_SYNTH_KEY in result
    assert result[EXPECTED_SYNTH_KEY] != 0.0
# end func

