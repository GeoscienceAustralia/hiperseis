#!/usr/bin/env python

import os
import copy
import pytest

from seismic.network_event_dataset import NetworkEventDataset
from seismic.stream_processing import swap_ne_channels, negate_channel


@pytest.fixture(scope='session')
def synth_event_file():
    return os.path.join(os.path.split(__file__)[0], 'data', 'synth_event_waveforms_for_rf_test.h5')


@pytest.fixture(scope='module')
def _master_event_dataset(synth_event_file):
    master_ned = NetworkEventDataset(synth_event_file)
    return master_ned


@pytest.fixture(scope='function')
def ned_original(_master_event_dataset):
    return copy.deepcopy(_master_event_dataset)


@pytest.fixture(scope='function')
def ned_channel_swapped(_master_event_dataset):
    ned = copy.deepcopy(_master_event_dataset)
    ned.apply(lambda stream: swap_ne_channels(None, stream))
    return ned


@pytest.fixture(scope='function', params=['N', 'E', 'both'])
def ned_channel_negated(_master_event_dataset, request):
    ned = copy.deepcopy(_master_event_dataset)
    if request.param == 'both':
        ned.apply(lambda stream: negate_channel(None, stream, channel='N'))
        ned.apply(lambda stream: negate_channel(None, stream, channel='E'))
    else:
        ned.apply(lambda stream: negate_channel(None, stream, channel=request.param))
    # end if
    ned.param = request.param
    return ned

