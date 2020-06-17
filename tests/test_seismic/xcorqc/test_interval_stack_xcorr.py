#!/bin/env python
"""
Description:
    Tests various aspects of x-correlation functionality

References:

CreationDate:   5/20/19
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     5/20/19   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
from seismic.xcorqc.xcorqc import IntervalStackXCorr
from obspy import read_inventory
import os, sys
import pytest
from netCDF4 import Dataset
import numpy as np
import tempfile
import shutil

# Prepare input
netsta1 = 'AU.ARMA'
netsta2 = 'AU.QLP'

# TODO: Fix resource management here so that asdf_files_dir gets deleted when tests finished/finalized.
path = os.path.dirname(os.path.abspath(__file__))

# Initialize input data
asdf_files_dir = tempfile.mkdtemp(suffix='_test')
asdf_file_list1 = os.path.join(asdf_files_dir, 'asdf_file_list1.txt')
asdf_file_list2 = os.path.join(asdf_files_dir, 'asdf_file_list2.txt')

f1 = open(asdf_file_list1, 'w+')
f2 = open(asdf_file_list2, 'w+')
f1.write('%s/data/test_data_ARMA.h5\n'%(path))
f2.write('%s/data/test_data_QLP.h5\n'%(path))
f1.close()
f2.close()

fds1 = FederatedASDFDataSet(asdf_file_list1)
fds2 = FederatedASDFDataSet(asdf_file_list2)

# Initialize input inventory
inv = None
inv = read_inventory('%s/data/response_inventory.fdsnxml'%(path))

# Unzip expected results
expected_folder = tempfile.mkdtemp()
cmd = 'tar -zxvf %s -C %s'%('%s/data/expected/expected.tar.gz'%path, expected_folder)
os.system(cmd)

@pytest.fixture(params=['BHZ', '00T'])
def cha(request):
    return request.param

@pytest.fixture(params=[86400, 3600])
def interval_seconds(request):
    return request.param

@pytest.fixture(params=[3600, 1800])
def window_seconds(request):
    return request.param

@pytest.fixture(params=[False, True])
def envelope_normalize(request):
    return request.param

@pytest.fixture(params=[False, True])
def ensemble_stack(request):
    return request.param

@pytest.fixture(params=[False, True])
def clip_to_2std(request):
    return request.param

@pytest.fixture(params=[False, True])
def one_bit_normalize(request):
    return request.param

@pytest.fixture(params=[False, True])
def whitening(request):
    return request.param

@pytest.fixture(params=[inv.select(network='AU', station='ARMA'), None])
def inv1(request):
    return request.param

@pytest.fixture(params=[inv.select(network='AU', station='QLP'), None])
def inv2(request):
    return request.param

def test_interval_stack_xcorr_(tmpdir, cha, inv1, inv2, interval_seconds, window_seconds,
                               clip_to_2std, whitening, one_bit_normalize,
                               envelope_normalize, ensemble_stack):
    start_time = '2011-03-11T00:00:00Z'
    end_time   = '2011-03-14T00:00:00Z'

    tag = '%s.%d.%d.cs%d.w%d.obn%d.en%d.es%d.resp%d.resp%d'%(cha, interval_seconds, window_seconds,
                                                          clip_to_2std, whitening, one_bit_normalize,
                                                          envelope_normalize, ensemble_stack,
                                                          1 if inv1 else 0,
                                                          1 if inv2 else 0)

    # skipping inconsistent parameterizations
    if (one_bit_normalize and clip_to_2std): return

    if isinstance(tmpdir, str):
        output_folder = os.path.join(tmpdir, 'output')
        os.makedirs(output_folder)
    else:
        output_folder = str(tmpdir.mkdir('output'))
    # end if

    IntervalStackXCorr(fds1, fds2,
                       start_time, end_time,
                       netsta1, netsta2,
                       inv1,
                       inv2,
                       'vel',
                       50,
                       cha,
                       cha,
                       50, 250,
                       resample_rate=0.25,
                       buffer_seconds=interval_seconds,
                       interval_seconds=interval_seconds,
                       window_seconds=window_seconds, flo=0.01, fhi=0.125,
                       clip_to_2std=clip_to_2std, whitening=whitening,
                       one_bit_normalize=one_bit_normalize, envelope_normalize=envelope_normalize,
                       ensemble_stack=ensemble_stack,
                       outputPath=output_folder, verbose=2, tracking_tag=tag)

    # Read result
    fn = os.path.join(output_folder, '%s.%s.%s.nc'%(netsta1, netsta2, tag))
    dc = Dataset(fn)
    xcorr_c = dc.variables['xcorr'][:]

    # Read expected
    fn = '%s/%s.%s.%s.nc'%(expected_folder, netsta1, netsta2, tag)
    de = Dataset(fn)
    xcorr_e = de.variables['xcorr'][:]

    rtol = 1e-2
    atol = 1e-2
    if(clip_to_2std or whitening or one_bit_normalize):
        rtol *= 50
        atol *= 50
    # end if

    assert np.allclose(xcorr_c, xcorr_e, rtol=rtol, atol=atol)

    shutil.rmtree(output_folder)
# end func


if __name__ == '__main__':
    test_dir = tempfile.mkdtemp()
    cha = 'BHZ'
    inv1 = inv.select(network='AU', station='ARMA')
    inv2 = inv.select(network='AU', station='QLP')
    interval_seconds = 86400
    window_seconds = 3600
    clip_to_2std = False
    whitening = False
    one_bit_normalize = False
    envelope_normalize = False
    ensemble_stack = False
    test_interval_stack_xcorr_(test_dir, cha, inv1, inv2, interval_seconds, window_seconds,
                               clip_to_2std, whitening, one_bit_normalize,
                               envelope_normalize, ensemble_stack)
# end if
