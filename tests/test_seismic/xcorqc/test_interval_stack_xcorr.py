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
location_code = ''

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
output_folder = str(tempfile.mkdtemp())

@pytest.fixture(params=['BHZ', '00T'])
def cha(request):
    return request.param

@pytest.fixture(params=[86400, 3600])
def interval_seconds(request):
    return request.param

@pytest.fixture(params=[0, 0.1])
def window_overlap(request):
    return request.param

@pytest.fixture(params=[3600, 1800])
def window_seconds(request):
    return request.param

@pytest.fixture(params=[False, True])
def ensemble_stack(request):
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

def test_interval_stack_xcorr(cha, inv1, inv2, interval_seconds, window_seconds,
                               window_overlap, whitening, ensemble_stack):
    """
    Note that these tests exclude the application of 'stacking', 'window_buffer_length' --
    they are tested in the test_xcorr.py
    """
    start_time = '2011-03-11T00:00:00Z'
    end_time   = '2011-03-14T00:00:00Z'

    tag = '%s_is%d_ws%d_wo%.2f_w%d_es%d_resp%d_resp%d'%(
            cha, interval_seconds, window_seconds,
            window_overlap, whitening, ensemble_stack,
            1 if inv1 else 0,
            1 if inv2 else 0)

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
                       read_ahead_window_seconds=interval_seconds,
                       interval_seconds=interval_seconds,
                       window_seconds=window_seconds,
                       window_overlap=window_overlap,
                       flo=0.01, fhi=0.125,
                       whitening=whitening,
                       ensemble_stack=ensemble_stack,
                       outputPath=output_folder, verbose=2, tracking_tag=tag)

    # Read result
    fn = os.path.join(output_folder, '%s.%s.%s.%s.%s.%s.%s.nc'%(netsta1, location_code, cha,
                                                                netsta2, location_code, cha, tag))
    dc = Dataset(fn)
    xcorr_c = dc.variables['xcorr'][:]

    # Read expected
    fn = '%s/%s.%s.%s.%s.%s.%s.%s.nc'%(expected_folder, netsta1, location_code, cha,
                                       netsta2, location_code, cha, tag)
    de = Dataset(fn)
    xcorr_e = de.variables['xcorr'][:]

    rtol = 1e-3
    atol = 1e-3

    assert np.allclose(xcorr_c, xcorr_e, rtol=rtol, atol=atol)
# end func

if __name__ == '__main__':
    cha = 'BHZ'
    inv1 = inv.select(network='AU', station='ARMA')
    inv2 = inv.select(network='AU', station='QLP')
    interval_seconds = 86400
    window_seconds = 3600
    window_overlap = 0.1
    whitening = False
    ensemble_stack = True

    test_interval_stack_xcorr(cha, inv1, inv2, interval_seconds, window_seconds,
                              window_overlap, whitening, ensemble_stack)
# end if
