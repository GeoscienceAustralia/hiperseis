#!/bin/env python
"""
Description:
    Tests various aspects of correlator.py

References:

CreationDate:   11/12/2021
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     11/12/2021   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
from seismic.xcorqc.correlator import process
import os
from netCDF4 import Dataset
import numpy as np
import tempfile

# Prepare input
netsta1 = 'II.WRAB'
netsta2 = 'AU.MAW'
location_code1 = '00'
location_code2 = '10'

path = os.path.dirname(os.path.abspath(__file__))

# Initialize input data
files_dir = tempfile.mkdtemp(suffix='_test_correlator')
#files_dir = '/tmp'
asdf_file_list1 = os.path.join(files_dir, 'asdf_file_list1.txt')
asdf_file_list2 = os.path.join(files_dir, 'asdf_file_list2.txt')
pref_file1 = os.path.join(files_dir, 'pref_file1.txt')
pref_file2 = os.path.join(files_dir, 'pref_file2.txt')

f1 = open(asdf_file_list1, 'w+')
f2 = open(asdf_file_list2, 'w+')
f1.write('%s/data/test_data_WRAB.h5\n'%(path))
f2.write('%s/data/test_data_MAW.h5\n'%(path))
f1.close()
f2.close()

f3 = open(pref_file1, 'w+')
f4 = open(pref_file2, 'w+')
f3.write('%s, %s\n'%(netsta1, location_code1))
f4.write('%s, %s\n'%(netsta1, location_code2))
f3.close()
f4.close()

fds1 = FederatedASDFDataSet(asdf_file_list1)
fds2 = FederatedASDFDataSet(asdf_file_list2)

expected_folder = '%s/data/expected/'%(path)
output_folder = str(tempfile.mkdtemp())
os.mkdir(os.path.join(output_folder, 'stacked'))
os.mkdir(os.path.join(output_folder, 'unstacked'))
#output_folder = '/tmp'

def test_correlator():
    start_time = '2006-11-03T00:00:00'
    end_time   = '2006-11-04T00:00:00'
    cha = 'BHZ'

    # Stacked
    for loc_pref, loc_code in zip([pref_file1, pref_file2], [location_code1, location_code2]):
        curr_output_folder = os.path.join(output_folder, 'stacked')
        curr_expected_folder = os.path.join(expected_folder, 'stacked')

        process(asdf_file_list1, asdf_file_list2,
                curr_output_folder,
                86400, 3600, 0.1,
                0.05, 86400, 4, 0.05, -1, None, None, 0.002, 2, netsta1,
                netsta2, None, start_time, end_time, None, 'vel',
                50, False, True, 0.02, True, loc_pref,
                '*Z', '*N', '*E', '*Z', '*N', '*E', 'z', False, False,
                True, False, False, True, None)


        # Read result
        fn = os.path.join(curr_output_folder, '%s.%s.%s.%s.%s.%s.nc'%(netsta1, loc_code, cha,
                                                                      netsta2, '', cha))
        dc = Dataset(fn)
        xcorr_c = dc.variables['xcorr'][:]

        # Read expected
        fn = '%s/%s.%s.%s.%s.%s.%s.nc'%(curr_expected_folder, netsta1, loc_code, cha,
                                        netsta2, '', cha)
        de = Dataset(fn)
        xcorr_e = de.variables['xcorr'][:]

        rtol = 1e-5
        atol = 1e-5

        assert np.allclose(xcorr_c, xcorr_e, rtol=rtol, atol=atol)
    # end for

    # Unstacked
    for loc_pref, loc_code in zip([pref_file1, pref_file2], [location_code1, location_code2]):
        curr_output_folder = os.path.join(output_folder, 'unstacked')
        curr_expected_folder = os.path.join(expected_folder, 'unstacked')

        process(asdf_file_list1, asdf_file_list2,
                curr_output_folder,
                86400, 3600, 0.1,
                0.05, 86400, 4, 0.05, -1, None, None, 0.002, 2, netsta1,
                netsta2, None, start_time, end_time, None, 'vel',
                50, False, True, 0.02, True, loc_pref,
                '*Z', '*N', '*E', '*Z', '*N', '*E', 'z', False, False,
                False, False, False, True, None)


        # Read result
        fn = os.path.join(curr_output_folder, '%s.%s.%s.%s.%s.%s.nc'%(netsta1, loc_code, cha,
                                                                      netsta2, '', cha))
        dc = Dataset(fn)
        xcorr_c = dc.variables['xcorr'][:]

        # Read expected
        fn = '%s/%s.%s.%s.%s.%s.%s.nc'%(curr_expected_folder, netsta1, loc_code, cha,
                                        netsta2, '', cha)
        de = Dataset(fn)
        xcorr_e = de.variables['xcorr'][:]

        rtol = 1e-5
        atol = 1e-5

        assert np.allclose(xcorr_c, xcorr_e, rtol=rtol, atol=atol)
    # end for
# end func

if __name__ == '__main__':
    test_correlator()
# end if
