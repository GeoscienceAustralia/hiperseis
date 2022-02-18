#!/usr/bin/env python
#
# Test that the user's setup on Gadi is capable of importing and using require libraries
# for x-correlation in mpi runtime environment.
#     module purge
#     module load pbs
#     module load python3-as-python
#     module load openmpi/3.1.4
#     module load hdf5/1.10.5p

import sys

def test_setup():
    try:
        print("Testing mpi4py import")
        from mpi4py import MPI
        print("SUCCESS!")
    except:
        print("Failed to import mpi4py")
        sys.exit(1)

    try:
        print("Testing h5py import")
        import h5py
        print("SUCCESS!")
    except:
        sys.exit(1)

    try:
        print("Testing obspy import")
        import obspy
        print("SUCCESS!")
    except:
        sys.exit(1)

    try:
        print("Testing pyasdf import")
        import pyasdf
        print("SUCCESS!")
    except:
        sys.exit(1)

    try:
        print("Testing netCDF4 import")
        from netCDF4 import Dataset
        print("SUCCESS!")
    except:
        sys.exit(1)

    try:
        print("Opening ASDF file for reading")
        afile = pyasdf.ASDFDataSet('../../tests/test_seismic/xcorqc/data/test_data_CMSA.h5', mode='r')
        print("SUCCESS!")
    except:
        print("Failed opening ASDF file")
        sys.exit(1)

    print(afile)

    try:
        print("Reading trace from ASDF file")
        stream, inv = afile.get_data_for_tag('AU.CMSA', 'raw_recording')
        print(inv)
        print(stream[0])
        print("SUCCESS!")
    except:
        print("Failed to read data from ASDF file")
        sys.exit(1)


if __name__ == '__main__':
    test_setup()

