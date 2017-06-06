import pytest
import numpy as np
from seismic.mpiops import rank, size, comm, run_once

try:
    import h5py
    H5PY = True
except:
    H5PY = False


def test_helloworld():
    ranks = comm.allgather(rank)
    assert len(ranks) == size


@pytest.mark.skipif(not H5PY, reason='Skipped as h5py is not installed')
def test_h5py(random_filename):
    hdf = run_once(random_filename, ext='.hdf5')
    f = h5py.File(hdf, 'w', driver='mpio', comm=comm)
    dset = f.create_dataset('test', (4,), dtype='i')
    # print('here3', rank)
    dset[rank] = rank
    f.close()

    f = h5py.File(hdf, 'r', libver='latest')
    b = f['test'][:]
    assert len(b) == size
    assert np.allclose(b, np.array(range(size)))
