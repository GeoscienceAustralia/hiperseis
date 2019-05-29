from __future__ import absolute_import

import numpy as np
import pytest

pytest.importorskip('mpi4py', 'mpi4py unavailable')
from seismic.traveltime.mpiops import rank, size, comm, run_once

try:
    import h5py
    H5PY = h5py.get_config().mpi
except:
    H5PY = False


def test_helloworld():
    ranks = comm.allgather(rank)
    assert len(ranks) == size


@pytest.mark.skipif(not H5PY, reason='Skipped as parallel h5py is not available')
def test_h5py(random_filename):
    hdf = run_once(random_filename, ext='.hdf5')
    f = h5py.File(hdf, 'w', driver='mpio', comm=comm)
    dset = f.create_dataset('test', (size,), dtype='i')
    dset[rank] = rank
    f.close()

    f = h5py.File(hdf, 'r', libver='latest')
    b = f['test'][:]
    f.close()
    assert len(b) == size
    assert np.allclose(b, np.array(range(size)))
