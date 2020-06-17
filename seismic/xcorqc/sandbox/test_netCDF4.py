# This works for h5py==2.8.0 + netCDF4==1.5.3, OR with h5py==2.9.0 (or later) and netCDF4==1.5.2,
# but not with h5py==2.9.0 (or later) and netCDF4==1.5.3

import h5py  # For h5py 2.9.0 or 2.10.0, comment this line to fix, or move it after netCDF4 import
from netCDF4 import Dataset


def test_main():
    root_grp = Dataset('result.nc', 'w', format='NETCDF4')
    root_grp.createDimension('dummy', 1)  # comment this line to fix
    x = root_grp.createVariable('x', 'f4')
    val = 1.0
    x[:] = val
    root_grp.close()


if __name__ == '__main__':
    print(h5py.version.info)
    test_main()
