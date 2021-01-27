# In VDI,  source this file to setup virtual env for hiperseis

git pull

module purge

module load pbs  python3/3.6.7  hdf5/1.8.14p  openmpi/3.1.3  geos/3.5.0  proj/4.9.3  netcdf/4.3.3.1p 


module list
# output will be 
# 1) pbs                    3) szip/2.1               5) netcdf/4.3.3.1p        7) python3/3.6.7          9) geos/3.5.0
# 2) openmpi/3.1.3          4) hdf5/1.8.14p           6) intel-mkl/2018.2.199   8) proj/4.9.3

which python3
python3 -V

source /g/data/ha3/Passive/Software/VENV/vdi_para_h5py/bin/activate


# Run unit-tests may take 10+ minutes to complete
# pytest -v tests/
