######  
# login to vdi  and source vdi_env.sh to setup virtual env for hiperseis
#####

module purge

# vdi tested 2021-02-01:
module load pbs  
module load python3/3.7.4 
module load hdf5/1.10.5p 
module load openmpi  
module load geos/3.8.0  
module load proj/6.2.1

module list

source /g/data/ha3/Passive/Software/VENV/para_h5py/bin/activate

which python3
python3 -V


# The output would be:
# Loading python3/3.7.4
#  Loading requirement: intel-mkl/2020.1.217
#Currently Loaded Modulefiles:
# 1) pbs/dummy   2) intel-mkl/2020.1.217   3) python3/3.7.4   4) hdf5/1.10.5p   5) openmpi/4.0.5-debug   6) geos/3.8.0   7) proj/6.2.1  
# /g/data/ha3/Passive/Software/VENV/para_h5py/bin/python3
# Python 3.7.4

# If you use this for the first time, try to run the test-suite below to see if the environment is OK:
# pytest -v tests/


