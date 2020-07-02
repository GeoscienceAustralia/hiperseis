######  login gadi  and source gadi_env.sh to setup virtual env for hiperseis

module purge

module load pbs  python3/3.7.4 hdf5/1.10.5p openmpi/3.1.4  geos/3.8.0  proj/6.2.1

module list

source /g/data/ha3/Passive/Software/VENV/para_h5py/bin/activate

which python3
python3 -V

# pytest -v tests/

#( para_h5py) fxz547@gadi-login-02 ~/github/hiperseis (develop) $ module list
#Currently Loaded Modulefiles:
# 1) pbs   2) intel-mkl/2019.3.199   3) python3/3.7.4   4) hdf5/1.10.5p   5) openmpi/3.1.4   6) geos/3.8.0   7) proj/6.2.1  


