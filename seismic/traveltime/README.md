## How to generate input data for 3D Travel Time Tomography inversion

### Install `HiPerSeis` software.



-  ``git clone https://github.com/GeoscienceAustralia/hiperseis.git``
- ``cd hiperseis``
- ``git submodule init``
- ``git submodule update``

- ``pip install -v --user -e .`` (into user's own home ~/.local/lib/python2.7/site-packages/)


If you are using the library without ``pip install``-ing it, then make sure the ``hiperseis`` folder
is in your ``PYTHONPATH`` environment variable.  E.g. in bash shell, if ``~/dev/hiperseis`` is where
``hiperseis`` was checked out:

- ``export PYTHONPATH=~/dev/hiperseis:$PYTHONPATH``


### The Test Events Arrivals CSV file

- tests/testdata    


### Run Sort Rays Program
    
    $ export PSTHOME=/g/data/ha3/fxz547/Githubz/hiperseis/
    
    $ export ELLIPCORR=$PSTHOME/ellip-corr/
    
    $ python $PSTHOME/seismic/traveltime/sort_rays.py /path2/p_arrivals.txt  p_arrivals_sorted1x1.csv P $PSTHOME/seismic/traveltime/param1x1 $PSTHOME/seismic/traveltime/csv_columns.json
    
    This will generate a file  `p_arrivals_sorted1x1.csv_inv.txt` as input for inversion program.
    
    $ python $PSTHOME/seismic/traveltime/sort_rays.py /path2/s_arrivals.txt  s_arrivals_sorted1x1.csv S $PSTHOME/seismic/traveltime/param1x1 $PSTHOME/seismic/traveltime/csv_columns.json
    
    This will generate a file `s_arrivals_sorted1x1.csv_inv.txt` as  input for inversion program.
    
### Examples:

     python $PSTHOME/seismic/traveltime/sort_rays.py $PSTHOME/tests/testdata/100K_ensemble.p.txt  p_arrivals_sorted1x1.csv P $PSTHOME/seismic/traveltime/param1x1 $PSTHOME/seismic/traveltime/csv_columns.json
     
     python $PSTHOME/seismic/traveltime/sort_rays.py $PSTHOME/tests/testdata/100K_ensemble.s.txt  s_arrivals_sorted1x1.csv S $PSTHOME/seismic/traveltime/param1x1 $PSTHOME/seismic/traveltime/csv_columns.jso
    
