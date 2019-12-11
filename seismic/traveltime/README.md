## How to generate input data for 3D Travel Time Tomography inversion

### Install `HiPerSeis` software.


`git clone https://github.com/GeoscienceAustralia/hiperseis`


### Get Events Arrivals CSV file
    
Assume this data is obtained from upstream pipeline modules: phase-picking programms, SeisComp3 iLoc and events dump, etc. 


### Run Sort Rays Program
    
    $ export PSTHOME=/g/data/ha3/fxz547/Githubz/hiperseis/
    
    $ export ELLIPCORR=$PSTHOME/ellip-corr/
    
    $ python $PSTHOME/seismic/traveltime/sort_rays.py /path2/p_arrivals.txt  p_arrivals_sorted1x1.csv P $PSTHOME/seismic/traveltime/param1x1 $PSTHOME/seismic/traveltime/csv_columns.json
    
    $ python $PSTHOME/seismic/traveltime/sort_rays.py /path2/s_arrivals.txt  s_arrivals_sorted1x1.csv S $PSTHOME/seismic/traveltime/param1x1 $PSTHOME/seismic/traveltime/csv_columns.json

The results will be in 

    p_arrivals_sorted1x1.csv_inv.txt
    
    AND
    
    s_arrivals_sorted1x1.csv_inv.txt
    
