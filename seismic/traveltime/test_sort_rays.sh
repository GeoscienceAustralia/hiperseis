
export PSTHOME=~/Githubz/hiperseis/

export ELLIPCORR=$PSTHOME/ellip-corr/

mkdir $PSTHOME/tmpdir

cd $PSTHOME/tmpdir

python $PSTHOME/seismic/traveltime/sort_rays.py $PSTHOME/tests/testdata/100K_ensemble.p.txt  p_arrivals_sorted1x1.csv P $PSTHOME/seismic/traveltime/param1x1 $PSTHOME/seismic/traveltime/csv_columns.json

python $PSTHOME/seismic/traveltime/sort_rays.py $PSTHOME/tests/testdata/100K_ensemble.s.txt  s_arrivals_sorted1x1.csv S $PSTHOME/seismic/traveltime/param1x1 $PSTHOME/seismic/traveltime/csv_columns.json
