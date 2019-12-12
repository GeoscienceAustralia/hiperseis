

# Before running this script, make sure you are in a python enviroment where ellipcorr,  seismic and dependency packages have been installed
# typically a virtualenv through source ...venvname/bin/activate; or conda activate venv_name

export PSTHOME=~/Githubz/hiperseis/
export PYTHONPATH=$PSTHOME:$PYTHONPATH
export ELLIPCORR=$PSTHOME/ellip-corr/

wkdir="$PSTHOME/tmpdir"


if [[ ! -e $wkdir ]]; then
  mkdir $wkdir
fi

cd $wkdir


python $PSTHOME/seismic/traveltime/sort_rays.py $PSTHOME/tests/testdata/100K_ensemble.p.txt  p_arrivals_sorted1x1.csv P $PSTHOME/seismic/traveltime/param1x1 $PSTHOME/seismic/traveltime/csv_columns.json

python $PSTHOME/seismic/traveltime/sort_rays.py $PSTHOME/tests/testdata/100K_ensemble.s.txt  s_arrivals_sorted1x1.csv S $PSTHOME/seismic/traveltime/param1x1 $PSTHOME/seismic/traveltime/csv_columns.json


ls -ltr $wkdir
