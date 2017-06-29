#!/usr/bin/env bash
#module load agdc-py2-prod/1.2.2
#module load pyqt
export PATH=/g/data/ha3/axc547/miniconda2/bin:$PATH
#export PYTHONPATH=/g/data/ha3/axc547/miniconda2/
source activate seismicpy27
python iris2asdf.py
 
