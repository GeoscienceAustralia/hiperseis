#!/bin/bash
#PBS -P vy72
#PBS -q normal
#PBS -l ncpus=9648
#PBS -l walltime=05:00:00
#PBS -l mem=256GB
#PBS -l wd
#PBS -WMail_Users=alexei.gorbatov@ga.gov.au
#PBS -m ae

mpirun ./run
