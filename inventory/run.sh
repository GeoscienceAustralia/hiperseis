#!/bin/bash
#PBS -l select=1:ncpus=1:mem=8gb:cput=7200s
#PBS -N StationMetadataCleanup
time python ./engd2stxml.py
