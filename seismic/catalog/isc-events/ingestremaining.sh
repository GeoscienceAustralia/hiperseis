#!/bin/bash

for sc3ml in `find 2015 -name *.xml | xargs`; do
	cutfilename=`echo $sc3ml | cut -d'/' -f4`
	if ! grep -q $cutfilename "already-ingested.txt"; then
		scdb -i $sc3ml -d sysop:sysop@localhost/seiscomp3
		#echo $sc3ml
	fi
done
