#!/bin/bash

marchEvents=`scevtls -d mysql://sysop:sysop@localhost/seiscomp3 --begin "2015-03-01 00:00:00" --end "2015-03-30 23:59:59" | xargs`

for event in $marchEvents; do
	./playback.sh $event
done
