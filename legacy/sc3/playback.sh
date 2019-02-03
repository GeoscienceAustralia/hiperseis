#!/usr/bin/env bash

if [ "$#"  -lt 1 ]; then
echo "Usage: $0 [mseed-volume]";
exit 0;
fi
set -x; seiscomp stop; seiscomp start spread scmaster
DBFLAG="mysql://sysop:sysop@localhost/seiscomp3"
VERBOSITY="-v"

evid=$1
if [[ $1 == *"/"* ]]; then
	evid=`echo $1 | cut -d"/" -f3`  
fi

scevtstreams -E $1 -d mysql://sysop:sysop@localhost/seiscomp3 -L 0 -m 3600 > $evid-stream.txt
net_sta_cha_exprs=`cat $evid-stream.txt | cut -d";" -f3 | xargs`
for expression in $net_sta_cha_exprs; do
	if [[ $expression == \.* ]]; then
		sta=`echo $expression | cut -d"." -f2`
		potLoc=`find /opt/seiscomp3/var/lib/archive/2015/ -name $sta -type d | xargs`
		if [[ ! -z "${potLoc// }" && ! $potLocs == *" "* ]]; then
			resSta=`echo $potLoc | rev | cut -d"/" -f1 | rev`
			resNet=`echo $potLoc | rev | cut -d"/" -f2 | rev`
			if [ $resSta == $sta ]; then
				echo "replacing $expression with $resNet$expression ... in the file $evid-stream.txt"
				sed -i -e "s/;$expression/;$resNet$expression/g" $evid-stream.txt
			fi
		fi
	fi
done

scart -dsvE --list $evid-stream.txt /opt/seiscomp3/var/lib/archive > $evid-perm.ms
beginTime=`head -1 $evid-stream.txt | cut -d";" -f1`
endTime=`head -1 $evid-stream.txt | cut -d";" -f2`
scart -dsvE -t "$beginTime~$endTime" -n 7G /opt/seiscomp3/var/lib/archive > $evid-temp.ms
python mergePermTemp.py $evid-perm.ms $evid-temp.ms $evid-merged.ms

scautopick --ep --playback -I file://$evid-merged.ms -d $DBFLAG > $evid-picks.xml
scautoloc --ep $evid-picks.xml -d $DBFLAG $VERBOSITY > $evid-origins.xml
grep origin $evid-origins.xml
if [ $? -eq 0 ]; then
	echo "Origins for event $1 generated successfully ..."
	echo $1 >> origin-success.txt
fi
scamp --ep $evid-origins.xml -I file://$evid-merged.ms -d $DBFLAG $VERBOSITY > $evid-amps.xml
scmag --ep $evid-amps.xml -d $DBFLAG $VERBOSITY > $evid-mags.xml
scevent --ep $evid-mags.xml -d $DBFLAG $VERBOSITY > $evid-events.xml

# scdb -i events.xml -d $DBFLAG $VERBOSITY
set -x; seiscomp stop; seiscomp start
