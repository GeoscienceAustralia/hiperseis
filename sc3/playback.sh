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
	if [[ $expression = \.* ]]; then
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

scautopick --ep --playback -I file://$evid-perm.ms -d $DBFLAG > picks.xml
scautoloc --ep picks.xml -d $DBFLAG $VERBOSITY > origins.xml
scamp --ep origins.xml -I file://$evid-perm.ms -d $DBFLAG $VERBOSITY > amps.xml
scmag --ep amps.xml -d $DBFLAG $VERBOSITY > mags.xml
scevent --ep mags.xml -d $DBFLAG $VERBOSITY > events.xml

# scdb -i events.xml -d $DBFLAG $VERBOSITY
COVERALLS_REPO_TOKEN=fzU0I2CIwPKYFCrkn0ZtsVzSTZUWI7KhJ2
set -x; seiscomp stop; seiscomp start
