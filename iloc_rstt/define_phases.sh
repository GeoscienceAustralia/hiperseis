#!/bin/bash
INPUT_DIR=/home/centos/mtIsa_events/2010
OUT_DIR=/home/centos/phase-defined-events

for ev in `find $INPUT_DIR -type f -name "*.xml" | xargs`
do
        scdb -i $ev -d mysql://sysop:sysop@localhost/seiscomp3
done

mkdir -p $OUT_DIR
grep "<event " $INPUT_DIR/*.xml | cut -d'"' -f2 | sed -e 's/\(.*\)/\1 update_db=1 do_gridsearch=0/' > instr_file_2010.txt
iloc sc3db < instr_file_2010.txt > log_2010.out
for id in `grep "<event " $INPUT_DIR/*.xml | cut -d'"' -f2 | xargs`; do scxmldump -fPAMF -E "$id" -o $OUT_DIR/${id//[:\/]/\.}.xml -d mysql://sysop:sysop@localhost/seiscomp3; done
