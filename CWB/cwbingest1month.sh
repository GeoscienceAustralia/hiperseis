#!/bin/bash

#java -jar -Xmx1600m ~/CWBQuery/CWBQuery.jar -h 54.153.144.205 -t dcc -s 'AUAS31.BHN' -b '2015/03/10 00:00:00' -d 31d
CURRENT_DIR=`pwd`
WORKING_DIR=~/CWBQuery-$(date +"%m-%d-%y-%H-%M-%S-%N")
TEMPDIR=$WORKING_DIR/tempDir-$(date +"%m-%d-%y-%H-%M-%S-%N")
START_TIME="2015/03/01 00:00:00"
SEISCOMP3_ARCHIVE=/opt/seiscomp3/var/lib/archive
CWB_QUERY_JAR_FILE=$WORKING_DIR/CWBQuery.jar
CWB_QUERY_SERVER_IP=54.153.144.205
TARBALL_FILE=EdgeCWBRelease.tar.gz

mkdir -p $WORKING_DIR

if ! type -p java; then
        echo "This is the first time this script is running ..."
        echo "Installing java ... "
        sudo yum install java-1.7.0-openjdk -y
fi

if ! type -p scart ; then
        echo "Seiscomp3 is not installed or there is some problem with the scart tool! Aborting."
        exit 1
fi

cd $WORKING_DIR

JARS=`ls *.jar | xargs`
if [ -z "${JARS// }" ]; then
        echo "Jars have not been downloaded yet. Downloading ..."
        echo "Creating temporary directory $TEMPDIR for downloading and uncompressing jar files..."
        mkdir $TEMPDIR
        echo "cd-ing into the temporary directory $TEMPDIR ..."
        cd $TEMPDIR

        if [ ! -f  $TARBALL_FILE ]; then
                wget ftp://hazards.cr.usgs.gov/CWBQuery/$TARBALL_FILE
                tar -xzvf $TARBALL_FILE
                echo "Copying jar files to working directory"
                cp `find . -name *.jar | xargs` $WORKING_DIR
        fi

        echo "Returning to working directory..."
        cd $WORKING_DIR
        echo "Deleting the temporary directory that was created for uncompressing jar files..."
        rm -rf $TEMPDIR
fi

while IFS= read -r net_sta_cha_expr
do
        miniseed_generation_success=0
        TEMPDIR=$WORKING_DIR/$net_sta_cha_expr-tempDir-$(date +"%m-%d-%y-%H-%M-%S-%N")
        echo "Creating temporary directory $TEMPDIR for generating the miniseed files..."
        mkdir $TEMPDIR
        echo "cd-ing into temporary directory $TEMPDIR ..."
        cd $TEMPDIR

        JAVA_CMD="java -jar -Xmx1600m $CWB_QUERY_JAR_FILE -h $CWB_QUERY_SERVER_IP -t dcc -s \"$net_sta_cha_expr\" -b \"2015/03/01 00:00:00\" -d 31d > /tmp/$net_sta_cha_expr-cwbquery-log.out 2>/tmp/$net_sta_cha_expr-cwbquery-log.err"
        echo "Executing the command $JAVA_CMD ... "

        eval $JAVA_CMD &
        child_pid=$!
        echo "The pid of the java command process just spawned is $child_pid ... "
	count=0
        while : ; do
		echo "Checking with kill -0 ..."
                kill -0 $child_pid
		retCode=$?
		echo "return code of kill -0 = $retCode ..."
                if [ $retCode -ne 0 ]; then
                        echo "Miniseed generation succeeded"
                        miniseed_generation_success=1
                        break
                fi
                sleep 30
                size=`du -h $TEMPDIR | xargs | cut -d" " -f1 | tr -dc '0-9.' | cut -f1 -d"."`
                unit=`du -h $TEMPDIR | xargs | cut -d" " -f1 | tr -d '0-9.'`
                if [[ ($unit == *"G"* && $size -gt 5) || $count -gt 30 ]]; then
	                echo "Killing with kill -9 ..."
                        kill -9 $child_pid
                        break
                fi
		count=`expr $count + 1`
        done

        if [ $miniseed_generation_success -eq 1 ]; then
                for ms_file in `ls *.msd | xargs`; do
                        if [ -f ${ms_file} ]; then
                                SCART_CMD="scart -I ${ms_file} $SEISCOMP3_ARCHIVE"
                                echo "Executing the command $SCART_CMD ..."
                                eval $SCART_CMD
                                if [ $? -ne 0 ]; then
                                        echo "The ingesting of the miniseed file ${ms_file} did not succeed! Aborting."
                                        exit 1
                                fi
                        fi
                echo "Finished ingesting the miniseeds for $net_sta_cha_expr"
                done
        fi

        echo "Returning to working directory..."
        cd $WORKING_DIR
        echo "Deleting the temporary directory that was created for generating miniseed file..."
        rm -rf $TEMPDIR
done < $1

cd $CURRENT_DIR

