#!/bin/bash

CURRENT_DIR=`pwd`
WORKING_DIR=~/CWBQuery
TEMPDIR=$WORKING_DIR/tempDir-$(date +"%m-%d-%y-%H-%M-%S")
START_TIME="2005/01/01 00:00:00"
BEGIN_TIME_VAR_FILE=/tmp/beginTime.txt
SEISCOMP3_ARCHIVE=/opt/seiscomp3/var/lib/archive
CWB_QUERY_JAR_FILE=$WORKING_DIR/CWBQuery.jar
CWB_QUERY_SERVER_IP=13.55.154.202
TARBALL_FILE=EdgeCWBRelease.tar.gz

if [ ! -d $WORKING_DIR ]; then
	mkdir $WORKING_DIR
fi

if ! type -p java; then
	echo "This is the first time this script is running ..."
	echo "Installing java ... "
	sudo yum install java-1.7.0-openjdk -y
fi

if [ ! -f $BEGIN_TIME_VAR_FILE ]; then
	echo "Initializing the BEGIN_TIME VARIABLE"
	echo $START_TIME > $BEGIN_TIME_VAR_FILE
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
	mkdir $TEMPDIR-$(date +"%m-%d-%y-%r")
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

TEMPDIR=$WORKING_DIR/tempDir-$(date +"%m-%d-%y-%H-%M-%S")
echo "Creating temporary directory $TEMPDIR for generating the miniseed files..."
mkdir $TEMPDIR
echo "cd-ing into temporary directory $TEMPDIR ..."
cd $TEMPDIR
BEGIN_TIME=`cat $BEGIN_TIME_VAR_FILE`
echo "BEGIN_TIME = $BEGIN_TIME"

JAVA_CMD="java -jar -Xmx1600m $CWB_QUERY_JAR_FILE -h $CWB_QUERY_SERVER_IP -t ms -s '............' -b '$BEGIN_TIME' -d 3600"
echo "Executing the command $JAVA_CMD ... "
eval $JAVA_CMD

for ms_file in `ls *.ms | xargs`; do
	if [ -f ${ms_file} ]; then
		SCART_CMD="scart -I ${ms_file} $SEISCOMP3_ARCHIVE"
		echo "Executing the command $SCART_CMD ..."
		eval $SCART_CMD
		if [ $? -ne 0 ]; then
			echo "The ingesting of the miniseed file ${ms_file} did not succeed! Aborting."
			exit 1
		fi
	fi
done

echo "Finished ingesting the miniseeds for the timewindow $BEGIN_TIME + 1 hour"

echo "Returning to working directory..."
cd $WORKING_DIR
echo "Deleting the temporary directory that was created for generating miniseed file..."
rm -rf $TEMPDIR
NEXT_BEGIN_TIME=$(date '+%Y/%m/%d %H:%M:%S' -d "$BEGIN_TIME 1 hour")
echo "Next begin time = $NEXT_BEGIN_TIME "
echo $NEXT_BEGIN_TIME > $BEGIN_TIME_VAR_FILE

cd $CURRENT_DIR
