#! /bin/env bash

# Purpose:  jupyter-notebook control script
#
# Usage:    ./jupyter-service-cntl.sh {start|stop|status|restart}
#
# Author:   Fei Zhang
# Date:     2018-09-25
#
# set environment vars for jupyter to use.
# u46 anaconda and gdal env
#export PATH=/g/data/u46/fxz547/anaconda2/bin:$PATH
#export GDAL_DATA=/g/data/u46/fxz547/anaconda2/share/gdal/
#export PYTHONPATH='/home/547/fxz547/Githubz/agdc-v2'

# Or module load
#module use /g/data/v10/public/modules/modulefiles --append
#module load agdc-py2-prod
#module load agdc-py2-dev

# Alternateivly, 
# dynamically modify your pythonpath in each notebook
# import sys; sys.path.insert("/pythonmodule/path/")


# jupyter-notebook --port=5999 &> jupyter_notebook.log  &
# jupyter-notebook --no-browser   --port=4999 &> jupyter_notebook.log  &

# log file
LOG_FILE=$HOME/tmp/jupyter_notebook.log

PROC_NAME='jupyter-notebook'

status(){ 
    echo Checking if jupyter-notebook is already running? 
    jpid=`pgrep -u $USER -f $PROC_NAME`
    echo The $PROC_NAME pid is: $jpid

    if [ ${jpid}1 -gt 1 ]; then
        echo "It ($jpid) is an integer, will return 1"
        return 1
        #return $jpid  # this number will be cast into [0 256) eg: 8625 -> 177
    else
        echo  It $jpid is empty or not defined, will return 0
        return 0
    fi
}


start() {
    
    echo first check if the service is started already
    echo "If not, do start"
    # code to start app comes here 
    # example: daemon program_name &

    echo `date` > $LOG_FILE

    echo Starting jupyter on $HOSTNAME >> $LOG_FILE 

    echo " --------------------------------------------------" >> $LOG_FILE

    jupyter-notebook &>> $LOG_FILE &

    echo jupyer has been started: `pgrep jupyter`

}

stop() {
    # code to stop app comes here 
    # example: killproc program_name
    echo 'Stopping jupyter now..'
    echo 'Running  kill -9 `pgrep jupyter` '
    kill -9 `pgrep -u $USER jupyter` 

    status
}

case "$1" in 
    start)
       status
       retval=$?

       if [ $retval -eq 0 ]; then
           start
       else
           echo  $PROC_NAME is already running
           echo $retval
       fi
       ;;
    stop)
       stop
       ;;
    restart)
       stop
       start
       ;;
    status)
       # code to check status of app comes here 
       # example: status program_name
       status
       retval=$?

       if [ $retval -eq 0 ]; then
          echo "Not running " 
       else
           echo  $PROC_NAME is already running
           echo $retval
       fi
       ;;
    *)
       echo "Usage: $0 {start|stop|status|restart}"
esac

exit 0

