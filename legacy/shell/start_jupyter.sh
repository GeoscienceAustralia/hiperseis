#! /bin/env bash

LOGFILE=tempworks/_jupyter.log
jupyter notebook &> $LOGFILE &

