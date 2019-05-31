#!/bin/env python
"""
Description:
    Small utility for running the demultiplex binary in parallel.
References:

CreationDate:   09/04/19
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     00/04/19   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

from mpi4py import MPI
import os, sys

from os.path import join, exists
from collections import defaultdict
import numpy as np
from obspy import Stream, Trace, UTCDateTime
from multiprocessing import Pool, TimeoutError
import pyasdf
from obspy.core.trace import Trace
import click
import glob
import subprocess
import random

def runprocess(cmd, get_results=False):
    results = []
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in p.stdout:
        if (get_results): results.append(line.strip())
        #else: print line
    # end for

    #for line in p.stderr:
    #    print line
    p.wait()

    return p.returncode, results
# end func

def split_list(lst, npartitions):
    random.shuffle(lst)
    k, m = divmod(len(lst), npartitions)
    return [lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(npartitions)]
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('input-folder', required=True,
                type=click.Path(exists=True))
@click.argument('output-folder', required=True,
                type=click.Path(exists=True))
@click.option('--extension', default='ms',
              help="File extension for miniseed files; default is 'ms'")
def process(input_folder, output_folder, extension):
    """
    INPUT_FOLDER: Path to input folder containing mseed files to be demultiplexed \n
    OUTPUT_FOLDER: Output folder \n

    Example usage:
    mpirun -np 2 python demultiplex.py ./ /tmp/output

    Note:
    """

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()

    msfiles = None
    if(rank == 0):
        msfiles = glob.glob(input_folder + '/*.%s'%(extension))
        msfiles = split_list(msfiles, nproc)
    # end if

    msfiles = comm.scatter(msfiles, root=0)

    path = os.path.dirname(os.path.abspath(__file__))

    for ifn, fn in enumerate(msfiles):
        cmd = '%s/demultiplex %s %s %s.%s'%(path, fn, output_folder, rank, ifn)
        
        ec, results = runprocess(cmd, get_results=True)

        print ((ec, results))
    # end for
# end func

if (__name__ == '__main__'):
    process()
# end if
