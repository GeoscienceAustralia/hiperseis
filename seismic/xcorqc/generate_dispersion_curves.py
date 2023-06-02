#!/bin/env python
"""
Description:
    Runs Rhys Hawkins' code in parallel to generate dispersion curves
    based on cross-correlations of station-pairs. Note that this script
    call shell scripts that are expected to be in the current working
    directory.

    todo: remove dependence on shell scripts.

References:

CreationDate:   10/01/20

Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     10/01/20   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

from mpi4py import MPI
import subprocess

import click
import psutil
from seismic.misc import split_list

def kill(proc_pid):
    process = psutil.Process(proc_pid)
    for proc in process.children(recursive=True):
        proc.kill()
    # end for
    process.kill()
# end func

def runprocess(cmd, get_results=False):
    results = []
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    try:
        pstdout, pstderr = p.communicate(timeout=6000)
        for line in pstdout.splitlines():
            if (get_results): results.append(line.strip())
            #else: print line
        # end for
    except subprocess.TimeoutExpired:
        print ('Command %s FAILED. KILLING Proc: %d'%(cmd, p.pid))
        kill(p.pid)
        p.returncode = 1
    # end try

    return p.returncode, results
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('station-pair',
                type=click.Path('r'))
@click.argument('output-dir',
                type=click.Path('r'))
def main(station_pair, output_dir):
    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()

    workload = []
    if(rank==0):
        f = open(station_pair, 'r')
        for line in f:
            sp = line.split(' ')[1]
            workload.append([sp, output_dir])
        # end for
        workload = split_list(workload, nproc)
    # end if
    workload = comm.bcast(workload, root=0)

    for item in workload[rank]:
        cmds = ['./01_create_initial_target_phase_rayleigh.sh %s %s'%(item[0], item[1]),
                './02_fit_initial_target_phase_rayleigh.sh %s %s'%(item[0], item[1]),
                './03_fit_bessel_rayleigh.sh %s %s'%(item[0], item[1])]

        for cmd in cmds:
            rc, r = runprocess(cmd, get_results=True)
            print ('command: %s, return code %d: '%(cmd,rc))
            #for line in r: print (line.strip())
        # end for
    # end for
# end func

if __name__=='__main__':
    main()
# end if
