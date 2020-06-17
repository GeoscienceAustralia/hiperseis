#!/bin/env python
"""
Description:
    Utility Program for running the demultiplex binary in parallel.
References:

CreationDate:   09/04/19
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     00/04/19   RH
    LastUpdate:     2020-04-10  Fei Zhang  Take over and clean-up the code.

"""

import glob
import os
import random
import subprocess

import click
from mpi4py import MPI


def runprocess(cmd, get_results=False):
    results = []
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in p.stdout:
        if (get_results): results.append(line.strip())
        # else: print line
    # end for

    # for line in p.stderr:
    #    print line
    p.wait()

    return p.returncode, results


# end func

def split_list(lst, npartitions):
    k, m = divmod(len(lst), npartitions)
    return [lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(npartitions)]


# end func

class ProgressTracker:
    def __init__(self, output_folder, restart_mode=False):
        self.output_folder = output_folder
        self.restart_mode = restart_mode

        self.comm = MPI.COMM_WORLD
        self.nproc = self.comm.Get_size()
        self.rank = self.comm.Get_rank()

        self.prev_progress = 0  # progress from a previous run
        self.progress = 0
        self.proc_fn = os.path.join(output_folder, 'prog.%d.txt' % (self.rank))

        if (self.restart_mode):
            if (not os.path.exists(self.proc_fn)):
                raise Exception('Progress file (%s) not found' % (self.proc_fn))
            # end if

            self.prev_progress = int(open(self.proc_fn).read())
        # end if

    # end func

    def increment(self):
        self.progress += 1
        if (self.restart_mode and (self.prev_progress > 0) and (self.progress < self.prev_progress)):
            return False
        else:
            tmpfn = self.proc_fn + '.tmp'
            f = open(tmpfn, 'w+')
            f.write(str(self.progress))
            f.close()
            os.rename(tmpfn, self.proc_fn)

            return True
        # end if
    # end func


# end class

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('input-folder', required=True,
                type=click.Path(exists=True))
@click.argument('output-folder', required=True,
                type=click.Path(exists=True))
@click.option('--extension', default='ms',
              help="File extension for miniseed files; default is 'ms'")
@click.option('--restart', is_flag=True, default=False, help='Restart job')
def process(input_folder, output_folder, extension, restart):
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
    if (rank == 0):
        msfiles = glob.glob(input_folder + '/*.%s' % (extension))

        random.Random(nproc).shuffle(msfiles)  # using nproc as seed so that shuffle produces the same
        # ordering when jobs are restarted.

        msfiles = split_list(msfiles, nproc)
    # end if

    msfiles = comm.scatter(msfiles, root=0)

    path = os.path.dirname(os.path.abspath(__file__))

    # Progress tracker
    progTracker = ProgressTracker(output_folder=output_folder, restart_mode=restart)

    for ifn, fn in enumerate(msfiles):
        if (progTracker.increment()):
            pass
        else:
            print(('Found results for mseed file %s. Moving along..' % (fn)))
            continue
        # end if

        cmd = '%s/demultiplex %s %s %s.%s' % (path, fn, output_folder, rank, ifn)

        ec, results = runprocess(cmd, get_results=True)

        print((ec, results))
    # end for


#######################################################################################################
# Example Cmdline run:  mpirun -np 2 python3 demultiplex.py /Datasets/ /tmp/Output_Miniseeds/
#######################################################################################################
if (__name__ == '__main__'):
    process()
