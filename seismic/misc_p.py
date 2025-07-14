"""
Description:
    Miscellaneous parallel functions that don't fit elsewhere

References:

CreationDate:   02/06/23
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     02/06/23   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import os
from mpi4py import MPI
from logging import Logger

class ProgressTracker:
    def __init__(self, output_folder, restart_mode=False):
        self.output_folder = output_folder
        self.restart_mode = restart_mode

        self.comm = MPI.COMM_WORLD
        self.nproc = self.comm.Get_size()
        self.rank = self.comm.Get_rank()

        self.prev_progress = 0 # progress from a previous run
        self.progress = 0
        self.proc_fn = os.path.join(output_folder, 'prog.%d.txt' % (self.rank))

        if(self.restart_mode):
            if(os.path.exists(self.proc_fn)):
                try:
                    self.prev_progress = int(open(self.proc_fn).read())
                except:
                    pass
                # end try
            # end if
        # end if
    # end func

    def increment(self):
        self.progress += 1
        if(self.restart_mode and (self.prev_progress > 0) and (self.progress < self.prev_progress)):
            return False
        else:
            try:
                f = open(self.proc_fn, 'w+')
                f.write(str(self.progress))
                f.close()
            except:
                pass
            # end try
            return True
        # end if
    # end func
# end class

def parallel_abort(msg: str, logger:Logger=None):
    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()

    if(logger is not None): logger.error('Aborting job from rank {}: {}'.format(rank, msg))
    comm.Abort()
# end func