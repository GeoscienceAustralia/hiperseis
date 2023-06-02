"""
Description:
    Miscellaneous functions that don't fit elsewhere

References:

CreationDate:   03/21/23
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     03/21/23   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import subprocess
import os

from mpi4py import MPI

def get_git_revision_hash() -> str:
    """
    Returns the current git hash, if this file is a part of the repository
    """
    prev_path = os.getcwd()
    path = os.path.dirname(os.path.realpath(__file__))
    os.chdir(path)
    result = ''
    try:
        result = subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode('ascii').strip()
    except Exception as e:
        pass
    # end try

    os.chdir(prev_path)
    return result
# end func

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
                self.prev_progress = int(open(self.proc_fn).read())
            # end if
        # end if
    # end func

    def increment(self):
        self.progress += 1
        if(self.restart_mode and (self.prev_progress > 0) and (self.progress < self.prev_progress)):
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

