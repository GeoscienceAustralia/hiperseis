import logging
from mpi4py import MPI

log = logging.getLogger(__name__)

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


def run_once(f, *args, **kwargs):
    """
    Run a function on one node and then broadcast result to all.

    :param str f: The function to be evaluated. Can take arbitrary arguments
                and return anything or nothing
    :param str args: Other positional arguments to pass on to f (optional)
    :param str kwargs: Other named arguments to pass on to f (optional)

    :return: The value returned by f.
    :rtype: unknown
    """
    if rank == 0:
        f_result = f(*args, **kwargs)
    else:
        f_result = None
    result = comm.bcast(f_result, root=0)
    return result
