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
import os, glob, fnmatch, sys
import numpy as np
import logging
logging.basicConfig()

def setup_logger(name, log_file, level=logging.INFO, propagate=False):
    """
    Function to setup a logger; adapted from stackoverflow
    """
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    handler = logging.FileHandler(log_file, mode='w')
    handler.setFormatter(formatter)

    logger = logging.getLogger(name+log_file)
    logger.setLevel(level)
    logger.addHandler(handler)
    logger.propagate = propagate
    return logger
# end func

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

def recursive_glob(treeroot, pattern):
    results = []
    for base, dirs, files in os.walk(treeroot):
        goodfiles = fnmatch.filter(files, pattern)
        results.extend(os.path.join(base, f) for f in goodfiles)
    return results
# end func

def split_list(lst, npartitions):
    k, m = divmod(len(lst), npartitions)
    return [lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(npartitions)]
# end func

def rtp2xyz(r, theta, phi):
    """
    @param r: radius
    @param theta: colat in radians
    @param phi: lon in radians
    @return: x,y,z coordinates on a sphere of radius r
    """
    xout = np.zeros((r.shape[0], 3))
    rst = r * np.sin(theta)
    xout[:, 0] = rst * np.cos(phi)
    xout[:, 1] = rst * np.sin(phi)
    xout[:, 2] = r * np.cos(theta)
    return xout
# end func