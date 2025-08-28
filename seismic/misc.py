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
import traceback
logging.basicConfig()

def setup_logger(name, log_file=None, level=logging.INFO, propagate=False):
    """
    Function to setup a logger; adapted from stackoverflow
    """
    class ConditionalFormatter(logging.Formatter):
        def format(self, record):
            if hasattr(record, 'simple') and record.simple:
                return record.getMessage()
            else:
                return logging.Formatter.format(self, record)
            # end if
        # end func
    # end class

    handler = None
    if(log_file):
        handler = logging.FileHandler(log_file, mode='w')
    # end if

    formatter = ConditionalFormatter('%(asctime)s %(levelname)s %(message)s')

    logger = logging.getLogger(name+log_file if log_file else '')
    logger.setLevel(level)
    if (handler is not None):
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    # end if
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

def read_key_value_pairs(file_path:str, keys:list, strict=False)->dict:
    """
    Reads a text file containing colon-separated key-value pairs and returns a dictionary with values for specified keys.
    Raises a ValueError if any of the specified keys are not found.

    :param file_path: Path to the text file.
    :param keys: A list of keys to search for in the file.
    :param strict: Ensures all keys are found in the file
    :return: A dictionary with the specified keys and their corresponding values.
    :raises ValueError: If any key is not found in the file, if strict is set to True.
    """
    result = {}
    keys_found = set()

    f = None
    try:
        f = open(file_path, 'r')
    except Exception as e:
        raise e
    else:
        with open(file_path, 'r') as f:
            for line in f:
                # Split the line by colon to get the key-value pair
                if ':' in line:
                    key, value = line.strip().split(':', 1)
                    key, value = key.strip(), value.strip()  # Clean up extra spaces
                    # If the key is in the provided list, add it to the result
                    if key in keys:
                        result[key] = value
                        keys_found.add(key)
                    # end if
            # end for
        # end with
    # end try

    # Check if any keys were not found
    missing_keys = set(keys) - keys_found
    if missing_keys:
        raise ValueError(f"The following keys were not found in the file: {', '.join(missing_keys)}")
    # end if

    return result
# end func

def print_exception(e: Exception):
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, e, fname, exc_tb.tb_lineno)
    print(traceback.format_exc())
# end func