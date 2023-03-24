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