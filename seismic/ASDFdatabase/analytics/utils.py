"""
Description:
    Implements miscellaneous functionalities for generating waveform-analytics

References:

CreationDate:   18/01/2024
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     07/02/23   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

from obspy import read_inventory
from seismic.inventory.response import ResponseFactory
from multiprocess import Manager, freeze_support
import os

class ProgressTracker(object):
    def __init__(self, manager: Manager):
        self.max = manager.Value('max', 0)
        self.i = manager.Value('i', 0)
        self.lock = manager.Lock()
    # end func

    def initialize(self, max_value, initial_value=0):
        with self.lock:
            self.max.value = max_value
            self.i.value = initial_value
            # end with
    # end func

    def increment(self):
        with self.lock:
            if (self.i.value < self.max.value):
                self.i.value += 1
            # end if
        # end with
    # end func

    def now(self):
        with self.lock:
            return self.i.value, self.max.value
        # end with
    # end func
# end class

def get_response(input_file, network=None, station=None, location=None, channel=None):
    def is_sqlite_db(db_fn):
        """Check if a file is a valid SQLite db"""
        if not os.path.isfile(db_fn):
            return False
        try:
            with open(db_fn, "rb") as f:
                header = f.read(16)
            return header == b"SQLite format 3\000"
        except Exception:
            return False
        # end try
    # end func

    result = None
    if('xml' in input_file.lower()):
        inv = None
        try:
            inv = read_inventory(input_file)
        except Exception as e:
            print('Failed to read inventory file {} with error: {}'.format(input_file, e))
        # end try

        if(inv is not None):
            inv = inv.select(network=network, station=station, location=location,
                             channel=channel)
            if(inv is not None):
                seedid = inv.get_contents()['channels'][0]
                resp_obj = inv.get_response(seedid,
                                            inv.networks[0].stations[0].channels[0].start_date)
                result = resp_obj
        # end if
    elif(is_sqlite_db(input_file)):
        nslc = None
        if( (network is not None and len(network) > 0) and
            (station is not None and len(station) > 0) and
            (location is not None) and
            (channel is not None and len(channel) > 0) ):

            nslc = '.'.join((network, station, location, channel))

            rf = ResponseFactory()
            rf.createFromDB(input_file)
            result = rf.getResponse(nslc)
        # end if
    else:
        resp_name = 'resp'
        # read resp file
        resp_inv = None
        try:
            resp_inv = read_inventory(input_file, format='RESP')
        except Exception as e:
            print('Failed to read RESP file {} with error: {}'.format(input_file, e))
        # end try

        if(resp_inv is not None):
            rf = ResponseFactory()
            rf.createFromInventory(resp_name, resp_inv)

            result = rf.getResponse(resp_name)
        # end if
    # end if

    return result
# end func
