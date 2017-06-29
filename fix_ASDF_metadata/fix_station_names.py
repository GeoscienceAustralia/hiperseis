import pyasdf
from os.path import basename, join
import time
import glob
import json
import sys
import numpy as np


code_start_time = time.time()

# =========================== User Input Required =========================== #

# Path to the data
data_path = '/g/data/ha3/Passive/'

# IRIS Virtual Ntework name
virt_net = '_GA_ANUtest'

# FDSN network identifier
FDSNnetwork = 'X5'

# =========================================================================== #


ASDF_path = join(data_path, virt_net, FDSNnetwork, 'ASDF')

#JSON filename for network
JSON_in = join(ASDF_path, FDSNnetwork + '_raw_dataDB.json')
#ASDF filename for network
ASDF_in = join(ASDF_path, FDSNnetwork + '.h5')

#load in asdf file
ds = pyasdf.ASDFDataSet(ASDF_in)

#read the JSON file
f = open(JSON_in, 'r')
_json_dict = json.load(f)

data_set_length = len(_json_dict.keys())

# iterate through all of the waveforms in the ASDF file via the JSON dict
# decode the asscoated logfile and add it into the ASDF file as auxillary data
# get the field notes and add to auxillary data
# then check is serial numbers match and fix station codes
for _i, (key,value) in enumerate(_json_dict.iteritems()):
    print "\r     Parsing Waveform ", _i + 1, ' of ', data_set_length, ' ....',
    sys.stdout.flush()

    # get network of waveform
    net = _json_dict[key]["new_network"]

    # get the associated logfile (full path)
    log_file = _json_dict[key]["log_filename"]

    # get the station of Logfile
    sta = basename(log_file).split('.')[0].replace("LOGFILE", '')

    #check if the logfile is already in the ASDF auxillary data
    data_type = "LogfileData"

    print(data_type)

    # path to logfile data in ASDF auxillary
    path = net + "/" + sta + "/" + "Temperature"

    #add the auxillary data
    ds.add_auxiliary_data(data=np.array([0]),
                          data_type=data_type,
                          path=path,
                          parameters={})


del ds
