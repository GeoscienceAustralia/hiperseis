import json
from os.path import join
import glob

# =========================== User Input Required =========================== #

# Path to the data
data_path = '/g/data/ha3/Passive/'

# IRIS Virtual Ntework name
virt_net = '_ANU'

# FDSN network identifier
FDSNnetwork = '7F(2013-2014)'

# =========================================================================== #

path_field_notes = join(data_path, virt_net, FDSNnetwork, 'network_metadata/field_notes/')

# get field notes text file
field_notes_txt = glob.glob(join(path_field_notes, "*.txt"))[0]

print(field_notes_txt)

with open(field_notes_txt, 'r') as f:
    notes = f.readlines()


field_notes_dict = {}

# go through lines
for _i, line in enumerate(notes):
    pass