#!/usr/bin/env python
"""
Load inventory file from hdf5 format and export station location plots to graphics files.

Usage:
    python generate_network_plots.py inventory_20190206.h5

where inventory_20190206.h5 here is an example file name.  The output folder for the
graphics files is inferred from the inventory file name, which would be "inventory_20190206"
in this example.
"""

import os
import sys
import pandas as pd

from seismic.inventory.plotting import save_network_local_plots, save_station_local_plots

try:
    import tqdm
    show_progress = True
except:
    show_progress = False
    print("Run 'pip install tqdm' to see progress bar.")

# Set up file names
infile_name = sys.argv[1]
print("Loading file {0}".format(infile_name))
basename, _ = os.path.splitext(infile_name)
_, basename = os.path.split(basename)
folder_name = "plots." + basename
print("Saving plots to {0}".format(folder_name))

# Read inventory from hdf5 file
db = pd.read_hdf(infile_name, mode='r', key='inventory')

if show_progress:
    pbar = tqdm.tqdm(total=len(db), ascii=True)
    progressor = pbar.update
else:
    progressor = None

# Save network plots to files
save_network_local_plots(db, folder_name, progressor=progressor)
if show_progress:
    pbar.close()
