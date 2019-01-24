#!/usr/bin/env python

import os
import sys
import pandas as pd

from plotting import saveNetworkLocalPlots, saveStationLocalPlots

try:
    import tqdm
    show_progress = True
except:
    show_progress = False
    print("Run 'pip install tqdm' to see progress bar.")


infile_name = sys.argv[1]
print("Loading file {0}".format(infile_name))
basename, _ = os.path.splitext(infile_name)
_, basename = os.path.split(basename)
folder_name = "plots." + basename
print("Saving plots to {0}".format(folder_name))

db = pd.read_hdf(infile_name, mode='r', key='inventory')

if show_progress:
    pbar = tqdm.tqdm(total=len(db), ascii=True)
    progressor = pbar.update
else:
    progressor = None
saveNetworkLocalPlots(db, folder_name, progressor=progressor)
if show_progress:
    pbar.close()

# if show_progress:
#     pbar = tqdm.tqdm(total=len(db), ascii=True)
#     progressor = pbar.update
# else:
#     progressor = None
# saveStationLocalPlots(db, folder_name, progressor=progressor)
# if show_progress:
#     pbar.close()
