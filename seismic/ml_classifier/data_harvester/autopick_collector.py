# script to collect waveforms corresponding to the automated S-wave picks and save them with an index to a directory
import sys
# sys.path.append('/g/data1a/ha3/rakib/seismic/pst/passive-seismic')
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
from seismic.ml_classifier.data_harvester.autopicks_approx import pickLoader
from obspy.core import UTCDateTime

fds = FederatedASDFDataSet('/g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt', variant='db', use_json_db=True,
                           logger=None)

import numpy as np

saveDir = '/g/data/ha3/rlt118/neural-datasets/autopicks/'
with open(saveDir + 'ensemble.s.txt', 'r') as f:
    picks = f.readlines()
    pickCtr = len(picks) - 1
    lastpick = picks[-2].split()
    mintime = UTCDateTime(float(lastpick[9]))

pickIt = pickLoader(fds, mintime=mintime)

with open(saveDir + 'ensemble.s.txt', 'w') as f:
    # write the existing picks
    f.writelines([item for item in picks[:-1]])

    for pick, data in pickIt:
        pickCtr += 1
        f.write(' '.join(pick) + '\n')
        np.save(saveDir + str(pickCtr) + '.npy', data)

# the waveform corresponding to the nth line of the pick file is saved as '${n}.npy'
