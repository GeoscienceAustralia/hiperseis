#script to collect waveforms corresponding to the automated S-wave picks and save them with an index to a directory
import sys
sys.path.append('/g/data1a/ha3/rakib/seismic/pst/passive-seismic')
from ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
from autopicks import pickLoader

from obspy.clients.fdsn.client import Client
ic=Client("IRIS")
fds=FederatedASDFDataSet('/g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt', variant='db', use_json_db=True, logger=None)
pickIt=pickLoader(fds,ic)

import numpy as np

pickCtr=0
saveDir='/g/data/ha3/rlt118/neural-datasets/autopicks/'
with open(saveDir+'ensemble.s.txt','w') as f:
    for pick,data in pickIt:
        pickCtr+=1
        f.write(' '.join(pick)+'\n')
        np.save(saveDir+str(pickCtr)+'.npy',data)
        

#the waveform corresponding to the nth line of the pick file is saved as '${n}.npy'
