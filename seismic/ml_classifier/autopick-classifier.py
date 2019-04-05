import sys
sys.path.append('/g/data1a/ha3/rakib/seismic/pst/passive-seismic')
from ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
from data_harvester.autopicks import pickLoader

from obspy.clients.fdsn.client import Client

import model
import numpy as np

ic=Client("IRIS")
fds=FederatedASDFDataSet('/g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt', variant='db', use_json_db=True, logger=None)

pickIt=pickLoader(fds,ic)

network=model.shakenet(pretrained_weights='shakenet-model.hdf5')

outfile='/g/data/ha3/rlt118/curatedpicks/ensemble.s.txt'
with open(outfile,'w') as f:
    #initialise the file with the same format as the input    
    f.write("#eventID originTimestamp mag originLon originLat originDepthKm net sta cha pickTimestamp phase stationLon stationLat az baz distance ttResidual snr qualityMeasureCWT domFreq qualityMeasureSlope bandIndex nSigma\n")

    for pick,data in pickIt:
        inArr=np.resize(data,(1,6002,1))
        print model.predict(inArr,batch_size=1)
