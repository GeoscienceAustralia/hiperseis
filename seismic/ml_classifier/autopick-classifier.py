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
rawoutfile='/g/data/ha3/rlt118/curatedpicks/ensemble.s.avail.txt'
poutfile='/g/data/ha3/rlt118/curatedpicks/ensemble.s.badp.txt'
with open(outfile,'w') as f, open(rawoutfile,'w') as g, open(poutfile,'w') as e:
    #initialise the file with the same format as the input    
    f.write("#eventID originTimestamp mag originLon originLat originDepthKm net sta cha pickTimestamp phase stationLon stationLat az baz distance ttResidual snr qualityMeasureCWT domFreq qualityMeasureSlope bandIndex nSigma\n")
    g.write("#eventID originTimestamp mag originLon originLat originDepthKm net sta cha pickTimestamp phase stationLon stationLat az baz distance ttResidual snr qualityMeasureCWT domFreq qualityMeasureSlope bandIndex nSigma\n")
    for pick,data in pickIt:
        print pick[15]+" degrees"
        inArr=np.resize(data,(1,6002,1))
        result= network.predict(inArr,batch_size=1)
        g.write(' '.join(pick)+'\n')
        if np.argmax(result[0])==0:
            f.write(' '.join(pick)+'\n')
            print "Is an S wave"
        elif np.argmax(result[0])==1:
            print 'false P wave'
            e.write(' '.join(pick)+'\n')
        else:
            print "Possibly not an S wave"
