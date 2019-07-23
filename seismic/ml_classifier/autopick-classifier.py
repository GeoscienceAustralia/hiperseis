from .data_harvester.autopicks import pickLoader

from . import model
import numpy as np
from . import autopick_data

network=model.shakenet(pretrained_weights='shakenet-model.hdf5')

datagen=autopick_data.dataGenerator(128)

outfile='/g/data/ha3/rlt118/autopicks/ensemble.s.verified.txt'



with open(outfile,'w') as f:
    #initialise the file with the same format as the input    
    f.write("#eventID originTimestamp mag originLon originLat originDepthKm net sta cha pickTimestamp phase stationLon stationLat az baz distance ttResidual snr qualityMeasureCWT domFreq qualityMeasureSlope bandIndex nSigma\n")

    results=network.predict_generator(datagen)

print(len(results))

