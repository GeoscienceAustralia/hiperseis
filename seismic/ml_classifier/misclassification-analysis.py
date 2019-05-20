from .model import *
from .data import *
import os.path as path
import matplotlib
matplotlib.use('Agg')
import obspy.core as oc


model=shakenet(pretrained_weights='shakenet-model.hdf5')


dataDir='/g/data/ha3/rlt118/neural-datasets/categoriser-teleseismic/'
SIDs=getIDs(0)
NIDs=getIDs(1)

for ID in SIDs:
    x=np.resize(np.load(path.join(dataDir,ID+'.npy')),(1,6002,1))
    result=model.predict(x,batch_size=1)
    trc=oc.stream.read(path.join(dataDir,ID+'.pkl'))
    if np.argmax(result[0])==0:
        outfname='correctS/'+ID+'.png'
    else:
        outfname='badS-N/'+ID+'.png'
    trc.plot(outfile=outfname)

for ID in NIDs:
    x=np.resize(np.load(path.join(dataDir,ID+'.npy')),(1,6002,1))
    result=model.predict(x,batch_size=1)
    trc=oc.stream.read(path.join(dataDir,ID+'.pkl'))
    if np.argmax(result[0])==0:
        outfname='badN-S/'+ID+'.png'
    else:
        outfname='correctN/'+ID+'.png'
    trc.plot(outfile=outfname)
