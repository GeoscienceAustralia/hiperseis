#Parallelised autopick harvester. We have like a million picks so this is the only way to go about it
import sys

from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
from seismic.ml_classifier.data_harvester.autopicks import pickLoaderRand


from obspy.clients.fdsn.client import Client
ic=Client("IRIS")
fds=FederatedASDFDataSet('/g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt', variant='db', use_json_db=True, logger=None)

import numpy as np

pl=pickLoaderRand(fds,ic)

import multiprocessing as mp

nproc=mp.cpu_count()
print(nproc)
def lockInit(l):
    global lock
    lock=l
l=mp.Lock()
pool=mp.Pool(processes=nproc,initializer=lockInit,initargs=(l,))

def pickSave(ind):
    [pick,data]=pl.getPick(ind)
    if data:
        lock.acquire()
        np.save(data,saveDir+str(ind)+'.npy')
        dbtxt.write(' '.join(pick)+' '+str(ind)+'\n')
        lock.release()

pickCtr=0
saveDir='/g/data/ha3/rlt118/neural-datasets/autopicks/'
dbtxt=open(saveDir+'ensemble.s.txt','w')
pool.map(pickSave,range(0,len(pl)))
dbtxt.close()

#the waveform corresponding to the nth line of the pick file is saved as '${n}.npy'


