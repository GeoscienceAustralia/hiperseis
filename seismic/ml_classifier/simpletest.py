#load the first test object and run it through the trained network to see what happens
from model import *
import numpy as np
import matplotlib.pyplot as plt

model=shakenet(pretrained_weights='shakenet-model.hdf5')

testin=np.load('/g/data/ha3/rlt118/neural-datasets/valset/771_trc.npy')
testin=np.reshape(testin,(1,1024,1))

testdists=np.load('/g/data/ha3/rlt118/neural-datasets/valset/771_msk.npy')

res=model.predict(testin,batch_size=1)
plt.plot(res[0,:,0])
plt.plot(res[0,:,1])
plt.plot(res[0,:,2])
plt.plot(testdists[:,0]/15)
plt.plot(testdists[:,1]/15)
plt.show()

testin=np.load('/g/data/ha3/rlt118/neural-datasets/valset/883_trc.npy')
testin=np.reshape(testin,(1,1024,1))

testdists=np.load('/g/data/ha3/rlt118/neural-datasets/valset/883_msk.npy')

res=model.predict(testin,batch_size=1)
plt.plot(res[0,:,0])
plt.plot(res[0,:,1])
plt.plot(res[0,:,2])
plt.plot(testdists[:,0]/15)
plt.plot(testdists[:,1]/15)
plt.show()
