from datagenerator import DataGenerator
import os
from random import shuffle,seed
#this is hardcoded, maybe make a more general version that gets the IDs from the folder in the future
datafolder='/g/data/ha3/rlt118/neural-datasets/autopicks/'

#build a list of IDs and dictionary of labels
with open(datafolder+'ensemble.s.txt','r') as f:
    pickdb=f.readlines()

IDs=range(1,len(pickdb)+1)

labels={}
for ID in IDs:
    labels[ID]=0

def autoGenerator(batch_size):
    

    return DataGenerator(IDs,labels,datafolder,batch_size=batch_size)

