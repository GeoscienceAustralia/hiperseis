from datagenerator import DataGenerator
import os
from random import shuffle,seed
#this is hardcoded, maybe make a more general version that gets the IDs from the folder in the future
datafolder='/g/data/ha3/rlt118/neural-datasets/categoriser-teleseismic/smallset/'
#datafolder='/home/ubuntu/seismic_waves_4_ml'  #aws test sample data dir

#build a list of IDs and dictionary of labels
files=os.listdir(datafolder)
IDs=[]
for fname in files:
    if fname.endswith('.npy'):
        IDs.append(fname.rstrip('.npy'))

labels={}
Sctr=0
Nctr=0
for ID in IDs:
    if ID.endswith('_S'):
        labels[ID]=0
        Sctr+=1
    else:
        labels[ID]=1
        Nctr+=1

seed(0)
shuffle(IDs)
trainLen=int(0.8*len(IDs))#use 80% of the data as training set
valLen=int(0.1*len(IDs))
partition={}
partition['train'],partition['val'],partition['test']=IDs[0:trainLen],IDs[trainLen:trainLen+valLen],IDs[trainLen+valLen:]


#create ID lists for S, N wave test cases

inv_labels=[[],[]]
for ID in partition['test']:
    inv_labels[labels[ID]].append(ID)

print('dataset contains '+str(Sctr)+' S waves and '+str(Nctr)+' noise waveforms.')
print('training set contains '+str(len(partition['train']))+' waveforms and test set contains '+str(len(partition['test']))+' waveforms.')



def trainGenerator(batch_size):
    

    return DataGenerator(partition['train'],labels,datafolder,batch_size=batch_size)

def valGenerator(batch_size):
    return DataGenerator(partition['val'],labels,datafolder,batch_size=batch_size)

def testGenerator(batch_size):
    

    return DataGenerator(partition['test'],labels,datafolder,batch_size=batch_size)


def SGenerator(batch_size):
    return DataGenerator(inv_labels[0],labels,datafolder,batch_size=batch_size)

def NGenerator(batch_size):
    return DataGenerator(inv_labels[1],labels,datafolder,batch_size=batch_size)

def getIDs(category):
    return inv_labels[category]
