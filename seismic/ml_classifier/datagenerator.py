import numpy as np
import keras
import os.path as path

#inspired by tutorial at https://stanford.edu/~shervine/blog/keras-how-to-generate-data-on-the-fly

class DataGenerator(keras.utils.Sequence):
    def __init__(self,list_IDs,labels,dataDir,batch_size=32,dim=(6002,),n_channels=1,n_classes=2,shuffle=True):
        #initialise data generator
        self.dim=dim
        self.batch_size=batch_size
        self.labels=labels
        self.list_IDs=list_IDs
        self.n_channels=n_channels
        self.n_classes=n_classes
        self.shuffle=shuffle
        self.dataDir=dataDir
        self.on_epoch_end()

    def __len__(self):
        #number of batches per epoch. Ensure batch_size is int
        return int(len(self.list_IDs)/self.batch_size)

    def __getitem__(self,index):
        #indices for this batch        
        indexes=self.indexes[index*self.batch_size:(index+1)*self.batch_size]
        #IDs for this batch
        list_IDs_temp=[self.list_IDs[k] for k in indexes]

        #get the batch arrays
        X,y=self.__data_generation(list_IDs_temp)

        return X,y

    def on_epoch_end(self):
        #reshuffle (or just reset) indices for the next batch
        self.indexes=np.arange(len(self.list_IDs))
        if self.shuffle:
            np.random.shuffle(self.indexes)

    def __data_generation(self,list_IDs_temp):
        #generate empty batch arrays. X is the input, y is the mask
        X=np.empty(self.__dimTup())
        y=np.empty(self.batch_size,dtype=int)
        #Generate the data
        for i,ID in enumerate(list_IDs_temp):
            #Store input
            X[i,]=np.resize(np.load(path.join(self.dataDir,ID+'.npy')),self.__dimTup()[1:])#some of the resampled traces don't have length of exactly 6002.
                                                                                      #np.resize() will truncate or pad them through repetition as necessary.
            y[i]=self.labels[ID]
        return X,keras.utils.to_categorical(y,num_classes=self.n_classes)

    def __dimTup(self):
        return (self.batch_size,)+self.dim+(self.n_channels,)
        
