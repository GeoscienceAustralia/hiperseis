#class for storing processed seismic traces with ground-truth probability distributions for P and S picks.

from obspy.core import *
import numpy as np

class LearnObj:
    def __init__(self,trace,pdist,sdist):
        if len(trace)!=len(pdist) or len(pdist)!=len(sdist): raise ValueError('Cannot create a learning point with mismatching data lengths!')
        if type(trace) is not Trace or type(pdist) is not np.ndarray or type(sdist) is not np.ndarray: raise TypeError('Learning point needs trace as ObsPy trace object, pdist and sdist as numpy ndarrays')
        self.tr=trace
        self.pd=pdist
        self.sd=sdist

    def __len__(self):
        return len(self.tr)

    def trvals(self):#ndarray for training
        return self.tr.data
    
    def getTrace(self): #return trace including all metadata
        return self.tr

    def groundTruth(self): #return the probability distributions for checking against output of picker
        return [self.pd,self.sd]
