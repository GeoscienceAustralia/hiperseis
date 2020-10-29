import glob
import os
import pickle

for fname in glob.glob('/g/data/ha3/rlt118/neural-datasets/categoriser-teleseismic/smallset/*.pkl'):
    with open(fname, "rb") as f:
        s = pickle.load(f)
    s.plot(outfile=os.path.splitext(fname)[0]+'.png')

