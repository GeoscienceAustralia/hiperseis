import sys
import glob
import pickle
import os

dirname = sys.argv[1]

for fname in glob.glob(dirname.rstrip('/')+'/*.pkl'):
    with open(fname, "rb") as f:
        s = pickle.load(f)
    basename = os.path.splitext(fname)[0]
    s.write(basename+".mseed")
