from obspy.core import stream
from pathlib import Path
import os

mseed_dir = '/g/data/ha3/rlt118/neural-datasets/autopicks-mseed/'

pathlist = Path(mseed_dir).glob('**/*.mseed')

for path in pathlist:
    s = stream.read(str(path))
    ppath = os.path.splitext(str(path))[0]+'.png'
    s.plot(outfile=ppath)
