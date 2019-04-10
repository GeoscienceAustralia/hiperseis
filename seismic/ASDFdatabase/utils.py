from obspy.core import UTCDateTime
import numpy as np

MAX_DATE = UTCDateTime(4102444800.0)
MIN_DATE = UTCDateTime(-2208988800.0)

def rtp2xyz(r, theta, phi):
    xout = np.zeros((r.shape[0], 3))
    rst = r * np.sin(theta)
    xout[:, 0] = rst * np.cos(phi)
    xout[:, 1] = rst * np.sin(phi)
    xout[:, 2] = r * np.cos(theta)
    return xout
# end func
