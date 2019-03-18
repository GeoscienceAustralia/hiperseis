import numpy
from functools import partial

# Try and use the faster Fourier transform functions from the pyfft module if
# available
try:
    import pyfftw

    pyfftw.interfaces.cache.enable()
    pyfftw.interfaces.cache.set_keepalive_time(100)
    from pyfftw.interfaces.scipy_fftpack import fftn as fftn
    from pyfftw.interfaces.scipy_fftpack import ifftn as ifftn

    from pyfftw.interfaces.numpy_fft import rfft as rfft
    from pyfftw.interfaces.numpy_fft import irfft as irfft

    #fftn  = partial(fftn_, planner_effort='FFTW_ESTIMATE')
    #ifftn = partial(ifftn_, planner_effort='FFTW_ESTIMATE')

    #rfft  = partial(rfft_, planner_effort='FFTW_ESTIMATE')
    #irfft = partial(irfft_, planner_effort='FFTW_ESTIMATE')
except ImportError:
    from scipy.fftpack import fftn, ifftn
    from numpy.fft import rfft, irfft
# end try

def ndflip(a):
    """Inverts an n-dimensional array along each of its axes"""
    ind = (slice(None,None,-1),)*a.ndim
    return a[ind]
# end func

