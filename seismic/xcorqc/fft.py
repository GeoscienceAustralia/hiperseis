import numpy as np

#TODO move below to fft lib
# Try and use the faster Fourier transform functions from the anfft module if
# available
try:
    import anfft as _anfft
    # measure == True for self-optimisation of repeat Fourier transforms of 
    # similarly-shaped arrays

    def fftn(A,shape=None):
        if shape != None:
            A = _checkffttype(A)
            A = procrustes(A,target=shape,side='after',padval=0)
        return _anfft.fftn(A,measure=True)

    def ifftn(B,shape=None):
        if shape != None:
            B = _checkffttype(B)
            B = procrustes(B,target=shape,side='after',padval=0)
        return _anfft.ifftn(B,measure=True)

    def _checkffttype(C):
        # make sure input arrays are typed correctly for FFTW
        if C.dtype == 'complex256':
            # the only incompatible complex type --> complex64
            C = np.complex128(C)
        elif C.dtype not in ['float32','float64','complex64','complex128']:
            # any other incompatible type --> float64
            C = np.float64(C)
        return C
		
# Otherwise use the normal scipy fftpack ones instead (~2-3x slower!)
except ImportError:
    """
    print \
    "Module 'anfft' (FFTW Python bindings) could not be imported.\n"\
    "To install it, try running 'easy_install anfft' from the terminal.\n"\
    "Falling back on the slower 'fftpack' module for ND Fourier transforms."
    """
    from scipy.fftpack import fftn, ifftn

def ndflip(a):
    """Inverts an n-dimensional array along each of its axes"""
    ind = (slice(None,None,-1),)*a.ndim
    return a[ind]

def procrustes(a,target,side='both',padval=0):
    """
    Forces an array to a target size by either padding it with a constant or
    truncating it

    Arguments:
        a     Input array of any type or shape
        target     Dimensions to pad/trim to, must be a list or tuple
    """

    try:
        if len(target) != a.ndim:
            raise TypeError('Target shape must have the same number of dimensions as the input')
    except TypeError:
        raise TypeError('Target must be array-like')

    try:
        b = np.ones(target,a.dtype)*padval
    except TypeError:
        raise TypeError('Pad value must be numeric')
    except ValueError:
        raise ValueError('Pad value must be scalar')

    aind = [slice(None,None)]*a.ndim
    bind = [slice(None,None)]*a.ndim

    # pad/trim comes after the array in each dimension
    if side == 'after':
        for dd in xrange(a.ndim):
            if a.shape[dd] > target[dd]:
                aind[dd] = slice(None,target[dd])
            elif a.shape[dd] < target[dd]:
                bind[dd] = slice(None,a.shape[dd])

    # pad/trim comes before the array in each dimension
    elif side == 'before':
        for dd in xrange(a.ndim):
            if a.shape[dd] > target[dd]:
                aind[dd] = slice(a.shape[dd]-target[dd],None)
            elif a.shape[dd] < target[dd]:
                bind[dd] = slice(target[dd]-a.shape[dd],None)

    # pad/trim both sides of the array in each dimension
    elif side == 'both':
        for dd in xrange(a.ndim):
            if a.shape[dd] > target[dd]:
                diff = (a.shape[dd]-target[dd])/2.
                aind[dd] = slice(np.floor(diff),a.shape[dd]-np.ceil(diff))
            elif a.shape[dd] < target[dd]:
                diff = (target[dd]-a.shape[dd])/2.
                bind[dd] = slice(np.floor(diff),target[dd]-np.ceil(diff))
    
    else:
        raise Exception('Invalid choice of pad type: %s' %side)

    b[bind] = a[aind]

    return b

