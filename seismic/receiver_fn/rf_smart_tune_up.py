import numpy as np
from scipy.spatial.distance import euclidean
from fastdtw import fastdtw
import itertools as iter
from matplotlib.pyplot import plot, show, figure, ylim, xlabel, ylabel, legend, subplot2grid, GridSpec
from copy import deepcopy
from sklearn.cluster import DBSCAN

# This program shows different techniques of selecting good RF functions


def compare_pairs(data):
    distance, _ = fastdtw(data[0], data[1], dist=euclidean)
    return distance

def crossSpectrum(x, y):

#-------------------Remove mean-------------------
    nperseg = x.size/20
    cross = np.zeros(nperseg, dtype='complex128')
    for ind in range(x.size / nperseg):

        xp = x[ind * nperseg: (ind + 1)*nperseg]
        yp = y[ind * nperseg: (ind + 1)*nperseg]
        xp = xp - np.mean(xp)
        yp = yp - np.mean(xp)

        # Do FFT
        cfx = np.fft.fft(xp)
        cfy = np.fft.fft(yp)

        # Get cross spectrum
        cross += cfx.conj() * cfy

    freq = np.fft.fftfreq(nperseg)
    return cross, freq


def coh(y,y2):
    p11, freq = crossSpectrum(y, y)
    p22, freq = crossSpectrum(y2, y2)
    p12, freq = crossSpectrum(y, y2)
    # coherence
    part1 = np.divide(np.abs(p12)**2, p11.real, out=np.zeros_like(np.abs(p12)**2), where=p11.real!=0)
    coh = np.divide(part1, p22.real, out=np.zeros_like(part1), where=p22.real!=0)

#   plot( freq[freq > 0], coh[freq > 0])
#   show()
#   return coh[freq > 0]

    return  freq[freq > 0], coh[freq > 0]


# -------------Main---------------------------------

if __name__ == '__main__':
    t=[]
    x=[]
    with open('rf_test.dat') as fin:
         for line in fin:
             a, b = line.split()
             t.append(a)
             x.append(b)

    t = np.array(t, dtype='float')
    x = np.array(x, dtype='float')
    swipe = []

    # Generating synthetic data and adding noise

    for i in xrange(10):
        noise = np.random.random(x.size) - 0.5
        if i < 4:
            swipe.append(x + noise/(i+1))
        else:
            swipe.append(x)

    swipe = np.array(swipe)

    # First trying to find median value

    average = np.median(swipe,axis=0)

    # Then see how coherent each signal is comparing to median

    fig1=figure(1)
    gs1=GridSpec(5,6)
    gs1.update(left=0.01,right=0.01)
    ax1=subplot2grid((5,6),(0,0),rowspan=4)
    for i in xrange(swipe.shape[0]):

         ax1.plot(t,swipe[i,:]+i)

    ax1.set_ylabel('RF number')

    ax2=subplot2grid((5,6),(4,0),colspan=1,rowspan=2)
    ax2.plot(t,average)
    ax2.set_title('Median RF')

    ax3=subplot2grid((5,6),(0,2),colspan=4,rowspan=3)

    for i in xrange(swipe.shape[0]):
          ax3.plot(coh(average,swipe[i,:])[0],coh(average,swipe[i,:])[1],label=str(i))
    ax3.legend()
    ax3.set_xlabel('Apparent frequency')
    ax3.set_ylabel('Coherence with median')


    # Next method is remove one trace and se its contribution to average, if not changed then is similar to main group
    ax4=subplot2grid((5,6),(3,2),colspan=2,rowspan=2)
    ind=np.ones((swipe.shape[0],),bool)
    dev=[]
    ch=[]
    knive=[]
    sn=[]
    pulse_ind=np.max(t[t<0])-1.
    for i in xrange(swipe.shape[0]):
        ch.append(np.amax(coh(average,swipe[i,:])))
        dev.append(np.sum((swipe[i,:]-average)**2)/(swipe.shape[0]-1))
        ind[i]=False
        knive.append(np.std(swipe[ind]))
        sn.append(np.std(swipe[i,t>pulse_ind])/np.std(swipe[i,t<pulse_ind]))

    ax4.plot(list(range(len(knive))),knive)
    ax4.set_xlabel('Number of Excluded RF')
    ax4.set_ylabel('RMS')

    ax6=subplot2grid((5,6),(0,1),rowspan=4)

    ax6.plot(sn,list(range(len(sn))))
    ax6.set_xlabel('S/N')


    # Another method is to find shape similarity between each RF

    ax5=subplot2grid((5,6),(3,4),colspan=2,rowspan=2)

    distance=map(compare_pairs,iter.combinations(swipe,2))
    index=list((i,j) for ((i,_),(j,_)) in iter.combinations(enumerate(swipe),2))
#   for i in xrange(len(index)):
#         print index[i],distance[i]
    # First check that distance betwen points
    index=np.array(index)
    distance=np.array(distance)
    matrix=np.zeros((np.amax(index)+1 ,1+np.amax(index)))+np.amax(distance)
#   print matrix[index].shape,distance.shape,index.shape
    matrix[index[:,0],index[:,1]]=distance[:]
    clustering=DBSCAN(eps=3,min_samples=2,metric='precomputed').fit(matrix)

    ax5.plot(list(range(len(list(clustering.labels_)))),list(clustering.labels_))
    ax5.set_xlabel('RF number')
    ax5.set_ylabel('Similarity group')
#   fig1.tight_layout()
    fig1.subplots_adjust(hspace=3.5,wspace=1.4)
    show()
