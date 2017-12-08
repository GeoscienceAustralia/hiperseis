import pyasdf
from pyasdf import ASDFDataSet
from obspy.core import Stream

class ASDFUtil(object):
    def __init__(self, asdfFilePath):
        self.asdfFile = asdfFilePath
        self.asdfDataSet = pyasdf.ASDFDataSet(self.asdfFile, mode='r')

    def get_waveforms(self, network, station, location, channel, starttime,
                      endtime, quality=None, minimumlength=None,
                      longestonly=None, filename=None, attach_response=False,
                      **kwargs):
        st = Stream()
        #ignoring channel for now as all the 7D network waveforms have only BH? channels
        filteredList = [i for i in self.asdfDataSet.waveforms[network+'.'+station].list() if 'raw_recording' in i and UTC(i.split("__")[1]) < starttime and UTC(i.split("__")[2]) > endtime]
        for t in filteredList:
            st += self.asdfDataSet.waveforms[network+'.'+station][t]
        return st
