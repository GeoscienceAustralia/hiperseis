import json
from obspy import read_inventory, read_events, UTCDateTime, Stream, read
from obspy.clients.fdsn.client import Client
from obspy.clients.fdsn.header import FDSNException
import pyasdf
import itertools
import sys
import os


class Client2ASDF(object):
    def __init__(self, client="IRIS", network="AU"):
        """Object to query IRIS or other client with a bounding box and time interval.
           Returns ASDF containing station information and data for interval."""
        self._client = client
        self._network = network
        self.ref_stations = []

    def queryByBBoxInterval(self, outputFileName, bbox, timeinterval, chan='*Z', bbpadding=2,
                            event_id=None, verbose=False):
        """ Time interval is a tuple (starttime,endtime)
        """
        assert len(timeinterval) == 2, "timeinterval must be a tuple of ascending timestamps. len=" + str(
            len(timeinterval)) + " " + str(timeinterval)

        query_ds = pyasdf.ASDFDataSet(outputFileName)

        client = Client(self._client)
        ref_inv = client.get_stations(network=self._network,
                                      starttime=UTCDateTime(timeinterval[0]),
                                      endtime=UTCDateTime(timeinterval[1]),
                                      minlongitude=bbox[0] - bbpadding,
                                      maxlongitude=bbox[1] + bbpadding,
                                      minlatitude=bbox[2] - bbpadding,
                                      maxlatitude=bbox[3] + bbpadding,
                                      level='channel')

        if verbose:
            print(ref_inv)

        ref_st = Stream()

        # go through inventory and request timeseries data
        for net in ref_inv:
            for stn in net:
                stime = UTCDateTime(timeinterval[0])
                etime = UTCDateTime(timeinterval[1])
                step = 3600*24*10
                while stime + step < etime:
                    try:
                        ref_st = client.get_waveforms(network=net.code, station=stn.code,
                                                      channel=chan, location='*',
                                                      starttime=stime,
                                                      endtime=stime+step)
                        print ref_st
                        self.ref_stations.append(net.code + '.' + stn.code)
                        st_inv = ref_inv.select(station=stn.code, channel=chan)
                        
                        query_ds.add_stationxml(st_inv)
                        for tr in ref_st:
                            query_ds.add_waveforms(tr, "reference_station")
                    except FDSNException:
                        print('Data not available from Reference Station: ' + stn.code)
                    # end try
                    stime += step
                #wend
        # end for

        #tr.write(os.path.join(os.path.dirname(outputFileName), tr.id + ".MSEED"),
        #         format="MSEED") # Don't write miniseed
        if verbose:
            print("Wrote Reference Waveforms to ASDF file: " + outputFileName)
            print('\nWaveform data query completed.')

        metaOutputFileName = os.path.join(os.path.dirname(outputFileName),
                                          'meta.%s.xml'%(os.path.basename(outputFileName)))
        ref_inv.write(metaOutputFileName, format="STATIONXML")
        del query_ds

if __name__ == "__main__":
    fn = '/g/data/ha3/rakib/_ANU/7G(2013-2015)/refData/stka.6m.h5'
    fn = '/tmp/stka.6m.h5'
    c = Client2ASDF()

    p = 0.0001
    #inka = (-27.741, 140.746)
    #c.queryByBBoxInterval(fn, [inka[1] - p, inka[1] + p, inka[0] - p, inka[0] + p],
    #                      ("2015-01-01T00:00:00", "2015-06-01T00:00:00"), 1, verbose=True)

    stka = (-31.8769, 141.5952)
    c.queryByBBoxInterval(fn, [stka[1] - p, stka[1] + p, stka[0] - p, stka[0] + p],
                          ("2015-01-01T00:00:00", "2015-01-02T00:00:00"), 1, verbose=True)

    # test that the file exists and was written
    assert os.path.isfile(fn), "ASDF file not written"
    #os.remove(fn)
