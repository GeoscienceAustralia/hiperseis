#! /usr/bin/env python
"""
Description:
    Retrieve waveform data from a federated ASDF database
    

CreationDate:   2021-02-09
Developer:      fei.zhang@ga.gov.au

"""

import os
import sys
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet


# Section to define functions or class
def get_wdata(indexfile, net, sta, loc, cha, start, end):
    """
    define a function
    :return:
    """

    fds = FederatedASDFDataSet(indexfile)

    d = fds.get_waveforms(net, sta, loc, cha, start, end)

    return d


# =============================================
# CLI entry point  of this script

# Prepare run environment:
# source vdi_env.sh
# export PYTHONPATH=/g/data/ha3/fxz547/Githubz/hiperseis
# ---------------------------------------------
if __name__ == "__main__":

    if len(sys.argv) < 2:
        index_txt = "/g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt"
        print("! WARN:  No index text file provided, will use the default ", index_txt)
        print("******** General USAGE: python3 %s %s **********" % (sys.argv[0], "index_txt_file"))
    else:
        index_txt = sys.argv[1]  # production data = "/g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt"

    # search criteria for the waveform data (cannot use wildcard *, sql statement)
    # A good example, which will get waveform data
    netcode = 'AU'
    stacode = 'QIS'
    location = ''
    channel = 'BHZ'
    start_dt = '2015-06-01T00:00:00'
    end_dt = '2015-06-02T00:06:00'

# ********** network and station ************* AU SDAN,
    # Number of rows found in table netsta: 0
    # Contents of the data set for station AU.SDAN:
    #     - Has no StationXML file
    #     - 1 Waveform Tag(s):
    #         raw_recording
# inconsistent ASDF vs FedASDF
# have waveform data in "/g/data/ha3/Passive/STRIPED_DATA/GA_PERM/2018-2019.h5" , but cannot found in fedasdf.
    netcode = 'AU'
    stacode = 'SDAN'
    location = '00'
    channel = 'BHZ'
    start_dt = '2018-01-01T04:00:00'
    end_dt = '2019-12-17T23:59:59'

# ********** network and station ************* AU STKA
    # Number of rows found 35
    # Contents of the data set for station AU.STKA:
    #     - Has a StationXML file
    #     - 1 Waveform Tag(s):
    #         raw_recording
# consistent with the ASDF direct access
#Found database: /g/data/ha3/Passive/SHARED_DATA/Index/778837537aa72d892df7b0ba22320f537c1d8f6a.db
# 1 Trace(s) in Stream:
# AU.STKA..BHE | 2018-01-01T00:00:00.000000Z - 2018-01-01T23:59:59.975000Z | 40.0 Hz, 3456000 samples
    netcode = 'AU'
    stacode = 'STKA'
    location = ''
    channel = 'BHE'
    start_dt = '2018-01-01T00:00:00'
    end_dt = '2018-01-01T23:59:59.999'


    fds = FederatedASDFDataSet(index_txt)
    stream = fds.get_waveforms(netcode, stacode, location, channel, start_dt, end_dt)

    # process the data stream

    print(type(stream))  # <class 'obspy.core.stream.Stream'>

    print(stream)

    # select * from wdb where net='AU' and sta='QIS' and loc='' and cha='BHZ' and et>=1433116800.000000 and st<=1433203560.000000
    # 2 Trace(s) in Stream:
    # AU.QIS..BHZ | 2015-05-31T23:59:59.994500Z - 2015-06-01T23:59:59.994500Z | 40.0 Hz, 3456001 samples
    # AU.QIS..BHZ | 2015-06-01T23:59:55.769500Z - 2015-06-02T00:05:59.994500Z | 40.0 Hz, 14570 samples

    if (stream.count() > 0):
        # stream.plot()
        atrace = stream[0]
        atrace.plot()  # the first Trace
