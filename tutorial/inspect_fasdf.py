#! /usr/bin/env python
"""
Inspect a Federated ASDF database

CreationDate:   2021-05-17
Developer:      fei.zhang@ga.gov.au

"""

import os
import pandas as pd
import pyasdf
import sqlite3


def get_netsta_pdf(dburl="/g/data/ha3/Passive/SHARED_DATA/Index/778837537aa72d892df7b0ba22320f537c1d8f6a.db"):
    """ Read sqlite table for network.station info into a pandas DataFrame
    """

    con = sqlite3.connect(dburl)

    # Read sqlite query results into a pandas DataFrame
    df = pd.read_sql_query("SELECT * from netsta", con)

    # Verify the result of SQL query stored in the dataframe
    # print("The size of the TABLE netsta=", df.size)
    # print(df.head(5))

    con.close()

    return df


def select_netsta_rows(net, sta):
    """
    select rows according the given (net,sta)
    Args:
        net: network code
        sta: station code

    Returns: pandas dataframe

    """
    netsta_pdf = get_netsta_pdf()

    filter1 = netsta_pdf["net"] == net
    filter2 = netsta_pdf["sta"] == sta  # 'ARMA' #'AXCOZ'

    a_pdf = netsta_pdf.loc[filter1 & filter2]

    if (a_pdf.size >0):
        print(a_pdf.size, "rows found for ", net,sta)
    else:
        print("!!! WARNING !!!: NO record found for ",net,sta)

    return a_pdf


def check_h5file(asdffile):
    """
    Check a particular asdf file listed in the index.txt of the federated ASDF database
    :return:
    """
    # asdffile = "/g/data/ha3/Passive/STRIPED_DATA/GA_PERM/2018-2019.h5"  # very large
    # asdffile = "/g/data/ha3/GASeisDataArchive/DevSpace/2020.h5"

    ds = pyasdf.ASDFDataSet(asdffile, mode="r")

    ##  This may take a few hours to complete run in VDI
    # ds.validate()  # PyASDF provide an validation function

    ###  This takes a few hours to complete run in VDI
    for net_station in ds.waveforms.list():
        net, sta = net_station.split(".")
        print("********** network and station *************", net, sta)
        select_netsta_rows(net, sta)
        print(ds.waveforms[net_station])


# =============================================
# Section for quick test of this script
# ---------------------------------------------
if __name__ == "__main__":

    netcode = 'AU'
    stacode = 'ARMA'  # 'AXCOZ'

    print("********* Query for the networ, station (%s,%s)" % (netcode, stacode))
    df = select_netsta_rows(netcode, stacode)
    print(df.head(5))

    path2_asdffile = "/g/data/ha3/Passive/STRIPED_DATA/GA_PERM/2018-2019.h5"
    check_h5file(path2_asdffile)