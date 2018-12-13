#!/usr/bin/env python
import os
import sys
import numpy as np
from obspy.core import utcdatetime
from obspy.geodetics.base import locations2degrees
from obspy.clients.fdsn.client import Client
# client = Client("IRIS")
from obspy import read_inventory
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
from obspy.clients.nrl import NRL
from obspy.core.util.obspy_types import (ComplexWithUncertainties,
                                         FloatWithUncertainties,
                                         FloatWithUncertaintiesAndUnit,
                                         ObsPyException,
                                         ZeroSamplingRate)
from collections import defaultdict

import pathlib2 as pathlib

def read_eng(fname):
    ''' We read Engdahl file in this subroutine

     AAI   Ambon             BMG, Indonesia, IA-Ne              -3.6870  128.1945      0.0   2005001  2286324  I
     AAII                                                       -3.6871  128.1940      0.0   2005001  2286324  I
     AAK   Ala Archa         Kyrgyzstan                         42.6390   74.4940      0.0   2005001  2286324  I
     ABJI                                                       -7.7957  114.2342      0.0   2005001  2286324  I
     APSI                                                       -0.9108  121.6487      0.0   2005001  2286324  I
     AS01  Alice Springs Arra                                  -23.6647  133.9508      0.0   2005001  2286324  I
    0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
              10        20        30        40        50        60        70        80
    '''
    eng_stn = []
    columns = ((0, 6), (59, 67), (68, 77), (78, 86))
    with open(fname) as stations:
        for string in stations:
            dataline = [string[c[0]:c[1]].strip() for c in columns]
            eng_stn.append(dataline)
    eng_stn = np.array(eng_stn)
    return eng_stn


def read_isc(fname):
    # Engdahl file is imported and stored, now lets find appropriate network codes based on coordinates
    # first we read catalogue supplied by ISC that inherited Engdahl work
    '''
    109C     32.8892 -117.1100     0.0 2006-06-01 04:11:18 2008-01-04 01:26:30
    109C     32.8882 -117.1050   150.0 2008-01-04 01:26:30
                 FDSN 109C   TA -- BHZ 2004-05-04 23:00:00 2005-03-03 23:59:59
                 FDSN 109C   TA -- LHZ 2004-05-04 23:00:00 2005-03-03 23:59:59
                 FDSN 109C   TA -- BHZ 2005-04-11 00:00:00 2006-01-25 22:31:10
                 FDSN 109C   TA -- LHZ 2005-04-11 00:00:00 2006-01-25 22:31:10
    0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
              10        20        30        40        50        60        70        80

    '''
    # Target format for each row of isc_inv is:
    #   [STA, NET, LAT, LONG, ELE, CHA, START, END]
    isc_inv = []
    with open(fname) as stations:
        for string in stations:
            if len(string[0:5].strip()) > 2:
                # Station code is valid. Treat this as the beginning of a new station record.
                header = True
                #           STA     LAT      LON       ELE       START      END
                columns = ((0, 5), (7, 16), (17, 26), (27, 34), (35, 54,), (55, 74))
                dataline = [string[c[0]:c[1]] for c in columns]
                lon = (dataline[2])
                lat = (dataline[1])
                ele = (dataline[3])
                stn = dataline[0].strip()
                dataline.insert(1, 'IR')    # default Network Code
                dataline.insert(5, 'SHZ')   # default Channel Code
                clean = [x.strip(' ') for x in dataline]
                isc_inv.append(clean)
            else:
                # No valid station code. This row is interpreted as channel information.
                #           STA       NET       CHA       START     END
                columns = ((17, 22), (23, 27), (30, 34), (35, 54), (55, 74))
                dataline = [string[c[0]:c[1]] for c in columns]
                if header:
                    # Replace header default Network Code with actual Network Code
                    isc_inv[-1][1] = dataline[1].strip()
                    # Replace default Channel Code with actual Channel Code
                    isc_inv[-1][5] = dataline[2].strip()
                    header = False
                # if dataline[0].strip() != stn:
                #     print dataline, " Station header mismatch"
                # Insert (lat, long, elevation) data for this stage
                dataline.insert(2, ele)
                dataline.insert(2, lon)
                dataline.insert(2, lat)
                clean = [x.strip(' ') for x in dataline]
                isc_inv.append(clean)

    isc_inv = np.array(isc_inv)
    stations.close()
    return isc_inv


def extract_unique_sensors_responses(inv):
    resultSensors = defaultdict(list)
    resultResponses = defaultdict(list)
    for net in inv.networks:
        for sta in net.stations:
            for cha in sta.channels:
                assert cha.code != 0
                try:
                    if (cha.code not in resultResponses.keys()):
                        if (cha.response and cha.response.get_paz()):

                            good = True
                            for rs in cha.response.response_stages:
                                if rs.decimation_delay is None:
                                    rs.decimation_delay = FloatWithUncertaintiesAndUnit(0)
                                if rs.decimation_correction is None:
                                    rs.decimation_correction = FloatWithUncertaintiesAndUnit(0)

                            if (good): # CONFUSED: why don't we index dict by network and station code as well, since cha.code is not unique
                                resultSensors[cha.code] = cha.sensor
                                resultResponses[cha.code] = cha.response
                except:
                    pass
    return resultSensors, resultResponses


# end func

def main(argv):
    with open("IRIS-ALL.xml", 'r', buffering=1024*1024) as f:
        inv = read_inventory(f)
    # if os.path.exists("IRIS-ALL.pkl"): # doesn't work on CentOS for some reason
    #     with open('IRIS-ALL.pkl', 'rb') as f:
    #         import cPickle as pkl
    #         inv = pkl.load(f)
    # else:
    #     inv = read_inventory("IRIS-ALL.xml")
    #     with open('IRIS-ALL.pkl', 'wb') as f:
    #         import pickle as pkl
    #         pkl.dump(inv, f, pkl.HIGHEST_PROTOCOL)
    sensorDict, responseDict = extract_unique_sensors_responses(inv)
    print('\nFound {0} response objects with keys: {1}'.format(len(responseDict.keys()), responseDict.keys()))

    # unknown stations in Indonesia are usually installed by Potsdam and we assume they have network name GE
    default_net = 'GE'
    ehb1 = read_eng('BMG.STN')
    ehb2 = read_eng('ISC.STN')
    ehb = np.unique(np.vstack((ehb1, ehb2)), axis=0)

    isc1 = read_isc('ehb.stn')
    isc2 = read_isc('iscehb.stn')
    isc = np.unique(np.vstack((isc1, isc2)), axis=0)

    catalogue = []
    for i in xrange(ehb.shape[0]):
        filed = False
        xml = False
        stn_found = isc[isc[:, 0] == ehb[i, 0], :]
        min_dist = 10e10
        if stn_found.shape[0] > 0:
            if stn_found.shape[0] > 1:
                for j in xrange(stn_found.shape[0]):
                    dist = locations2degrees(np.float(stn_found[j, 2]), np.float(stn_found[j, 3]), np.float(ehb[i, 1]),
                                             np.float(ehb[i, 2]))
                    if dist < min_dist:
                        min_dist = dist
                        record = stn_found[j, :]
            else:
                min_dist = locations2degrees(np.float(stn_found[0, 2]), np.float(stn_found[0, 3]), np.float(ehb[i, 1]),
                                             np.float(ehb[i, 2]))
                record = stn_found[0, :]

            #                Now we try to find the same station in XML file
            #                if min_dist > 1. or stn_found.shape[0]==0:

        xstn_found = inv.select(station=ehb[i, 0], channel="*HZ")

        if len(stn_found) == 0 and len(xstn_found) == 0:
            # we filed to find station anywhere and assign dummy values
            record = [ehb[i, 0], default_net, ehb[i, 1], ehb[i, 2], ehb[i, 3], 'SHZ', '1964-1-1 00:00:00',
                      '2599-12-31 23:59:59']
            min_dist = 0.
            filed = True
        else:
            # if station is found somehwere we try to iterate and see if XML has data giving it preference through adding extra value to min_dist found in ISC
            if len(xstn_found) > 0:
                #                        print "----------",len(xstn_found)
                #                        print xstn_found[0][0].latitude
                min_dist = min_dist + 0.1
                for j in xrange(len(xstn_found)):
                    dist = locations2degrees(xstn_found[j][0].latitude, xstn_found[j][0].longitude, np.float(ehb[i, 1]),
                                             np.float(ehb[i, 2]))
                    if min_dist > dist:
                        min_dist = dist
                        record = xstn_found[j]
                        #                                print record
                        xml = True

                    # last defence if stations have been done but distance between declared and found locations are more than 1 degree
        if min_dist > 1:
            record = [ehb[i, 0], default_net, ehb[i, 1], ehb[i, 2], ehb[i, 3], 'SHZ', '1964-1-1 00:00:00',
                      '2599-12-31 23:59:59']
            filed = True
        if xml:
            xml = False

        else:
            if filed:

                if len(record[7]) < 5:
                    record[7] = '2599-12-31 23:59:59'
                catalogue.append(record)

            else:

                stn_found = isc[(isc[:, 0] == record[0]) & (isc[:, 1] == record[1]), :]

                for k in xrange(stn_found.shape[0]):
                    net = Network(code=stn_found[k, 1], stations=[], description=' ')
                    if len(stn_found[k, 7]) < 5:
                        stn_found[k, 7] = '2599-12-31 23:59:59'
                    catalogue.append(stn_found[k, :])

    stn_found = np.unique(np.array(catalogue), axis=0)
    if len(stn_found[stn_found == '']) > 0 or len(stn_found[stn_found == ' ']) > 0:
        print
        "Some elements are empty, check the list"

    # we composed our inventory. However some stations from ISC list can be left behind. We check if some stations in ISC are forgotten
    lost = []
    for j in xrange(isc.shape[0]):
        # is there any common station name?
        common_st = stn_found[isc[j, 0] == stn_found[:, 0]]
        if common_st.shape[0] > 0:
            # is network code the same?
            common_net = common_st[common_st[:, 1] == isc[j, 1]]
            if common_net.shape[0] < 1:
                # ok we found forgotten one, check the XML
                if len(inv.select(station=isc[j, 0], network=isc[j, 1])) <= 0:
                    # Bingo...
                    lost.append(isc[j, :])
        else:
            if len(inv.select(station=isc[j, 0], network=isc[j, 1])) <= 0:
                # Bingo...
                lost.append(isc[j, :])

    stn_found = np.vstack((stn_found, np.array(lost)))

    netDict = defaultdict(list)
    for k in xrange(stn_found.shape[0]):
        result = inv.select(network=stn_found[k, 1])
        if (len(result.networks)):
            net = result.networks[0]
            net.stations = []
        else:
            net = Network(code=stn_found[k, 1], stations=[], description=' ')

        # print stn_found[k, 1]

        if len(stn_found[k, 7]) < 5:
            stn_found[k, 7] = '2599-12-31 23:59:59'
        catalogue.append(stn_found[k, :])
        sta = Station(code=stn_found[k, 0], creation_date=utcdatetime.UTCDateTime(stn_found[k, 6]), \
                      termination_date=utcdatetime.UTCDateTime(stn_found[k, 7]), \
                      site=Site(name=' '), \
                      latitude=np.float(stn_found[k, 2]), \
                      longitude=np.float(stn_found[k, 3]), \
                      elevation=np.float(stn_found[k, 4]))

        if (stn_found[k, 5] in responseDict.keys()):
            r = responseDict[stn_found[k, 5]]

            cha = Channel(code=stn_found[k, 5], \
                          depth=0., \
                          azimuth=0., \
                          dip=-90., \
                          location_code='', \
                          latitude=np.float(stn_found[k, 2]), \
                          longitude=np.float(stn_found[k, 3]), \
                          elevation=np.float(stn_found[k, 4]), \
                          # sensor=sensorDict[stn_found[k,5]], \
                          response=r)

            sta.channels.append(cha)

            if (type(netDict[stn_found[k, 1]]) == Network):
                netDict[stn_found[k, 1]].stations.append(sta)
            else:
                net.stations.append(sta)
                netDict[stn_found[k, 1]] = net

            #                 print 'np',stn_found[k,:]
            # end if

    our_xml = Inventory(networks=netDict.values(), source='EHB')

    print 'Writing output files..'
    pathlib.Path("output").mkdir(exist_ok=True)
    for inet, net in enumerate(our_xml.networks):
        currInv = Inventory(networks=[net], source='EHB')
        fname = "station.%d.xml" % (inet)
        try:
            currInv.write(os.path.join("output", fname), format="stationxml", validate=True)
        except Exception as e:
            print("FAILED writing file {0} for network {1}, continuing".format(fname, net.code))
            continue

    # our_xml.write("station.xml",format="stationxml", validate=True)
    our_xml.write("station.txt", format="stationtxt")


if __name__ == "__main__":
    main(sys.argv)

