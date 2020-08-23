"""
Use Json file to store station extra metadata

Fei Zhang
2020-07-08

# TODO: Design a proper json metadata structure to encaptulate all the metadata items. store in the extra tag of stationxml file
"""

import os
import sys
import json


def reformat_orientation_corrections(in_json_file):
    """ Read in a json file produced by AndrewMedlin, which looks like
    {
    "OA.BS24": {
        "date_range": [
            "2017-09-28T11:37:39.060000Z",
            "2018-08-17T22:17:32.680000Z"
        ],
        "azimuth_correction": 5.0
    },
    "OA.BS25": {
        "date_range": [
            "2017-09-28T11:37:43.100000Z",
            "2018-08-17T22:17:36.848000Z"
        ],
        "azimuth_correction": 1.0
    }
    }
    Reformat into
    {
    ORIENT_CORRECTIONS": [
    {
        "network_code":"OA",
        "station_code":"BS24",
        "start_dt": "2017-11-07T09:07:34.930000Z",
        "end_dt":   "2018-08-23T03:52:29.528000Z",
        "azimuth_correction": -5.0
    },
        {
        "network_code":"OA",
        "station_code":"CF28",
        "start_dt": "2017-11-07T09:07:34.930000Z",
        "end_dt":   "2018-08-23T03:52:29.528000Z",
        "azimuth_correction": -5.0
    }
    ]
    }

    :param in_json_file:
    :return:
    """
    # load a file

    with open(in_json_file) as f:
        medata = json.load(f)

    # print(type(medata))

    print(medata.keys())  # dict_keys(['7X.MA01', '7X.MA11', '7X.MA12', '7X.MA13', '7X.MA14', '7X.MA21', '7X.MA22'.....)]


    orient_corr= [ ]
    for netsta in medata.keys():
        # netsta = '7X.MA01'
        network, station = netsta.split('.')

        # print(network, station)

        # print(medata[netsta])

        try:
            az_corr = medata[netsta].get('azimuth_correction')
        except KeyError:
            print("KeyError", az_corr)

        start = medata[netsta]['date_range'][0]
        end = medata[netsta]['date_range'][1]
        # ['2009-09-10T02:57:07.160000Z', '2010-06-05T05:35:00.623400Z'] <class 'list'>
        # print(data[netsta]['azimuth_correction'], type(data[netsta]['azimuth_correction'])) # 38.0 <class 'float'>

        #print(start, end, az_corr)
        # Make a dictionary:
        dicorr = {"network":network, "station":station, "start_dt":start, "end_dt":end, "azimuth_correction":az_corr}
        orient_corr.append(dicorr)

    ORI = {"ORIENT_CORRECTIONS": orient_corr}
    return ORI


def read_my_metajson(metajson_file):
    """
    Read in a new json file proposed to store metadata
    :param metajson_file: a json file
    :return:
    """
    with open(metajson_file) as jsonfile:
        mdata = json.load(jsonfile)

    print(mdata)

    print(mdata.keys())

    print(type(mdata['gps_clock_corrections']))

    for corr in mdata['gps_clock_corrections']:
        print(corr["date"], corr["seconds"])

    return mdata


if __name__ == "__main__":
    json_file = sys.argv[1] # "/Datasets/Orientation_Correction_json/AQT_ori_error_estimates.json" # OA, OA2, 7X_, AQT_
    # json_file = "OA_ori_error_estimates.json" #"7X_ori_error_estimates.json"
    ori=reformat_orientation_corrections(json_file)

    print(json.dumps(ori, indent=2))

    with open("reformatted_orient_corr.json", "w") as f:
        json.dump(ori, f, indent=2)


    # my_meta_json = "extra_station_metadata.json"
    # mdata = read_my_metajson(my_meta_json)

    # json_files = sys.argv[1:]
    # for json_file in json_files:
    #     mdata = read_my_metajson(json_file)
    #
    #     print (mdata)
    #
    # print(json.dumps(mdata, indent=2))
