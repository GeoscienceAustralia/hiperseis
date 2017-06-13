"""
Written by Ashby Cooper 2017
"""

import struct
from obspy import UTCDateTime

def decode_gps(block):
    x = struct.unpack('>i', block[0:4])
    day = x[0]
    x = struct.unpack('>i', block[4:8])
    month = x[0]
    x = struct.unpack('>i', block[8:12])
    year = x[0]
    x = struct.unpack('>i', block[12:16])
    hour = x[0]
    x = struct.unpack('>i', block[16:20])
    minute = x[0]
    x = struct.unpack('>i', block[20:24])
    sec = x[0]
    x = struct.unpack('>d', block[24:32])
    lat = x[0]
    x = struct.unpack('>d', block[32:40])
    lng = x[0]
    x = struct.unpack('>d', block[40:48])
    alt = x[0]
    x = struct.unpack('>i', block[48:52])
    clock_error = x[0]

    print(year, month, day, hour, minute, sec)

    try:
        lock_time = UTCDateTime(year, month, day, hour, minute, sec)
    except ValueError:
        # must be something wrong with that entry
        return "EntryError"

    return(lock_time, lat, lng, alt, clock_error)


def decode_logfile(data):

    # dict to store values
    decode_dict = {}

    #read in the data and get all gps syncs
    counter = 0
    iter = True
    while iter == True:
        block = data[89 + (counter * 63):]
        try:
            gps_decode = decode_gps(block)
            if not gps_decode == "EntryError":
                decode_dict[gps_decode[0].timestamp] = {"lat": gps_decode[1],
                                                        "lng": gps_decode[2],
                                                        "elev": gps_decode[3],
                                                        "clock_error":gps_decode[4]}
        except struct.error:
            iter = False

        counter += 1


    return decode_dict