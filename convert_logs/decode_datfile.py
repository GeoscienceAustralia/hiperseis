#!/usr/bin/python

"""
Mondo script
"""

from __future__ import print_function
import datetime
import math
import os
import struct
import sys
import click
from obspy.core import UTCDateTime


def REC_START(fp):
    block = fp.read(24)
    x = struct.unpack('>i', block[0:4])
    year = x[0]
    x = struct.unpack('>i', block[4:8])
    doy = x[0]
    x = struct.unpack('>i', block[8:12])
    hour = x[0]
    x = struct.unpack('>i', block[12:16])
    minute = x[0]
    x = struct.unpack('>i', block[16:20])
    sec = x[0]
    start = UTCDateTime(1970, 1, 1, 0, 0)


def test_time_fields(year, month, day, hour, minute, sec, seedyear, strtime,
                     bad_time):
    if seedyear > 0:
        val = year - seedyear
        if (val > 1) or (val < -1):
            bad_time = 1
    if year < 2010 or year > 3000:
        bad_time = 1
    if month < 1 or month > 12:
        bad_time = 1
    if day < 1 or month > 31:
        bad_time = 1
    if (hour < 0 or hour > 23):
        bad_time = 1
    if (minute < 0 or minute > 59):
        bad_time = 1
    if (sec < 0 or sec > 59):
        bad_time = 1
    #    if not bad_time ==1:
    str1 = str(year) + "-" + str(month) + "-" + str(day)
    str2 = "%02d:%02d:%02d.000" % (hour, minute, sec)
    strtime = str1 + "T" + str2
    #   else:
    #       strtime=""
    return bad_time, strtime


def decode_gps(fp, bytes, print_bad_temperture, print_bad_gps, mylat, mylng,
               myalt, seedyear, strtime, bad_time, good_gps, bad_gps):

    block = fp.read(60)
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
    x = struct.unpack('>i', block[52:56])
    battery_percent = x[0]
    if battery_percent > 100:
        battery_percent = 0
    x = struct.unpack('>i', block[56:60])
    temperature = x[0]
    bad_time = 0
    if (temperature < -50) or (temperature > 150):
        if print_bad_temperture == 1:
            print("  @@@@ Bad temperture ", temperature, "  at ", fp.tell())
        temperature = 0
    bad_time, strtime = test_time_fields(year, month, day, hour, minute, sec,
                                         seedyear, strtime, bad_time)

    if bad_time == 1:
        if print_bad_gps == 1:
            print("    !!!GPS BAD TIME ", strtime, " at ", fp.tell())
        bad_gps = bad_gps + 1
        return "", mylat, mylng, myalt
    else:
        good_gps = good_gps + 1
        mylat.append(lat)
        mylng.append(lng)
        myalt.append(alt)
        return strtime, mylat, mylng, myalt


def decode_bsn(fp, bytes):
    outcome, block = get_block(fp, bytes)
    if outcome < 1:
        return outcome
    if len(block) < bytes:
        return -2
    x = struct.unpack('>i', block[0:bytes])
    print("Board Serial Number : ", x[0])
    return outcome


def decode_spr(fp, bytes):
    outcome, block = get_block(fp, bytes)
    if outcome < 1:
        return (outcome)
    if len(block) < bytes:
        return (-2)
    x = struct.unpack('>i', block[0:bytes])
    print("Sample Period       : ", x[0])
    return (outcome)


def decode_sms(fp, bytes):
    x1 = fp.tell()
    outcome, block = get_block(fp, bytes)
    if outcome < 1:
        return (outcome)
    if len(block) < bytes:
        return (-2)
    #    block = fp.read(bytes)
    x2 = fp.tell()
    #    x = struct.unpack('>s',block[0:bytes])
    print("Seismometer No.     : ", block)
    return (outcome)


def decode_fwv(fp, bytes):
    outcome, block = get_block(fp, bytes)
    if outcome < 1:
        return (outcome)
    if len(block) < bytes:
        return (-2)
    print("Firmware Version    : ", block)
    return (outcome)


def decode_smm(fp, bytes, typesseis):
    outcome, block = get_block(fp, bytes)
    if outcome < 1:
        return (outcome)
    if len(block) < bytes:
        return (-2)
    #    block = fp.read(bytes)
    x = struct.unpack('>i', block[0:bytes])
    print("Seismometer Type    : ", typesseis[x[0]])
    return (outcome)


def decode_rcs(fp, bytes):
    outcome = get_block(fp, bytes)
    if outcome < 1:
        return (outcome)
    if len(block) < bytes:
        return (-2)
    #    block = fp.read(bytes)
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
    if (year < 2010 or year > 3000):
        bad_time = 1
    elif (month < 0 or month > 12):
        bad_time = 1
    elif (day < 1 or month > 31):
        bad_time = 1
    elif (hour < 0 or hour > 23):
        bad_time = 1
    elif (minute < 0 or minute > 59):
        bad_time = 1
    elif (sec < 0 or sec > 59):
        bad_time = 1
    str1 = str(year) + "-" + str(month) + "-" + str(day)
    str2 = "%02d:%02d:%02d.000" % (hour, minute, sec)
    strtime = str1 + "  " + str2
    if (bad_time == 1):
        print("    !!!RCS BAD TIME ", strtime, " at ", fp.tell())
    else:
        print("RCS Record Start Time   :", strtime)
    return outcome


def decode_rce(fp, bytes):
    outcome, block = get_block(fp, bytes)
    if outcome < 1:
        return outcome
    if len(block) < bytes:
        return -2
    #    block = fp.read(bytes)
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


    if (year < 2010 or year > 3000):
        bad_time = 1
    elif (month < 0 or month > 12):
        bad_time = 1
    elif (day < 1 or month > 31):
        bad_time = 1
    elif (hour < 0 or hour > 23):
        bad_time = 1
    elif (minute < 0 or minute > 59):
        bad_time = 1
    elif (sec < 0 or sec > 59):
        bad_time = 1
    str1 = str(year) + "-" + str(month) + "-" + str(day)
    str2 = "%02d:%02d:%02d.000" % (hour, minute, sec)
    strtime = str1 + "  " + str2
    if (bad_time == 1):
        print("    !!!RCE BAD TIME ", strtime, " at ", fp.tell())
    else:
        print("Record END Time   :", strtime)
    return (outcome)


def decode_udf(fp, bytes):
    outcome, block = get_block(fp, bytes)
    if outcome < 1:
        return outcome
    if len(block) < bytes:
        return -2
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
    if (year < 2010 or year > 3000):
        bad_time = 1
    elif (month < 0 or month > 12):
        bad_time = 1
    elif (day < 1 or month > 31):
        bad_time = 1
    elif (hour < 0 or hour > 23):
        bad_time = 1
    elif (minute < 0 or minute > 59):
        bad_time = 1
    elif (sec < 0 or sec > 59):
        bad_time = 1
    str1 = str(year) + "-" + str(month) + "-" + str(day)
    str2 = "%02d:%02d:%02d.000" % (hour, minute, sec)
    strtime = str1 + "  " + str2
    if bad_time == 1:
        print("    !!!UDF BAD TIME ", strtime, " at ", fp.tell())
    else:
        print("GPS UPDATE FAILED Time   :", strtime)
    return (outcome)


def print_bad_strings(bad_strs, bad_strs_pos):
    print("BAD_STRINGS DUMP")
    print("      string   file-pos")
    for i in range(bad_str_id):
        print(i + 1, "  --->", end=' ')
        for k in range(len(bad_strs[i])):
            sys.stdout.write(bad_strs[i][k])
        print("<---   ", end=' ')
        sys.stdout.write(str(bad_strs_pos[i]))
        print(" ")


def decode_message(outcome, size):
    if outcome == -1:
        print(" *** IOERROR in reading file")
    elif outcome == 0:
        print(" *** End of File reached")
    else:
        print(" Requested ", size, " bytes. Only ", outcome,
              " bytes were available ")


def get_block(fp, bytes):
    try:
        block = fp.read(bytes)
        # ids = struct.unpack('{}s'.format(bytes), block[0:bytes])
        # print(ids, 'in get block')
        if len(block) > 0:
            return len(block), block
        else:
            return 0, None
    except IOError:
        return -1, None


def try_recover_file(fp, x1):
    # rewind the file o last correct code postion + 3 bytes
    fp.seek(3, 1)
    valid_code = 0
    while valid_code < 1:
        x2 = fp.tell()
        outcome, block = get_block(fp, 1)
        if outcome < 1:
            return (outcome)
        if block[0] == 'B':
            outcome = get_block(fp, 2)
            if outcome < 1:
                return (outcome)
            if block[0] == 'S' and block[1] == 'N':
                valid_code = 1
        elif block[0] == 'U':
            outcome = get_block(fp, 2)
            if outcome < 1:
                return (outcome)
            if block[0] == 'D' and block[1] == 'F':
                valid_code = 1
        elif block[0] == 'G':
            outcome = get_block(fp, 2)
            if outcome < 1:
                return (outcome)
            if block[0] == 'P' and block[1] == 'S':
                valid_code = 1
        elif block[0] == 'F':
            outcome = get_block(fp, 2)
            if outcome < 0:
                return (outcome)
            if block[0] == 'W' and block[1] == 'V':
                valid_code = 1
        elif block[0] == 'S':
            outcome = get_block(fp, 1)
            if outcome < 1:
                return (-1)
            if block[0] == 'P':
                outcome = get_block(fp, 1)
                if outcome < 1:
                    return (outcome)
                if block[0] == 'R':
                    valid_code = 1
            elif block[0] == 'M':
                outcome = get_block(fp, 2)
                if outcome < 0:
                    return (outcome)
                if block[0] == 'S' or block[0] == 'M':
                    valid_code = 1
        elif block[0] == 'R':
            outcome = get_block(fp, 1)
            if outcome < 1:
                return (outcome)
            if block[0] == 'C':
                outcome = get_block(fp, 1)
                if outcome < 1:
                    return (outcome)
                if block[0] == 'S' or block[0] == 'E':
                    valid_code = 1
    fp.seek(-3, 1)
    return (10)


def cal_median_value(mylist):
    loclist = sorted(mylist)
    if len(loclist) % 2 == 1:  # odd num of values
        val_lower = loclist[(len(loclist) + 1) / 2 - 1]
        val_upper = val_lower
    else:
        val_lower = loclist[len(loclist) / 2 - 1]
        val_upper = loclist[len(loclist) / 2]
    return val_lower, val_upper


def cal_statistic(mylist):
    if len(mylist) < 1:
        return (-9999, -9999)
    mean = sum(mylist) / len(mylist)
    total = 0
    for i in range(len(mylist)):
        val = mylist[i] - mean
        total = total + (val * val)
    degree = len(mylist) - 1
    valsqr = total / degree
    val = math.sqrt(valsqr)
    return (mean, val)


def set_bit(value, bit):
    value = (value | (1 << bit))
    return value


def test_fileformat_start(fp, filename):
    blk = []
    recstring = ""
    mseedheader = 1
    blk = fp.read(20)
    if len(blk) < 20:
        print("Failed to read first 20 chars")
        mseedheader = -1
        return mseedheader, recstring
    recstring = str(blk[0:5])
    if not recstring.isdigit():
        mseedheader = 0
    if not blk[6] == 'D':
        mseedheader = 0
    if not blk[7] == ' ':
        mseedheader = 0
    # rewind the file
    fp.seek(0, 0)
    return mseedheader, blk


@click.command()
@click.option('-b', '--bad_gps',
              default=False, type=bool, help='Print bad gps info')
@click.option('-u', '--gps_update',
              default=False, type=bool, help='Print the gps update info')
@click.option('-i', '--id_str',
              default=False, type=bool, help='Print bad id strings')
@click.option('-t', '--temperature',
              default=False, type=bool, help='Print bad temperature info')
@click.option('-a', '--all',
              default=False, type=bool, help='Print all')
@click.option('-y', '--year',
              default=False, type=bool, help='Gpsyear.\n'
                                             'max(Gpsyear - year) == 1')
@click.argument('datfile')
def anulog(datfile, bad_gps, id_str, gps_update, temperature, all, year):
    """Program to display contents of the logfile <datfile>.dat"""
    print(datfile)

    current_time = datetime.datetime.now()
    seedyear = current_time.year
    first_rec = 1
    debug = 0
    VALIDCODES = ["BSN", "FWV", "SPR", "SMM", "SMS", "RCS", "RCE", "UDF", "GPS"]
    VALIDSIZES = [4, 22, 4, 4, 10, 24, 24, 24, 60]
    SEISTYPES = ["LennartzLE-3DlITE", "GurupCMG-40T", "GuralpCMG-3ESP",
                 "TrilliumCompact", "MarkL4", "MarkL4C"]
    print_bad_gps = 0
    print_gps_update = 0
    gps_update_failed = 0
    print_bad_temperture = 0
    print_badid_update = 0
    bad_str_id = 0
    strtime = ""
    recoder_restarted_pos = []
    good_gps = 0
    bad_gps = 0
    flagmarker = int(0)
    mylat = []
    mylng = []
    myalt = []
    bad_strs = []
    bad_strs_pos = []
    headertest = 0
    recstr = []
    com_args = sys.argv[1:]
    max_args = len(com_args)
    # j = 0
    # while (j < max_args):
    #     if (com_args[j] == '-badgps' or com_args[j] == '-b'):
    #         if com_args[j + 1] == "yes":
    #             print_bad_gps = 1
    #     elif (com_args[j] == '-update' or com_args[j] == '-u'):
    #         if com_args[j + 1] == "yes":
    #             print_gps_update = 1
    #     elif (com_args[j] == '-idbad' or com_args[j] == '-i'):
    #         if com_args[j + 1] == "yes":
    #             print_badid_update = 1
    #     elif (com_args[j] == '-all' or com_args[j] == '-a'):
    #         if com_args[j + 1] == "yes":
    #             print_badid_update = 1
    #             print_gps_update = 1
    #             print_bad_gps = 1
    #             print_bad_temperture = 1
    #     elif (com_args[j] == '-temp' or com_args[j] == '-t'):
    #         if com_args[j + 1] == "yes":
    #             print_bad_temperture = 1
    #     elif (com_args[j] == '-year' or com_args[j] == '-y'):
    #         seedyear = int(com_args[j + 1])
    #     elif (com_args[j] == '-filename' or com_args[j] == '-f'):
    #         filename = com_args[j + 1]
    #     elif (com_args[j] == '-help' or com_args[j] == '-h'):
    #         pass
    #         sys.exit(3)
    #     j = j + 2

    if not os.path.isfile(datfile):
        print("File ", datfile, "  DOES NOT EXIST")
        sys.exit(4)

    try:
        fp = open(datfile, "rb")
    except IOError:
        print("Failed to open file ", datfile)
        sys.exit(10)
    loop_counter = 0
    x = 0
    # quick test of first lines
    headertest, recstr = test_fileformat_start(fp, datfile)
    if headertest < 0:
        fp.close()
        sys.exit(100)
    elif headertest == 1:
        print("********* WARNNING WARNINIG")
        print(
            "********* WARNNING WARNNING Indication first 8 chars match a minseed file -->",
            recstr)
        print("********* WARNNING WARNNING")

    out_d = {}
    for v in VALIDCODES:
        out_d[v] = 'Not Read or Not available'


    while 1:
        bad_time = 0
        file_postion = fp.tell()
        strtime = ""
        try:
            block = fp.read(3)
            if len(block) < 3:
                break
        except IOError:
            break
        (ids) = struct.unpack('3s', block[0:3])
        id_str = ids[0]

        if str(id_str) == 'BSN':
            recoder_restarted_pos.append(file_postion)
            print(flagmarker, flagmarker & 0x0001)
            if not (flagmarker & 0x0001):
                print("##########################################################")
                print("   Seedyear was ", seedyear,
                      "                    To change use -y year option")
                print("START OF RECORDER\n")
                first_rec = 0
            else:
                print("GPS  TOTAL    ", good_gps + bad_gps)
                print("Good          ", good_gps)
                print("Bad           ", bad_gps)
                print("UPDATE FAILED ", gps_update_failed)
                good_gps = 0
                bad_gps = 0
                gps_update_failed = 0
                fault = 'RECORDER RESTARTED OR EXCHANGED'
                print("**** {} *******************".format(fault))
                out_d[id_str] = fault
            outcome = decode_bsn(fp, 4)

            if outcome < 1:
                decode_message(outcome, 4)
                break
            loop_counter = 0
            flagmarker = set_bit(flagmarker, 0)
            out_d[id_str] = outcome

        elif id_str == 'SPR':
            flagmarker = set_bit(flagmarker, 2)
            outcome = decode_spr(fp, 4)
            if outcome < 1:
                decode_message(outcome, 4)
                break
            loop_counter = 0
            out_d[id_str] = outcome

        elif id_str == 'SMM':
            flagmarker = set_bit(flagmarker, 3)
            outcome = decode_smm(fp, 4, SEISTYPES)
            if outcome < 1:
                decode_message(outcome, 4)
                break
            loop_counter = 0
            out_d[id_str] = outcome
        elif id_str == 'FWV':
            flagmarker = set_bit(flagmarker, 1)
            outcome = decode_fwv(fp, 22)
            if outcome < 1:
                decode_message(outcome, 22)
                break
            loop_counter = 0
            out_d[id_str] = outcome
        elif id_str == 'SMS':
            flagmarker = set_bit(flagmarker, 4)
            outcome = decode_sms(fp, 4)
            if outcome < 1:
                decode_message(outcome, 4)
                break
            loop_counter = 0
            out_d[id_str] = outcome
        elif id_str == 'RCS':
            flagmarker = set_bit(flagmarker, 5)
            outcome = decode_rcs(fp, 24)
            if outcome < 1:
                decode_message(outcome, 24)
                break
            loop_counter = 0
            out_d[id_str] = outcome
        elif id_str == 'RCE':
            flagmarker = set_bit(flagmarker, 6)
            outcome = decode_rce(fp, 24)
            if outcome < 1:
                decode_message(outcome, 24)
                break
            loop_counter = 0
            out_d[id_str] = outcome
        elif id_str == 'UDF':
            flagmarker = set_bit(flagmarker, 7)
            outcome = decode_udf(fp, 24)
            if outcome < 1:
                decode_message(outcome, 24)
                break
            gps_update_failed = gps_update_failed + 1
            loop_counter = 0
            out_d[id_str] = outcome
        elif id_str == 'GPS':
            flagmarker = set_bit(flagmarker, 8)
            first_gps = 0
            (strtime, mylat, mylng, myalt) = decode_gps(fp, 60,
                                                        print_bad_temperture,
                                                        print_bad_gps, mylat, mylng,
                                                        myalt, seedyear, strtime,
                                                        bad_time)
            if len(strtime) > 0:
                if (not strtime[1] == 0) and (not strtime[2] == 0) and \
                        (not strtime[2] == 0):
                    if print_gps_update == 1:
                        print("GPS UPDATED ", strtime)
                    dt = UTCDateTime(strtime) - UTCDateTime(1970, 1, 1, 0, 0)
                    if not x:
                        x = dt
            loop_counter = 0
            out_d[id_str] = 'Not defined'
        else:
            bad_str_id = bad_str_id + 1
            bad_strs.append(id_str)
            bad_strs_pos.append(file_postion)
            if loop_counter == 3:
                fp.close()
                sys.exit(1)
            #        if print_badid_update == 1:
            #            print "Bad id str ",id_str," at pos ",fp.tell(),"  file marker ",file_postion
            debug = 1
            outcome = try_recover_file(fp, file_postion)
            if outcome < 0:
                print("Error reading file")
            elif outcome == 0:
                print("End of File reached")
                break
            #        if print_badid_update == 1:
            #            print "       Recoverd at pos ",fp.tell()
            loop_counter = loop_counter + 1
    print("GPS:  TOTAL      Good        Bad      Update failed")
    print("     ", good_gps + bad_gps, "     ", good_gps, "      ", bad_gps,
          "       ", gps_update_failed)
    print("BAD_STR_ID\'S  ", bad_str_id)
    mean, val = cal_statistic(mylat)
    strlat = str('%.5f' % mean) + " +/- " + str('%.5f' % val)
    mean, val = cal_statistic(mylng)
    strlng = str('%.5f' % mean) + " +/- " + str('%.5f' % val)
    mean, val = cal_statistic(myalt)
    stralt = str('%d' % mean) + " +/- " + str('%d' % val)
    if good_gps > 0:
        #    print "Geo Coords:   ",'%.5f'%(sum(mylat)/len(mylat)),"  ",'%.5f'%(sum(mylng)/len(mylng)),"   ",int(sum(myalt)/len(myalt))
        print("GEO COORDS")
        print("       MEANS VALS:   ", strlat, "  ", strlng, "   ", stralt)
        val1, val2 = cal_median_value(mylat)
        median_lat = (val1 + val2) / 2.0
        val1, val2 = cal_median_value(mylng)
        median_lng = (val1 + val2) / 2.0
        val1, val2 = cal_median_value(myalt)
        median_alt = (val1 + val2) / 2.0
        print("      MEDIAN VALS:", '   %.5f' % median_lat, "              ",
              '%.5f' % median_lng, "               ", '%.2f' % median_alt)
    print("######################################################################")
    if print_badid_update == 1:
        print_bad_strings(bad_strs, bad_strs_pos)
    restarts = len(recoder_restarted_pos)
    if restarts > 0:
        if not recoder_restarted_pos[0] == 0:
            print("Hit of corruption at start of file. The BSN entry")
            print("is not at the start of the file ")
        if restarts > 1:
            print("WARNING File ", datfile, " restarted!. Marks at location ")
            for i in range(1, restarts):
                print(i, ":     ", recoder_restarted_pos[i], " byte")
    # the 2 allowed values for flagmarker!!!
    x1 = 0xb111111101
    x2 = 0xb111111111

    if not flagmarker == x1 and not flagmarker == x2:
        print("LOGFILE HAS MISSING FLAGS:")
    for i in range(0, 9):
        x = (1 << i)
        if not (flagmarker & x):
            if i == 0:
                print("  Board Serial Number          NOT RECORDED")
            elif i == 1:
                print("  Firmware Version Number      NOT RECORDED")
            elif i == 2:
                print("  Sample Rate                  NOT RECORDED")
            elif i == 3:
                print("  Seismometer Mode Code        NOT RECORDED")
            elif i == 4:
                print("  Seismometer Serial Number    NOT RECORDED")
            elif i == 5:
                print("  Record Start Time            NOT RECORDED")
            elif i == 6:
                print("  Record End Time              NOT RECORDED")
            elif i == 8:
                print("  Gps Updates                  NOT RECORDED")
    fp.close()

    print(out_d)


if __name__ == '__main__':
    anulog()