#!/usr/bin/python
"""
Use this script to convert the binary .dat ANU log files into jsons.
"""

from __future__ import print_function
import os
from os.path import split, splitext, join, abspath, basename
import math
import struct
import sys
import json
import logging
import click
import datetime
import glob
from joblib import Parallel, delayed


logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


# constants
YEAR = datetime.datetime.now().year
VALIDCODES = ["BSN", "FWV", "SPR", "SMM", "SMS", "RCS", "RCE", "UDF", "GPS"]
VALIDSIZES = [4, 22, 4, 4, 10, 24, 24, 24, 60]
SEISTYPES = ["LennartzLE-3DlITE", "GurupCMG-40T", "GuralpCMG-3ESP",
             "TrilliumCompact", "MarkL4", "MarkL4C"]


def test_time_fields(year, month, day, hour, minute, sec):
    try:
        gps_time = datetime.datetime(year=year, month=month, day=day,
                                     hour=hour, minute=minute, second=sec)
        return gps_time.strftime('%d/%m/%Y %H:%M:%S')
    except ValueError:
        return False


def decode_gps(fp, bytes, print_bad_temperture, print_bad_gps, mylat, mylng,
               myalt, clock, battery, temp, seedyear, good_gps, bad_gps):

    block = fp.read(bytes)
    day = struct.unpack('>i', block[0:4])[0]
    month = struct.unpack('>i', block[4:8])[0]
    year = struct.unpack('>i', block[8:12])[0]
    good_year = True if abs(year-seedyear) <= 1 else False

    hour = struct.unpack('>i', block[12:16])[0]
    minute = struct.unpack('>i', block[16:20])[0]
    sec = struct.unpack('>i', block[20:24])[0]
    lat = struct.unpack('>d', block[24:32])[0]
    lng = struct.unpack('>d', block[32:40])[0]
    alt = struct.unpack('>d', block[40:48])[0]
    clock_error = struct.unpack('>i', block[48:52])[0]
    battery_percent = struct.unpack('>i', block[52:56])[0]
    if battery_percent > 100:
        battery_percent = 0
    temperature = struct.unpack('>i', block[56:60])[0]

    if (temperature < -50) or (temperature > 150):
        if print_bad_temperture:
            log.warning("  @@@@ Bad temperture {} at {}".format(
                temperature, fp.tell()))
        temperature = 0

    strtime = test_time_fields(year, month, day, hour, minute, sec)

    if good_year and strtime:
        good_gps += 1
        mylat.append(lat)
        mylng.append(lng)
        myalt.append(alt)
        clock.append(clock_error)
        battery.append(battery_percent)
        temp.append(temperature)
    else:
        bad_gps += 1
        bad_time_str = "{}/{}/{} {}:{}:{}".format(day, month, year, hour,
                                                  minute, sec)
        if print_bad_gps:
            log.warning("    !!!GPS BAD TIME {strtime} at {location}".format(
                strtime=bad_time_str, location=fp.tell()))

    return strtime, mylat, mylng, myalt, clock, battery, temp, \
        good_gps, bad_gps


def decode_bsn(fp, bytes):
    outcome, block = get_block(fp, bytes)
    if outcome >= 1:
        x = struct.unpack('>i', block[0:bytes])
        log.info("BSN: Board Serial Number: {}".format(x[0]))
        return outcome, x[0]
    else:
        return outcome, None


def decode_spr(fp, bytes):
    outcome, block = get_block(fp, bytes)
    if outcome >= 1:
        x = struct.unpack('>i', block[0:bytes])
        log.info("SPR: Sample Period: {}".format(x[0]))
        return outcome, x[0]
    else:
        return outcome, None


def decode_sms(fp, bytes):
    x1 = fp.tell()
    outcome, block = get_block(fp, bytes)
    if outcome >= 1:
        x2 = fp.tell()
        x = struct.unpack('>i', block[0:bytes])
        log.info("SMS: Seismometer No: {}".format(block))
        return outcome, x[0]
    else:
        return outcome, None


def decode_fwv(fp, bytes):
    outcome, block = get_block(fp, bytes)
    if outcome >= 1:
        log.info("FWV: Firmware Version    :{}".format(block))
        return outcome, block
    else:
        return outcome, None


def decode_smm(fp, bytes):
    outcome, block = get_block(fp, bytes)
    if outcome >= 1:
        x = struct.unpack('>i', block[0:bytes])
        log.info("SMM: Seismometer Type: {}".format(SEISTYPES[x[0]]))
        return outcome, SEISTYPES[x[0]]
    else:
        return outcome, None


def decode_rcs(fp, bytes):
    outcome, block = get_block(fp, bytes)
    str_start_time = None
    if outcome >= 1:
        try:
            str_start_time = _unpack_time(block)
            log.info("RCS: RCS Record Start Time: {}".format(str_start_time))
        except:
            str_start_time = 'Bad Recode START Time/RCS time'
            log.warning(str_start_time)

    return outcome, str_start_time


def decode_rce(fp, bytes):
    outcome, block = get_block(fp, bytes)
    str_stop_time = None
    if outcome >= 1:
        try:
            str_stop_time = _unpack_time(block)
            log.info("RCE: Record END Time: {}".format(str_stop_time))
        except ValueError:
            str_stop_time = 'Bad Recode END Time/RCE time'
            log.warning(str_stop_time)

    return outcome, str_stop_time


def _unpack_time(block):
    x = struct.unpack('>i', block[0:4])
    year = x[0]
    x = struct.unpack('>i', block[4:8])
    days = x[0]
    x = struct.unpack('>i', block[8:12])
    hours = x[0]
    x = struct.unpack('>i', block[12:16])
    minutes = x[0]
    x = struct.unpack('>i', block[16:20])
    secs = x[0]
    x = struct.unpack('>i', block[20:24])
    usec = x[0]
    start_time = datetime.datetime(year, 1, 1)
    start_time += datetime.timedelta(days=days - 1,
                                     hours=hours,
                                     minutes=minutes,
                                     seconds=secs,
                                     microseconds=usec)
    return start_time.strftime('%d/%m/%Y %H:%M:%S:%f')


def decode_udf(fp, bytes):
    outcome, block = get_block(fp, bytes)
    str_udf_time = None
    if outcome >= 1:
        try:
            str_udf_time = _unpack_time(block)
            log.info("UDF: GPS UPDATE FAILED Time: {}".format(str_udf_time))
        except ValueError:
            str_udf_time = "   !!!UDF BAD TIME at {}".format(fp.tell())
            log.warning(str_udf_time)

    return outcome, str_udf_time


def print_bad_strings(bad_strs, bad_strs_pos, bad_str_id):
    log.info("BAD_STRINGS DUMP")
    log.info("      string   file-pos")
    for i in range(bad_str_id):
        print(i + 1, "  --->", end=' ')
        for k in range(len(bad_strs[i])):
            sys.stdout.write(bad_strs[i][k])
        print("<---   ", end=' ')
        sys.stdout.write(str(bad_strs_pos[i]))
        print(" ")


def decode_message(outcome, size):
    if outcome == -1:
        log.warning(" *** IOERROR in reading file")
    elif outcome == 0:
        log.warning(" *** End of File reached")
    else:
        log.warning(" Requested {size} bytes. Only {outcome} bytes were "
                    "available".format(size=size, outcome=outcome))


def get_block(fp, bytes):
    try:
        block = fp.read(bytes)
        # ids = struct.unpack('{}s'.format(bytes), block[0:bytes])
        # print(ids, 'in get block')
        if len(block) > 0:
            return len(block), block
        elif len(block) < bytes:
            return -2, None
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
            return outcome
        if block[0] == 'B':
            outcome, block = get_block(fp, 2)
            if outcome < 1:
                return outcome
            if block[0] == 'S' and block[1] == 'N':
                valid_code = 1
        elif block[0] == 'U':
            outcome, block = get_block(fp, 2)
            if outcome < 1:
                return outcome
            if block[0] == 'D' and block[1] == 'F':
                valid_code = 1
        elif block[0] == 'G':
            outcome, block = get_block(fp, 2)
            if outcome < 1:
                return outcome
            if block[0] == 'P' and block[1] == 'S':
                valid_code = 1
        elif block[0] == 'F':
            outcome, block = get_block(fp, 2)
            if outcome < 0:
                return outcome
            if block[0] == 'W' and block[1] == 'V':
                valid_code = 1
        elif block[0] == 'S':
            outcome, block = get_block(fp, 1)
            if outcome < 1:
                return -1
            if block[0] == 'P':
                outcome, block = get_block(fp, 1)
                if outcome < 1:
                    return outcome
                if block[0] == 'R':
                    valid_code = 1
            elif block[0] == 'M':
                outcome, block = get_block(fp, 2)
                if outcome < 0:
                    return outcome
                if block[0] == 'S' or block[0] == 'M':
                    valid_code = 1
        elif block[0] == 'R':
            outcome, block = get_block(fp, 1)
            if outcome < 1:
                return outcome
            if block[0] == 'C':
                outcome, block = get_block(fp, 1)
                if outcome < 1:
                    return outcome
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


def test_fileformat_start(fp):
    blk = []
    recstring = ""
    mseedheader = 1
    blk = fp.read(20)
    if len(blk) < 20:
        log.warning("Failed to read first 20 chars")
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
@click.option('-a', '--all_print',
              default=False, type=bool, help='Print all')
@click.option('-y', '--year',
              default=YEAR,
              type=click.IntRange(2000, YEAR+1),
              help='Gpsyear. max(Gpsyear - year) == 1')
@click.option('-o', '--output_dir',
              default=None,
              type=click.Path(writable=True, file_okay=False),
              help='Output dir name. If no output dir is provided, \n'
                   'input dir will be used.')
@click.argument('datfile', type=click.Path(exists=True))
def anulog(datfile, bad_gps, id_str, gps_update,
           temperature, all_print, year, output_dir):
    """Program to display contents of the logfile <datfile>.dat"""

    # if a dir if provided, get all .dat files from that dir
    if os.path.isdir(datfile):
        datfiles = glob.glob(os.path.join(datfile, '*.dat'))
    else:  # otherwise if file if provided
        assert os.path.isfile(datfile)
        datfiles = [datfile]

    # n_jobs=-1 uses all cpus
    dicts = Parallel(n_jobs=-1)(delayed(decode_anulog)(
        d, bad_gps, id_str, gps_update, temperature, all_print, year)
                               for d in datfiles)

    # use input dir as default output dir
    basedir = os.path.abspath(output_dir) if output_dir \
        else os.path.split(datfiles[0])[0]

    if not os.path.exists(basedir):
        log.info('Supplied output dir does not exist. It will be created.')
        os.makedirs(basedir)

    for d, dd in zip(datfiles, dicts):
        output_file = join(basedir, splitext(basename(d))[0]) + '.json'
        json.dump(dd, open(output_file, 'w'))


def decode_anulog(datfile, bad_gps, id_str, gps_update,
                  temperature, all_print, year):
    """
    Parameters
    ----------
    datfile: str or file pointer
        path to the binary datafile that needs to be converted
    bad_gps: bool, optional, default False
        whether to report problems in data conversion during processing
    id_str: bool, optional, default False
        whether to report problematic id strings during processing
    gps_update:
    temperature: bool, optional
        whether to report temperature conversion problems during processing
    all_print: bool, optional
        whether to print all log messages
    year: int, optional
        optional GPS year
    Returns
    -------
    out_d: dict
        dictionary containing all the log data

    """

    # if a file path is sent instead of a file pointer
    if os.path.isfile(datfile):
        datfile = open(datfile, 'r')

    gps_update_failed = 0
    bad_str_id = 0
    recoder_restarted_pos = []
    good_gps_count = 0
    bad_gps_count = 0
    flagmarker = int(0)
    mylat = []
    mylng = []
    myalt = []
    clock = []
    battery = []
    temp = []
    bad_strs = []
    bad_strs_pos = []
    print_bad_gps = bad_gps
    print_gps_update = gps_update
    print_badid_update = id_str
    print_bad_temperture = temperature
    if all_print:
        print_badid_update = True
        print_gps_update = True
        print_bad_gps = True
        print_bad_temperture = True
    loop_counter = 0
    # quick test of first lines
    headertest, recstr = test_fileformat_start(datfile)
    if headertest < 0:
        datfile.close()
        sys.exit(100)
    elif headertest == 1:
        print("********* WARNNING WARNINIG")
        print("********* WARNNING WARNNING Indication first 8 chars match "
              "a minseed file -->", recstr)
        print("********* WARNNING WARNNING")
    out_d = {}
    for v in VALIDCODES:
        out_d[v] = 'Not Read or Not available'
    while True:
        file_postion = datfile.tell()
        strtime = ""
        try:
            block = datfile.read(3)
            if len(block) < 3:
                break
        except IOError:
            break
        (ids) = struct.unpack('3s', block[0:3])
        id_str = ids[0]

        if str(id_str) == 'BSN':
            recoder_restarted_pos.append(file_postion)
            if not (flagmarker & 0x0001):
                print("######################################################")
                print("   Seedyear was ", year,
                      "                    To change use -y year option")
                print("START OF RECORDER\n")
            else:
                print("GPS  TOTAL    ", good_gps_count + bad_gps_count)
                print("Good          ", good_gps_count)
                print("Bad           ", bad_gps_count)
                print("UPDATE FAILED ", gps_update_failed)
                good_gps_count = 0
                bad_gps_count = 0
                gps_update_failed = 0
                fault = 'RECORDER RESTARTED OR EXCHANGED'
                print("**** {} *******************".format(fault))
            outcome, out_d[id_str] = decode_bsn(datfile, 4)

            if outcome < 1:
                decode_message(outcome, 4)
                break
            loop_counter = 0
            flagmarker = set_bit(flagmarker, 0)

        elif id_str == 'SPR':
            flagmarker = set_bit(flagmarker, 2)
            outcome, out_d[id_str] = decode_spr(datfile, 4)
            if outcome < 1:
                decode_message(outcome, 4)
                break
            loop_counter = 0

        elif id_str == 'SMM':
            flagmarker = set_bit(flagmarker, 3)
            outcome, out_d[id_str] = decode_smm(datfile, 4)
            if outcome < 1:
                decode_message(outcome, 4)
                break
            loop_counter = 0

        elif id_str == 'FWV':
            flagmarker = set_bit(flagmarker, 1)
            outcome, out_d[id_str] = decode_fwv(datfile, 22)
            if outcome < 1:
                decode_message(outcome, 22)
                break
            loop_counter = 0

        elif id_str == 'SMS':
            flagmarker = set_bit(flagmarker, 4)
            outcome, out_d[id_str] = decode_sms(datfile, 4)
            if outcome < 1:
                decode_message(outcome, 4)
                break
            loop_counter = 0

        elif id_str == 'RCS':
            flagmarker = set_bit(flagmarker, 5)
            outcome, out_d[id_str] = decode_rcs(datfile, 24)

            if outcome < 1:
                decode_message(outcome, 24)
                break
            loop_counter = 0

        elif id_str == 'RCE':
            flagmarker = set_bit(flagmarker, 6)
            outcome, out_d[id_str] = decode_rce(datfile, 24)
            if outcome < 1:
                decode_message(outcome, 24)
                break
            loop_counter = 0

        elif id_str == 'UDF':
            flagmarker = set_bit(flagmarker, 7)
            outcome, out_d[id_str] = decode_udf(datfile, 24)
            if outcome < 1:
                decode_message(outcome, 24)
                break
            gps_update_failed = gps_update_failed + 1
            loop_counter = 0

        elif id_str == 'GPS':
            flagmarker = set_bit(flagmarker, 8)
            strtime, mylat, mylng, myalt, clock, battery, temp, \
            good_gps_count, bad_gps_count = \
                decode_gps(datfile, 60, print_bad_temperture, print_bad_gps,
                           mylat, mylng, myalt, clock, battery, temp, year,
                           good_gps_count, bad_gps_count)
            if strtime:
                if (not strtime[1] == 0) and (not strtime[2] == 0) and \
                        (not strtime[2] == 0):
                    if print_gps_update == 1:
                        log.info("GPS UPDATED {}".format(strtime))
            loop_counter = 0
        else:
            bad_str_id = bad_str_id + 1
            bad_strs.append(id_str)
            bad_strs_pos.append(file_postion)
            if loop_counter == 3:
                datfile.close()
                sys.exit(1)

            outcome = try_recover_file(datfile, file_postion)
            if outcome < 0:
                print("Error reading file")
            elif outcome == 0:
                print("End of File reached")
                break
            loop_counter = loop_counter + 1
    print("GPS:  TOTAL      Good        Bad      Update failed")
    print("     ", good_gps_count + bad_gps_count, "     ",
          good_gps_count, "      ", bad_gps_count,
          "       ", gps_update_failed)
    print("BAD_STR_ID\'S  ", bad_str_id)
    mean, val = cal_statistic(mylat)
    strlat = str('%.5f' % mean) + " +/- " + str('%.5f' % val)
    mean, val = cal_statistic(mylng)
    strlng = str('%.5f' % mean) + " +/- " + str('%.5f' % val)
    mean, val = cal_statistic(myalt)
    stralt = str('%d' % mean) + " +/- " + str('%d' % val)
    if good_gps_count > 0:
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
    print("##################################################################")
    if print_badid_update == 1:
        print_bad_strings(bad_strs, bad_strs_pos, bad_str_id)
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
    # x1 = 0xb111111101
    # x2 = 0xb111111111
    #
    # if not flagmarker == x1 and not flagmarker == x2:
    #     print("LOGFILE HAS MISSING FLAGS:")
    # TODO: fix flagmarker check
    """
        Those are really large hex numbers (760495542529 and 760495542545). 
        I have not seens the flagmaker equal 383, 447 (later being when 'UDF' 
        flag is present) 
        """
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
                print("  Seismometer Type             NOT RECORDED")
            elif i == 4:
                print("  Seismometer Serial Number    NOT RECORDED")
            elif i == 5:
                print("  Record Start Time            NOT RECORDED")
            elif i == 6:
                print("  Record End Time              NOT RECORDED")
            elif i == 8:
                print("  Gps Updates                  NOT RECORDED")
    datfile.close()
    out_d['GPS'] = {'ALTITUDE': myalt, 'LATITUDE': mylat, 'LONGITUDE': mylng,
                    'CLOCK': clock, 'BATTERY': battery, 'TEMPERATURE': temp}
    return out_d


if __name__ == '__main__':
    anulog()
