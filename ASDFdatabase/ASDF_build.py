import time
from os.path import join, exists, basename, isdir, dirname
from os import remove, mkdir
from struct import error as StructError
import json
import pyasdf
from pyasdf import ASDFWarning
import warnings

import glob

import numpy as np

from obspy.core import read

import sys
from query_input_yes_no import query_yes_no

warnings.filterwarnings("error")

code_start_time = time.time()

# =========================== User Input Required =========================== #

# Path to the data
data_path = '/g/data/ha3/Passive/'

# IRIS Virtual Ntework name
virt_net = '_GA_ANUtest'

# FDSN network identifier
FDSNnetwork = 'XX'

# =========================================================================== #

path_XML = join(data_path, virt_net, FDSNnetwork, 'network_metadata/stnXML', FDSNnetwork + '.xml')
path_DATA = join(data_path, virt_net, FDSNnetwork, 'raw_DATA/')
ASDF_path_out = join(data_path, virt_net, FDSNnetwork, 'ASDF')

if not exists(ASDF_path_out):
    mkdir(ASDF_path_out)

#JSON filename for network
JSON_out = join(ASDF_path_out, FDSNnetwork + '_raw_dataDB.json')
#ASDF filename for network
ASDF_out = join(ASDF_path_out, FDSNnetwork + '.h5')
#Logfile output
ASDF_log_out = join(ASDF_path_out, FDSNnetwork + '.log')

keys_list = []
info_list =[]

# remove log file if it exists
if exists(ASDF_log_out):
    remove(ASDF_log_out)

# query the user to overwrite JSON database file or not
if exists(JSON_out):
    delete_queary = query_yes_no("Remove Existing JSON database file?")
    if delete_queary == 'yes':
        # removing existing SQLdb
        remove(JSON_out)
    elif delete_queary == 'no':
        sys.exit(0)

# query the user to overwrite the ASDF database or not
if exists(ASDF_out):
    delete_queary = query_yes_no("Remove Existing ASDF File?")
    if delete_queary == 'yes':
        # removing existing ASDF
        remove(ASDF_out)
    elif delete_queary == 'no':
        sys.exit(0)


# create the log file
ASDF_log_file = open(ASDF_log_out, 'w')

# Create/open the ASDF file
ds = pyasdf.ASDFDataSet(ASDF_out, compression="gzip-3")

# Add the station XML data to the ASDF file
ds.add_stationxml(path_XML)


# function to create the ASDF waveform ID tag
def make_ASDF_tag(tr, tag):
# def make_ASDF_tag(ri, tag):
    data_name = "{net}.{sta}.{loc}.{cha}__{start}__{end}__{tag}".format(
        net=tr.stats.network,
        sta=tr.stats.station,
        loc=tr.stats.location,
        cha=tr.stats.channel,
        start=tr.stats.starttime.strftime("%Y-%m-%dT%H:%M:%S"),
        end=tr.stats.endtime.strftime("%Y-%m-%dT%H:%M:%S"),
        tag=tag)
    return data_name


waveforms_added = 0

# Get a list of service directories
service_dir_list = glob.glob(path_DATA + '*')


# iterate through service directories
for service in service_dir_list:

    if not isdir(service):
        continue

    print '\r Processing: ', basename(service)

    station_dir_list = glob.glob(service + '/*')

    # iterate through station directories
    for station_path in station_dir_list:
        station_name = basename(station_path)

        # get the logfile
        anu_logfile = glob.glob(join(station_path, '*.dat'))

        if len(anu_logfile) == 0:
            continue

        # get miniseed files
        seed_files = glob.glob(join(station_path, '*miniSEED/*'))  # '*miniSEED/*.mseed*'))

        if len(seed_files) == 0:
            continue

        print '\r Working on station: ', station_name

        # Iterate through the miniseed files, fix the header values and add waveforms
        for _i, filename in enumerate(seed_files):
            print "\r     Parsing miniseed file ", _i + 1, ' of ', len(seed_files), ' ....',
            sys.stdout.flush()

            try:
                # Read the stream
                st = read(filename)

            except (TypeError, StructError) as e:
                # the file is not miniseed
                ASDF_log_file.write(filename + '\t' + "TypeError\n")

            # iterate through traces in st (there will usually be only one trace per stream,
            # however if there are problems with the miniseed files - like much of the ANU data - there
            # can be more than one trace in a miniseed file (seperated by a large time gap)
            for tr in st:
                waveforms_added += 1

                # do some checks to make sure that the network, station, channel, location information is correct
                # Station Name: for now just assume that what is in the traces is correct
                orig_station = tr.stats.station
                new_station = orig_station
                # Network Code: the network code in the miniseed header is prone to user error
                # (i.e. whatever the operator entered into the instrument in the field)
                orig_net = tr.stats.network
                # use the first two characters as network code. Temporary networks have start and end year as well
                new_net = FDSNnetwork[:2]
                # overwrite network code in miniseed header
                tr.stats.network = new_net

                starttime = tr.stats.starttime.timestamp
                endtime = tr.stats.endtime.timestamp
                orig_chan = tr.stats.channel
                new_chan = orig_chan
                orig_loc = tr.stats.location
                new_loc = orig_loc

                # The ASDF formatted waveform name [full_id, station_id, starttime, endtime, tag]
                ASDF_tag = make_ASDF_tag(tr, "raw_recording").encode('ascii')

                # make a dictionary for the trace that will then be appended to a larger dictionary for whole network
                temp_dict = {"tr_starttime": starttime,
                             "tr_endtime": endtime,
                             "orig_network": str(orig_net),
                             "new_network": str(new_net),
                             "orig_station": str(orig_station),
                             "new_station": str(new_station),
                             "orig_channel": str(orig_chan),
                             "new_channel": str(new_chan),
                             "orig_location": str(orig_loc),
                             "new_location": str(new_loc),
                             "seed_path": str(dirname(filename)),
                             "seed_filename": str(basename(filename)),
                             "log_filename": str(anu_logfile[0])}

                try:
                    # Add waveform to the ASDF file
                    ds.add_waveforms(tr, tag="raw_recording")
                except ASDFWarning:
                    #trace already exist in ASDF file!
                    ASDF_log_file.write(filename + '\t' + ASDF_tag + '\t' + "ASDFDuplicateError\n")

            keys_list.append(str(ASDF_tag))
            info_list.append(temp_dict)




big_dictionary = dict(zip(keys_list, info_list))



with open(JSON_out, 'w') as fp:
    json.dump(big_dictionary, fp)


del ds
print '\n'

exec_time  = time.time() - code_start_time

exec_str = "--- Execution time: %s seconds ---" % exec_time
added_str = '--- Added ' + str(waveforms_added) + ' waveforms to ASDF and JSON database files ---'

print exec_str
print added_str

ASDF_log_file.write(exec_str + '\n')
ASDF_log_file.write(added_str + '\n')

ASDF_log_file.close()