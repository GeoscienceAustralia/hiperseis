import time
from os.path import join, exists, basename, isdir, dirname
from os import remove, mkdir
from struct import error as StructError
import json
import pyasdf
from pyasdf import ASDFWarning
import warnings
from collections import Counter, defaultdict

import glob

from obspy import read_inventory
from obspy.core import inventory, read, UTCDateTime

import sys
import subprocess
from query_input_yes_no import query_yes_no

warnings.filterwarnings("error")

code_start_time = time.time()

# =========================== User Input Required =========================== #

# Path to the data
data_path = '/g/data/ha3/Passive/'

# path to IRIS dataless seed to stationxml converter
seed_xml_conv_path = "/g/data/ha3/axc547/IRIS_SEED_XML/stationxml-converter-1.0.9.jar"

# IRIS Virtual Ntework name
virt_net = '_AusArray'

# FDSN network identifier2
FDSNnetwork = 'OA'

# survey starttime override
# date/time of the overall start for a deployment, this will override the start date/time stored within the dataless SEED metadata
# i.e. useful to exclude old data on SD cards recorded before the deployment if the clear SD command was not sent during instrument deployemnt
# default = None
# UTC Time!!!!
# format = "2017-09-10T00:00:00"
deployment_starttime_override = "2017-09-10T00:00:00"

# =========================================================================== #

# XML_in = join(data_path, virt_net, FDSNnetwork, 'network_metadata', FDSNnetwork + ".xml")
XML_path = join(data_path, virt_net, FDSNnetwork, 'network_metadata')
path_DATA = join(data_path, virt_net, FDSNnetwork, 'raw_DATA/')
ASDF_path_out = join(data_path, virt_net, FDSNnetwork, 'ASDF')

if not exists(ASDF_path_out):
    mkdir(ASDF_path_out)

# JSON filename for network
JSON_out = join(ASDF_path_out, FDSNnetwork + '_raw_dataDB.json')
# ASDF filename for network
ASDF_out = join(ASDF_path_out, FDSNnetwork + '.h5')
# Logfile output
ASDF_log_out = join(ASDF_path_out, FDSNnetwork + '.log')


keys_list = []
info_list = []
station_name_counter = Counter()
station_name_paras = {}

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

# open the station XML into obspy inventory
# inv = read_inventory(XML_in)




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

# function to make a number into a 4 digit string with leading zeros
def make_fourdig(a):
    if len(a) == 1:
        return '000' + a
    elif len(a) == 2:
        return '00' + a
    elif len(a) == 3:
        return '0' + a
    return a


waveforms_added = 0

# Get a list of service directories
service_dir_list = glob.glob(path_DATA + '*')

# counter for number of logfiles
logfile_counter = 0

# dictionary to keep inventory for all stations (default dict)
nscl_inventory_dict = {}

# dictionary to keep end date/start date for each netwrok,station, channel,location
nscl_start_end_dict ={}

#set with stations that are in ASDF file
station_start_end_dict = {}

# iterate through service directories
for service in service_dir_list:

    if not isdir(service):
        continue

    print '\r Processing: ', basename(service)

    station_dir_list = glob.glob(service + '/*')

    # iterate through station directories
    for station_path in station_dir_list:
        station_name = basename(station_path).split("_")[1]

        #get dataless seed file
        xseed_file = glob.glob(join(station_path, '*.dataless'))[0]
        xml_out = join(XML_path, FDSNnetwork[0:2]+station_name+'.xml')

        # shell command to decode the dataless seed
        decode_str = 'java -jar {0} -x -so "Geoscience Australia AusArray" -o {1} {2}'.format(seed_xml_conv_path,
                                                                                              xml_out, xseed_file)

        # decode the dataless seed file
        subprocess.call(decode_str, shell=True)

        # open up the recently created xml file
        station_inv = read_inventory(xml_out)

        # print(station_inv)
        # get station name for metadata
        meta_station_name = station_inv[0][0].code
        meta_network_code = station_inv[0].code

        # get miniseed files
        seed_files = glob.glob(join(station_path, '*miniSEED*/*'))  # '*miniSEED/*.mseed*'))

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

                if len(tr) == 0:
                    continue



                # do some checks to make sure that the network, station, channel, location information is correct
                # Station Name: assign station name in the metadata as correct
                orig_station = tr.stats.station
                # fix the station name
                new_station = meta_station_name
                tr.stats.station = new_station


                # Network Code: assign network name in the metadata as correct
                orig_net = tr.stats.network
                new_net = meta_network_code
                tr.stats.network = new_net



                starttime = tr.stats.starttime.timestamp
                endtime = tr.stats.endtime.timestamp
                orig_chan = tr.stats.channel
                new_chan = orig_chan
                orig_loc = tr.stats.location
                new_loc = orig_loc


                # get the inventory for the station
                sta_sel_inv = station_inv.select(network=new_net, station=new_station)

                # check that the starttime in metadata is within the trace timespan
                # i.e. the data is recorded during or after the clear SD card command was sent
                if not UTCDateTime(endtime) > sta_sel_inv[0][0].start_date:
                    print("trace is outside")
                    continue

                #check if there is a deployment starttime override set
                if not deployment_starttime_override == None:
                    if not UTCDateTime(endtime) > UTCDateTime(deployment_starttime_override):
                        print("trace is outside")
                        continue



                # if starttime < UTCDateTime("2017-08-08T21:06:06.500000Z"):
                #     continue

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
                             "log_filename": ""}



                # get the inventory object for the channel
                select_inv = station_inv.select(network=new_net, station=new_station, channel=new_chan, location=new_loc)

                # print("===================")
                # print(tr)
                # print(select_inv)




                nscl = new_net+"."+new_station+"."+new_chan+"."+new_loc

                # see if network_station_channel_loction is already in start_end dict
                if nscl in nscl_start_end_dict.keys():
                    # compare time to start and end times in dict and see if it is earlier/later
                    stored_starttime = nscl_start_end_dict[nscl][0]
                    stored_endtime = nscl_start_end_dict[nscl][1]
                    if starttime < stored_starttime:
                        nscl_start_end_dict[nscl][0] = starttime
                    elif endtime > stored_endtime:
                        nscl_start_end_dict[nscl][1] = endtime

                else:
                    nscl_start_end_dict[nscl] = [starttime, endtime]

                # see if station is already in start_end dict
                if new_station in station_start_end_dict.keys():
                    # compare time to start and end times in dict and see if it is earlier/later
                    stored_starttime = station_start_end_dict[new_station][0]
                    stored_endtime = station_start_end_dict[new_station][1]
                    if starttime < stored_starttime:
                        station_start_end_dict[new_station][0] = starttime
                    elif endtime > stored_endtime:
                        station_start_end_dict[new_station][1] = endtime

                else:
                    station_start_end_dict[new_station] = [starttime, endtime]


                try:
                    # Add waveform to the ASDF file
                    ds.add_waveforms(tr, tag="raw_recording", labels=[basename(service)])
                except ASDFWarning:
                    # trace already exist in ASDF file!
                    ASDF_log_file.write(filename + '\t' + ASDF_tag + '\t' + "ASDFDuplicateError\n")
                    continue

                waveforms_added += 1
                # try:
                #     #add inventory to ASDF
                #     ds.add_stationxml(select_inv)
                # except:
                #     print("couldnt add inventory")
                #     print select_inv

                # add inventory to dictionary overwrite if there is more than one (i.e. if there are multiple service intervals)
                nscl_inventory_dict[nscl] = select_inv




                keys_list.append(str(ASDF_tag))
                info_list.append(temp_dict)


# list of station level inventories
station_inventories_list =[]

station_inventories_default_dict = defaultdict(list)
station_inv_dict = {}

for nscl, inv in nscl_inventory_dict.iteritems():
    # get the channel level inventory
    print ""
    print nscl

    #modify the start and end dates for the channeel/location level inventory
    inv[0][0][0].start_date = UTCDateTime(nscl_start_end_dict[nscl][0])
    inv[0][0][0].end_date = UTCDateTime(nscl_start_end_dict[nscl][1])
    inv[0][0][0].elevation = inv[0][0].elevation

    print(inv[0][0][0])
    print(inv[0][0][0].start_date)
    print(inv[0][0][0].end_date)

    # add the channels to station_inv_dict
    station_inventories_default_dict[nscl.split(".")[1]].append(inv[0][0][0])

    station_inv_dict[nscl.split(".")[1]] = inv[0][0]

# go through stations
for station, chan_inv_list in station_inventories_default_dict.iteritems():

    print ""
    print station

    #make a new station inventory object with all of the channels and updated start and end dates
    # get the start/end dates from dict
    start_date = UTCDateTime(station_start_end_dict[station][0])
    end_date = UTCDateTime(station_start_end_dict[station][1])

    site = station_inv_dict[station].site

    # make the station_level inventory
    station_inv = inventory.Station(code=station, creation_date=start_date, termination_date=end_date,
                                    start_date=start_date, end_date=end_date,
                                site=site,
                                latitude=station_inv_dict[station].latitude,
                                longitude=station_inv_dict[station].longitude,
                                elevation=station_inv_dict[station].elevation,
                                vault="Transportable Array",
                                channels=station_inventories_default_dict[station],
                                    total_number_of_channels=len(station_inventories_default_dict[station]))

    print(station_inventories_default_dict[station][0])

    print(station_inv)
    station_inventories_list.append(station_inv)


network_start_end = False
# go through station start/end date dict and get the overall start_end date
for key, (start, end) in station_start_end_dict.iteritems():
    if not network_start_end:
        network_start_end = [start, end]
        continue

    if start < network_start_end[0]:
        network_start_end[0] = start
    elif end > network_start_end[1]:
        network_start_end[1] = end


# now make network level inventory
network_inv = inventory.Network(code=FDSNnetwork[0:2], start_date=UTCDateTime(network_start_end[0]),
                                end_date=UTCDateTime(network_start_end[1]),
                                stations=station_inventories_list,
                                total_number_of_stations=len(station_inventories_list))

# create the inventory
inv = inventory.Inventory(networks=[network_inv], source= "Geoscience Australia")


print "+==============================+++"
print ""

print(inv)
print(inv[0])
print(inv[0][0])
print(inv[0][0][0])
print(inv[0][0][1])

XML_file = join(XML_path, FDSNnetwork+'.xml')

if exists(XML_file):
    remove(XML_file)

# write the inventory into the default path
inv.write(path_or_file_object=XML_file, format='STATIONXML')

inv = read_inventory(XML_file)
print "---------"
print(inv)
print(inv[0])
print(inv[0][0])
print(inv[0][0][0])
print(inv[0][0][1])


# add it to ASDF file
ds.add_stationxml(inv)


big_dictionary = dict(zip(keys_list, info_list))

with open(JSON_out, 'w') as fp:
    json.dump(big_dictionary, fp)

del ds
print '\n'

exec_time = time.time() - code_start_time

exec_str = "--- Execution time: %s seconds ---" % exec_time
added_str = '--- Added ' + str(waveforms_added) + ' waveforms to ASDF and JSON database files ---'

print exec_str
print added_str

ASDF_log_file.write(exec_str + '\n')
ASDF_log_file.write(added_str + '\n')

ASDF_log_file.close()
