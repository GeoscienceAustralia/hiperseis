import time
from os.path import join, exists, basename, isdir, dirname
from os import remove, mkdir
from struct import error as StructError
import json
import pyasdf
from pyasdf import ASDFWarning
import warnings
from collections import Counter, defaultdict
from convert_logs.decode_datfile import decode_anulog

import glob

import numpy as np

from obspy.core import inventory, read, UTCDateTime

import sys
from query_input_yes_no import query_yes_no

warnings.filterwarnings("error")

code_start_time = time.time()

# =========================== User Input Required =========================== #

# Path to the data
data_path = '/g/data/ha3/Passive/'

# IRIS Virtual Ntework name
virt_net = '_ANU'

# FDSN network identifier2
FDSNnetwork = '7G(2013-2015)'

# rough year of survey
rough_year = 2013

# tolerence (degrees) for latitude and longitude when working out if two stations should have the same name
# i.e. 100m = 0.001 degrees - if two stations are seperated by less than this then they are the same station
tol = 0.001

# =========================================================================== #

grid_file = join(data_path, 'AUS_Seismic_MT_grid/AUS_seismic_MT_grid.txt')
XML_path_out = join(data_path, virt_net, FDSNnetwork, 'network_metadata')
path_DATA = join(data_path, virt_net, FDSNnetwork, 'raw_DATA/')
ASDF_path_out = join(data_path, virt_net, FDSNnetwork, 'ASDF')

if not exists(XML_path_out):
    mkdir(XML_path_out)

if not exists(ASDF_path_out):
    mkdir(ASDF_path_out)

# JSON filename for network
JSON_out = join(ASDF_path_out, FDSNnetwork + '_raw_dataDB.json')
# ASDF filename for network
ASDF_out = join(ASDF_path_out, FDSNnetwork + '.h5')
# Logfile output
ASDF_log_out = join(ASDF_path_out, FDSNnetwork + '.log')

# open up the AUS_Seismic_MT grid text file
with open(grid_file, 'r') as f:
    data = f.readlines()[1:]

AUS_grid_names = []
AUS_grid_latitude = []
AUS_grid_longitude = []

for line in data:
    fields = line.rstrip('\r\n').split(',')

    AUS_grid_names.append(fields[1])
    AUS_grid_longitude.append(float(fields[2]))
    AUS_grid_latitude.append(float(fields[3]))

lat_array = np.array(AUS_grid_latitude)
lon_array = np.array(AUS_grid_longitude)

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
station_inventory_dict = defaultdict(list)

# dictionary to keep end date/start date for each station
station_start_end_dict ={}

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

        # if the logfile doesnt exist then we have no idea on coordinates or station name
        if len(anu_logfile) == 0:
            ASDF_log_file.write(station_path + '\t' + 'NoLogfile\n')
            continue

        # logfile Unique ID so that they can be differentiated in ASDF auxillary data
        logfile_id = "UID" + make_fourdig(str(logfile_counter))
        # print(logfile_id)

        logfile_counter += 1

        # decode the logfile into dictionary
        try:
            logfile_dict = decode_anulog(anu_logfile[0], year=rough_year)
        except TypeError:
            # there was an error decoding the logfile
            # skip
            ASDF_log_file.write(anu_logfile[0] + '\t' + 'LogfileDecodeError\n')
            continue


        lat_list = logfile_dict['GPS']['LATITUDE']
        lng_list = logfile_dict['GPS']['LONGITUDE']
        alt_list = logfile_dict['GPS']['ALTITUDE']

        if lat_list == [] or lng_list == [] or alt_list == []:
            ASDF_log_file.write(anu_logfile[0] + '\t' + 'LogfileDecodeError\n')
            continue


        # remove outliers and then get mean
        def drop_outliers(x, mean, std_dev_one):
            if abs(x - mean) <= std_dev_one:
                return x


        # now avearge the lists to get coords
        av_lat = np.mean(filter(lambda x: drop_outliers(x, np.mean(lat_list), np.std(lat_list)), lat_list))
        av_lng = np.mean(filter(lambda x: drop_outliers(x, np.mean(lng_list), np.std(lng_list)), lng_list))
        av_alt = np.mean(filter(lambda x: drop_outliers(x, np.mean(alt_list), np.std(alt_list)), alt_list))

        # print(av_lat, av_lng, av_alt)

        # difference array between AUS Seismic GRID and av coordinates from logfile
        diff_array = np.absolute(lat_array - av_lat) + np.absolute(lon_array - av_lng)
        min_index = np.argmin(diff_array)

        aus_station = AUS_grid_names[min_index]
        aus_lat = AUS_grid_latitude[min_index]
        aus_lng = AUS_grid_longitude[min_index]

        found_match = False

        # if the dictionary is empty, no stations have been added in yet
        if not station_name_paras:
            new_station = aus_station
            station_name_paras[aus_station] = {'stored_lat': av_lat, 'stored_lng': av_lng, 'stored_alt': av_alt}
            station_latitude = av_lat
            station_longitude = av_lng
            station_altitude = av_alt

            found_match = True

        else:
            # check if the station has already been analysed
            for key in station_name_paras.keys():
                # print(aus_station, key, key[:-1])
                if aus_station == key[:-1] or aus_station == key:
                    # print("station in dict")
                    # check if it is within coordinates tolerance
                    # i.e. roughly within 100m same station
                    if abs(av_lat - station_name_paras[key]['stored_lat']) <= tol and abs(
                                    av_lng - station_name_paras[key]['stored_lng']) <= tol:
                        # print('station matches')
                        found_match = True
                        new_station = key
                        station_latitude = station_name_paras[key]['stored_lat']
                        station_longitude = station_name_paras[key]['stored_lng']
                        station_altitude = station_name_paras[key]['stored_alt']
                else:
                    # print("station not in dict")
                    # the station is not in the dictionary yet
                    new_station = aus_station
                    station_name_paras[aus_station] = {'stored_lat': av_lat, 'stored_lng': av_lng, 'stored_alt': av_alt}
                    station_latitude = av_lat
                    station_longitude = av_lng
                    station_altitude = av_alt

                    found_match = True

        # print(found_match)

        if not found_match:
            station_name_counter[aus_station] += 1
            new_station = aus_station + chr(ord('A') + (station_name_counter[aus_station] - 1))
            # print("no matching station", new_station)
            station_name_paras[new_station] = {'stored_lat': av_lat, 'stored_lng': av_lng, 'stored_alt': av_alt}
            station_latitude = av_lat
            station_longitude = av_lng
            station_altitude = av_alt

        # print(new_station)

        # add the logfile information into auxillary data
        data_type = "LogfileData"

        # path to logfile data in ASDF auxillary unique for each logfile for service interval
        overview_path = FDSNnetwork[0:2] + "_" + new_station + "/" + logfile_id + "/" + "Overview"
        temperature_path = FDSNnetwork[0:2] + "_" + new_station + "/" + logfile_id + "/" + "Temperature"
        lock_time_path = FDSNnetwork[0:2] + "_" + new_station + "/" + logfile_id + "/" + "LockTime"
        clock_drift_path = FDSNnetwork[0:2] + "_" + new_station + "/" + logfile_id + "/" + "ClockDrift"
        battery_path = FDSNnetwork[0:2] + "_" + new_station + "/" + logfile_id + "/" + "Battery"
        lat_path = FDSNnetwork[0:2] + "_" + new_station + "/" + logfile_id + "/" + "Latitude"
        lng_path = FDSNnetwork[0:2] + "_" + new_station + "/" + logfile_id + "/" + "Longitude"
        elev_path = FDSNnetwork[0:2] + "_" + new_station + "/" + logfile_id + "/" + "Elevation"

        # overview aux data:
        parameters = {'BSN': logfile_dict['BSN'],
                      'FWV': logfile_dict['FWV'],
                      'SPR': logfile_dict['SPR'],
                      'SMM': logfile_dict['SMM'],
                      'SMS': logfile_dict['SMS'],
                      'RCS': logfile_dict['RCS'],
                      'RCE': logfile_dict['RCE'],
                      'UDF': logfile_dict['UDF'],
                      'av_lat': av_lat, 'av_lng': av_lng, 'av_elev': av_alt}

        # add the auxillary data
        ds.add_auxiliary_data(data=np.array([0]),
                              data_type=data_type,
                              path=overview_path,
                              parameters=parameters)

        # temperature aux data
        ds.add_auxiliary_data(data=np.array(logfile_dict['GPS']['TEMPERATURE']),
                              data_type=data_type,
                              path=temperature_path,
                              parameters={"units": "Degrees Celsius"})

        # lock time auxillary data
        ds.add_auxiliary_data(data=np.array(logfile_dict['GPS']['LOCK_TIME']),
                              data_type=data_type,
                              path=lock_time_path,
                              parameters={"units": "UTC Time"})

        # clock_drift auxillary data
        ds.add_auxiliary_data(data=np.array(logfile_dict['GPS']['CLOCK']),
                              data_type=data_type,
                              path=clock_drift_path,
                              parameters={"units": "micro Seconds"})

        # battery voltage percentage auxillary data
        ds.add_auxiliary_data(data=np.array(logfile_dict['GPS']['BATTERY']),
                              data_type=data_type,
                              path=battery_path,
                              parameters={"units": "Voltage Percentage"})

        # latitude auxillary data
        ds.add_auxiliary_data(data=np.array(logfile_dict['GPS']['LATITUDE']),
                              data_type=data_type,
                              path=lat_path,
                              parameters={"units": "Decimal Degrees (Geographic)"})

        # longitude auxillary data
        ds.add_auxiliary_data(data=np.array(logfile_dict['GPS']['LONGITUDE']),
                              data_type=data_type,
                              path=lng_path,
                              parameters={"units": "Decimal Degrees (Geographic)"})

        # elevation auxillary data
        ds.add_auxiliary_data(data=np.array(logfile_dict['GPS']['ALTITUDE']),
                              data_type=data_type,
                              path=elev_path,
                              parameters={"units": "Meters"})

        # get miniseed files
        seed_files = glob.glob(join(station_path, '*miniSEED/*'))  # '*miniSEED/*.mseed*'))

        if len(seed_files) == 0:
            continue

        print '\r Working on station: ', new_station

        # dictionary for channel_location (keys) so that we can create an inventory to location level later
        channel_loc_dict = {}

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
                # Network Code: the network code in the miniseed header is prone to user error
                # (i.e. whatever the operator entered into the instrument in the field)
                orig_net = tr.stats.network
                # use the first two characters as network code. Temporary networks have start and end year as well
                new_net = FDSNnetwork[:2]
                # overwrite network code in miniseed header
                tr.stats.network = new_net
                tr.stats.station = new_station

                starttime = tr.stats.starttime.timestamp
                endtime = tr.stats.endtime.timestamp
                orig_chan = tr.stats.channel
                new_chan = orig_chan
                orig_loc = tr.stats.location
                new_loc = orig_loc

                # add channel_loc to dict
                channel_loc_dict[new_chan+'_'+new_loc] = {"samp": tr.stats.sampling_rate}

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
                    ds.add_waveforms(tr, tag="raw_recording", labels=[logfile_id, basename(service)])
                except ASDFWarning:
                    # trace already exist in ASDF file!
                    ASDF_log_file.write(filename + '\t' + ASDF_tag + '\t' + "ASDFDuplicateError\n")

                keys_list.append(str(ASDF_tag))
                info_list.append(temp_dict)

        # list for channel inventories
        channel_inventory_list = []

        for channel_loc in channel_loc_dict.keys():
            chan = channel_loc.split('_')[0]
            loc = channel_loc.split('_')[1]

            # set the azimuth and dip of the component. This is where you would record if the channels are not
            # aligned with grid north/east. i.e. from magnetometer data like OBS. Or you can change it later
            # also different sensors might use different conventions for the azimuth and dip when channels are
            # aligned properly
            if 'Z' in chan:
                az = '0'
                dip = '90'
            elif 'N' in chan:
                az = '0'
                dip = '0'
            elif 'E' in chan:
                az = '90'
                dip = '0'


            channel_inventory_list.append(inventory.Channel(code=chan, location_code=loc, depth=0, azimuth=az, dip=dip,
                              sample_rate=channel_loc_dict[channel_loc]['samp'],
                              clock_drift_in_seconds_per_sample=0, latitude=station_latitude,
                              longitude=station_longitude, elevation=station_altitude))

        # add the channels to station_inv_dict
        station_inventory_dict[new_station].append([channel_inventory_list, service])




# list of station level inventories
station_inventories_list =[]

# Now go through the station inventory dictionary and see if there are any mismatches for channel inventories
# for a station
for station, value in station_inventory_dict.iteritems():
    prev_channel_inv_set = False
    matched_channel = False
    service_mismatch_list = []
    # convert to channel list to set
    channel_inv_set = set(value[0][0])
    if prev_channel_inv_set:
        # see if any of the prev channels match
        matched_channel = bool(channel_inv_set & prev_channel_inv_set)

        if not matched_channel:
            # record which service interval does not match
            service_mismatch_list.append(value[1])

    prev_channel_inv_set = channel_inv_set

    # if some of the channels dont match then add it to the log file
    # there is a problem, i.e. sampling rates change for a given station
    if not len(service_mismatch_list) == 0:
        for service_mismatch in service_mismatch_list:
            ASDF_log_file.write("%s.%s_%s" % (FDSNnetwork, station, service_mismatch) + '\t' + "MismatchingChannels\n")


    # get the start/end dates from dict
    start_date = UTCDateTime(station_start_end_dict[station][0])
    end_date = UTCDateTime(station_start_end_dict[station][1])

    site = inventory.Site(name=station)

    # make the station_level inventory
    station = inventory.Station(code=station, creation_date=start_date, termination_date=end_date,
                                site=site,
                                latitude=value[0][0][0].latitude,
                                longitude=value[0][0][0].longitude,
                                elevation=value[0][0][0].elevation,
                                vault="Transportable Array",
                                channels=value[0][0])

    station_inventories_list.append(station)

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
                                stations=station_inventories_list)

# create the inventory
inv = inventory.Inventory(networks=[network_inv], source= "Geoscience Australia")

XML_file = join(XML_path_out, FDSNnetwork+'.xml')

if exists(XML_file):
    remove(XML_file)

# write the inventory into the default path
inv.write(path_or_file_object=XML_file, format='STATIONXML')

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
