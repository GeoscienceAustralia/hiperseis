import os
from query_input_yes_no import query_yes_no
import glob
from os.path import join
import subprocess
from obspy import read_inventory
import sys
import shutil


# =========================== User Input Required =========================== #

# Path to the data
data_path = '/Users/ashbycooper/Desktop/Passive/'

# path to IRIS dataless seed to stationxml converter
seed_xml_conv_path = "/Users/ashbycooper/SEEDtoXML/stationxml-converter-1.0.9.jar"

# IRIS Virtual Ntework name
virt_net = '_AusArray'

# FDSN network identifier
FDSNnetwork = 'OA'

service_name = "service1"

# =========================================================================== #

path_DATA = join(data_path, virt_net, FDSNnetwork, 'raw_DATA/')

service_path = join(path_DATA, service_name)

#create service directory if it doesnt exist
if not os.path.exists(service_path):
    os.mkdir(service_path)



# get external SD card
sd_card = glob.glob("/Volumes/*NO NAME*")[0]


# get the dataless seed file

xseed_file = glob.glob(join(sd_card, "*.dataless"))[0].replace(" ", "\ ")


xml_out = "/Users/ashbycooper/Desktop/temp.xml"

#shell command to decode the dataless seed
decode_str = 'java -jar {0} -x -so "Geoscience Australia AusArray" -o {1} {2}'.format(seed_xml_conv_path, xml_out, xseed_file)

# decode the dataless seed file
subprocess.call(decode_str, shell=True)

# open up the recently created xml file
station_inv = read_inventory(xml_out)

network_code = station_inv[0].code
station_code = station_inv[0][0].code


print "Network Code: ", network_code
print "Station Code: ", station_code

copy_query = query_yes_no("Copy Data for " + network_code+"."+station_code + "?")

if copy_query == 'no':
    sys.exit(0)

# create directory for station
station_path = join(service_path, network_code+"_"+station_code)

if not os.path.exists(station_path):
    os.mkdir(station_path)

# copy over files
data_contents = glob.glob(join(sd_card, '*.*'))

for data_file in data_contents:
    shutil.copy(data_file, station_path)

# create miniSeed directory
station_seed_path = join(station_path, "all_miniSEED_files_are_in_here")

if not os.path.exists(station_seed_path):
    os.mkdir(station_seed_path)


def listdir_nohidden(path):
    seed_list = glob.glob(path)
    ret_list = []
    for f in seed_list:
        if not os.path.basename(f).startswith('\xe2'):
            ret_list.append(f)

    return ret_list

# get generator for all non-hidden miniseed files
ret_seed_files = listdir_nohidden(join(sd_card, "all_miniSEED_files_are_in_here", "*"))

seed_file_count = 0

for _i, seed_file in enumerate(ret_seed_files):
    print "\r     Parsing miniseed file ", _i + 1, ' of ', len(ret_seed_files), ' ....'
    sys.stdout.flush()
    shutil.copy(seed_file, station_seed_path)
    seed_file_count += 1


print "\n"
print "All Done!"
print "-------- Retrieved {0} MiniSEED files from External Card --------".format(seed_file_count)








