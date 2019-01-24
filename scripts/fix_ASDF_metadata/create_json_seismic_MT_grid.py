import json
import numpy as np

with open('/g/data1/ha3/Passive/AUS_Seismic_MT_grid/AUS_seismic_MT_grid.txt', 'r') as f:
    data = f.readlines()[1:]


grid_dict = {}

names = []
latitude = []
longitude = []

for line in data:
    fields = line.rstrip('\r\n').split(',')
    grid_dict[fields[1]] = (float(fields[2]), float(fields[3]))

    names.append(fields[1])
    longitude.append(float(fields[2]))
    latitude.append(float(fields[3]))

lat_array = np.array(latitude)
lon_array = np.array(longitude)

test_lat = -45.2
test_lon = 122.1

diff_array = np.absolute(lat_array - test_lat) + np.absolute(lon_array - test_lon)

print(names[np.argmin(diff_array)])



