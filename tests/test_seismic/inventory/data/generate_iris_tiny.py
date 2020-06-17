#!/usr/bin/env python

import requests as req

from obspy import read_inventory

from seismic.inventory.iris_query import form_channel_request_url

outfile = 'IRIS-ALL_tiny'

url = form_channel_request_url('AU,ES', 'M*,E*', 'BHZ')
query = req.get(url)
with open(outfile + '.xml', 'w') as f:
    f.write(query.text)

inv = read_inventory(outfile + '.xml')
inv.write(outfile + '.txt', 'stationtxt')
