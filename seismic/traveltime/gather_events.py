"""
Parse multiple events xml files to gather all info about seismic events, picks, arrivals, etc.

How to Run:
python seismic/traveltime/gather_events.py /g/data/ha3/events_xmls_test # a quick test
# multiple input dirs of event xml files:
python seismic/traveltime/gather_events.py /g/data/ha3/events_xmls_sc3ml /g/data/ha3/fxz547/travel_time_tomography/new_events20180516
"""
import seismic.cluster.cluster as cluster
import os
import sys

if __name__ =="__main__":

    xmldir_list = sys.argv[1:]
    xmlfiles = cluster.recursive_glob(xmldir_list)

    print ("number of events files = ", len(xmlfiles))

    print ("First file: ", xmlfiles[0])
    print ("Last file: ", xmlfiles[-1])

    cluster.gather(xmldir_list,  "out_put_file_base")