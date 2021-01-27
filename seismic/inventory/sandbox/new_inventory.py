# Ref: https://docs.obspy.org/tutorial/code_snippets/stationxml_file_from_scratch.html

import os
import sys
import obspy
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
from obspy.clients.nrl import NRL


def create_new_skeleton_inventory_file(path2xmlfile):
    """
    write a NEW skeleton inventory xml file
    :param path2xmlfile: path to a new xml file.
    :return:
    """
    # We'll first create all the various objects. These strongly follow the
    # hierarchy of StationXML files.
    inv = Inventory(
        # We'll add networks later.
        networks=[],
        # The source should be the id whoever create the file.
        source="ObsPy-Tutorial")

    net = Network(
        # This is the network code according to the SEED standard.
        code="XX",
        # A list of stations. We'll add one later.
        stations=[],
        description="A test stations.",
        # Start-and end dates are optional.
        start_date=obspy.UTCDateTime(2016, 1, 2))

    sta = Station(
        # This is the station code according to the SEED standard.
        code="ABC",
        latitude=1.0,
        longitude=2.0,
        elevation=345.0,
        creation_date=obspy.UTCDateTime(2016, 1, 2),
        site=Site(name="First station"))

    cha = Channel(
        # This is the channel code according to the SEED standard.
        code="HHZ",
        # This is the location code according to the SEED standard.
        location_code="",
        # Note that these coordinates can differ from the station coordinates.
        latitude=1.0,
        longitude=2.0,
        elevation=345.0,
        depth=10.0,
        azimuth=0.0,
        dip=-90.0,
        sample_rate=200)

    # By default this accesses the NRL online. Offline copies of the NRL can
    # also be used instead
    nrl = NRL()
    # The contents of the NRL can be explored interactively in a Python prompt,
    # see API documentation of NRL submodule:
    # http://docs.obspy.org/packages/obspy.clients.nrl.html
    # Here we assume that the end point of data logger and sensor are already
    # known:
    response = nrl.get_response( # doctest: +SKIP
        sensor_keys=['Streckeisen', 'STS-1', '360 seconds'],
        datalogger_keys=['REF TEK', 'RT 130 & 130-SMA', '1', '200'])


    # Now tie it all together.
    cha.response = response
    sta.channels.append(cha)
    net.stations.append(sta)
    inv.networks.append(net)

    # And finally write it to a StationXML file. We also force a validation against
    # the StationXML schema to ensure it produces a valid StationXML file.
    #
    # Note that it is also possible to serialize to any of the other inventory
    # output formats ObsPy supports.
    inv.write(path2xmlfile, format="stationxml", validate=True)


# Usage: this_script  path2_newfile.xml
if __name__ == "__main__":

    if len(sys.argv)<2:
        print("Usage: %s path2_newfile.xml"% sys.argv[0])
        print("Example:  python new_inventory.py new_station.xml")
        sys.exit(1)
    else:
        newfile_xml = sys.argv[1]
        if os.path.exists(newfile_xml):
            print("The input file exists already, please provide a newfile.xml")
            sys.exit(2)
        else:
            create_new_skeleton_inventory_file(newfile_xml)
