#!/bin/env python
"""
Description:
    Small utility for exporting data from FedASDF for salvus

References:

CreationDate:   07/02/2023
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     07/02/23   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import os, shutil

from collections import defaultdict
import numpy as np
from obspy import Stream, UTCDateTime, read_inventory, Inventory
import pyasdf
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
import click
from shapely.geometry.polygon import Polygon, Point
from tqdm import tqdm
import json

class Domain():
    def __init__(self, domain_json_file:str):
        dom = json.load(open(domain_json_file, 'r'))
        lat_center, lat_extent, lon_center, lon_extent = \
            dom['lat_lon_box']['latitude_center'], \
            dom['lat_lon_box']['latitude_extent'], \
            dom['lat_lon_box']['longitude_center'], \
            dom['lat_lon_box']['longitude_extent']

        c = [lon_center, lat_center]
        e = [lon_extent, lat_extent]

        coords = ((c[0] - e[0], c[1] + e[1]),
                  (c[0] - e[0], c[1] - e[1]),
                  ((c[0] + e[0]), c[1] - e[1]),
                  ((c[0] + e[0]), c[1] + e[1]))
        self.coords = np.array(coords)
        self.bounding_polygon = Polygon(coords)
    # end func

    def contains(self, lon:float, lat:float):
        p = Point((lon, lat))

        return self.bounding_polygon.contains(p)
    # end func
# end class

def get_validated_waveform(fds: FederatedASDFDataSet,
                           net, sta, loc, cha, st, et):
    stream = fds.get_waveforms(net, sta, loc, cha,
                               st, et)
    if stream:
        try:
            stream = Stream([tr for tr in stream \
                             if tr.stats.asdf.tag == \
                             'raw_recording'])
        except Exception as e:
            pass
        # end try
    # end if

    stream.trim(st, et)
    stream.merge()

    if (len(stream) > 1):
        stream = Stream([])
        return stream
    # end if

    if any(isinstance(tr.data, np.ma.masked_array) for tr in stream):
        def has_masked_values(data_stream):
            rval = False
            for tr in data_stream:
                if (isinstance(tr.data, np.ma.masked_array)):
                    if (np.any(tr.data.mask)):
                        rval = True
                        break
                    # end if
                # end if
            # end for
            return rval

        # end func

        if (has_masked_values(stream)):
            pass
        else:
            for tr in stream: tr.data = np.array(tr.data)
        # end if
    # end for

    return stream
# end func

DEBUG = True
def extract_data_for_event(fds:FederatedASDFDataSet,
                           domain:Domain, event:dict,
                           inventory:Inventory,
                           output_folder, data_name='raw_data',
                           receiver_fields='displacement',
                           seconds_before=240, seconds_after=3600):
    """
    @param fds:
    @param domain:
    @param event:
    @param inventory:
    @param output_folder:
    @param data_name:
    @param seconds_before:
    @param seconds_after:
    @return:
    """

    def filter_channels(rows, include_sband=False):
        """
        @param rows: output of FederatedASDFDataset.get_stations
        @param include_sband:
        @return: filtered rows
        """

        def is_preferred_component(cha: str):
            pcomps = ['Z', 'N', 'E']  # preferred components

            for pc in pcomps:
                if (pc == cha[-1]): return True
            # end for

            return False

        # end func

        itype_dict = defaultdict(lambda: defaultdict( \
            lambda: defaultdict(lambda: defaultdict(list))))

        for irow, row in enumerate(rows):
            net, sta, loc, cha, lon, lat, elev = row

            if (cha[0] == 'S' and not include_sband): continue
            if (not is_preferred_component(cha)): continue

            comp = cha[-1]
            itype_dict[net][sta][loc][comp].append([cha, irow])
        # end for

        bands = ['H', 'B', 'S'] if include_sband else ['H', 'B']
        orows = []
        for nk, nv in itype_dict.items():
            for sk, sv in itype_dict[nk].items():
                for lk, lv in itype_dict[nk][sk].items():
                    for ck, cv in itype_dict[nk][sk][lk].items():

                        # More than one band-codes found for component
                        if (len(cv) > 1):
                            # select H, B or S in order of preference
                            for band in bands:
                                found = False
                                for item in cv:
                                    if (band == item[0][0]):
                                        found = True
                                        orows.append(rows[item[1]])
                                        break
                                    # end if
                                # end for
                                if (found): break
                            # end for
                            # print(nk,sk,lk, ck, cv)
                        else:
                            orows.append(rows[cv[0][1]])
                        # end if
                # end for
            # end for
        # end for

        return orows
    # end func

    def remove_comments(iinv: Inventory) -> Inventory:
        oinv = iinv.copy()
        for net in oinv.networks:
            net.comments = []
            for sta in net.stations:
                sta.comments = []
                for cha in sta.channels:
                    cha.comments = []
                # end for
            # end for
        # end for

        return oinv
    # end func

    assert len(event) == 1, 'Invalid event dict {}. Must contain a single event. Aborting..'.format(event)

    tag = 'raw_recording'
    nets = set({'AU'})

    # make folders
    ename = list(event.keys())[0]
    edict = event[ename]
    h5_dir = os.path.join(output_folder,
                            'EVENTS/{}/WAVEFORM_DATA/EXTERNAL/{}'.format(ename, data_name))
    os.makedirs(h5_dir, exist_ok=True)

    h5_ofn = os.path.join(h5_dir, 'receivers.h5')
    receivers_ofn = os.path.join(output_folder, 'EVENTS/{}/receivers.json'.format(ename))

    ods = pyasdf.ASDFDataSet(h5_ofn, mode='w')

    oinv = None
    otime = UTCDateTime(edict[0]['arguments']['reference_time_utc_string'])

    st, et = otime - seconds_before, otime + seconds_after
    rows = fds.get_stations(st, et)
    rows = filter_channels(rows) # filter out unwanted channels and band-codes

    receivers_dict = defaultdict(dict)
    pbar = tqdm(total=len(rows), desc=ename)
    for row in rows:
        net, sta, loc, cha, lon, lat, elev = row

        if (not domain.contains(lon, lat)): continue
        if (DEBUG and net not in nets): continue

        stream = get_validated_waveform(fds, net, sta, loc,
                                        cha, st, et)

        pbar.update()
        pbar.set_description(desc='{}: [{}.{}.{}.{}]'.format(ename, net, sta, loc, cha))
        if (len(stream)):
            # check instrument-response availability
            seed_id = '{}.{}.{}.{}'.format(net, sta, loc, cha)

            resp = None
            try:
                resp = inventory.get_response(seed_id, stream[0].stats.starttime)
            except Exception as e:
                continue
            # end try

            if(resp):

                try:
                    oinv = inventory.select(network=net, station=sta,
                                            location=loc, channel=cha)
                    oinv = remove_comments(oinv)
                    ods.add_stationxml(oinv)

                    ods.add_waveforms(stream, tag)

                    receivers_dict['{}.{}.{}'.format(net, sta, loc)] = \
                        {"class_name": "salvus.flow.simple_config.receiver.seismology.SideSetPoint3D",
                         "salvus_version": "0.12.8",
                         "arguments": {
                             "depth_in_m": 0.0,
                             "radius_of_sphere_in_m": 6371000.0,
                             "network_code": net,
                             "location_code": loc,
                             "side_set_name": "r1",
                             "latitude": lat,
                             "longitude": lon,
                             "station_code": sta,
                             "fields": [receiver_fields]}}
                except Exception as e:
                    print('Failed to add inventory/waveform with error: {}. Moving along..'.format(str(e)))
                # end try
            # end if

            #if (DEBUG): break
        # end if
    # end for
    pbar.close()

    json.dump(receivers_dict, open(receivers_ofn, 'w+'), indent=4)
    del ods
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('asdf-source', required=True,
                type=click.Path(exists=True))
@click.argument('salvus-domain-file', required=True,
                type=click.Path(exists=True))
@click.argument('salvus-events-file', required=True,
                type=click.Path(exists=True))
@click.argument('response-stationxml', required=True,
                type=click.Path(exists=True))
@click.argument('output-folder', required=True,
                type=click.Path(exists=True))
@click.option('--seconds-before', type=float, default=240.0, show_default=True,
              help="Start of data-window before origin-time")
@click.option('--seconds-after', type=float, default=3600.0, show_default=True,
              help="End of data-window after origin-time")
@click.option('--data-name', type=str, default='raw_data', show_default=True,
              help="Name of data folder within Salvus' expected folder hierarchy")
@click.option('--receiver-fields', type=str, default='displacement', show_default=True,
              help="Name of data folder within Salvus' expected folder hierarchy")
def process(asdf_source, salvus_domain_file, salvus_events_file,
            response_stationxml, output_folder,
            seconds_before, seconds_after, data_name, receiver_fields):
    """
    ASDF_SOURCE: Path to text file containing paths to ASDF files\n
    SALVUS_DOMAIN_FILE: Path to Salvus domain file in json format
    SALVUS_EVENTS_FILE: Path to Salvus events file in json format
    RESPONSE_STATIONXML: Path to stationXML file containing instrument responses
    OUTPUT_FOLDER: Output folder \n

    """

    fds = None
    dom = None
    events = None
    inv = None

    try:
        fds = FederatedASDFDataSet(asdf_source)
    except Exception as e:
        print(str(e))
        assert 0, 'Failed to instantiate FederatedASDFDataSet with input file: {}'.format(asdf_source)
    # end try

    try:
        print('Reading Salvus domain file: {}'.format(salvus_domain_file))
        dom = Domain(salvus_domain_file)
    except Exception as e:
        print(str(e))
        assert 0, 'Failed to instantiate Domain with input file: {}'.format(salvus_domain_file)
    # end try

    try:
        print('Reading Salvus events file: {}'.format(salvus_events_file))
        events = json.load(open(salvus_events_file, 'r'))
    except Exception as e:
        print(str(e))
        assert 0, 'Failed to load Salvus events file: {}'.format(salvus_events_file)
    # end try

    try:
        print('Reading inventory file with responses: {}'.format(response_stationxml))
        inv = read_inventory(response_stationxml)
    except Exception as e:
        print(str(e))
        assert 0, 'Failed to read inventory file: {}'.format(response_stationxml)
    # end try

    # extract data for all events
    for ek, e in tqdm(events.items(), desc='Events: '):
        extract_data_for_event(fds, dom, {ek:e}, inv, output_folder, data_name,
                               receiver_fields,
                               seconds_before, seconds_after)
        #if DEBUG: break
    # end for


    # copy salvus-events to output folder
    shutil.copyfile(salvus_events_file, os.path.join(output_folder, 'EVENTS/event_store.json'))
# end func

if (__name__ == '__main__'):
    process()
# end if