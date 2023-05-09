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
from obspy.clients.fdsn.client import Client
import pyasdf
from mpi4py import MPI
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
from seismic.ASDFdatabase._FederatedASDFDataSetImpl import split_list
import click
from shapely.geometry.polygon import Polygon, Point
from tqdm import tqdm
import json
import tempfile

class StandardDomain():
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

class Domain():
    def __init__(self, domain_json_file:str):
        dom = json.load(open(domain_json_file, 'r'))

        coords = dom['features'][0]['geometry']['coordinates'][0]
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

    try:
        stream.merge()
    except:
        stream = Stream([])
        return stream
    # end try

    if (len(stream) > 1):
        stream = Stream([])
        return stream
    elif((stream[0].stats.endtime - stream[0].stats.starttime) < (et - st)):
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
            stream = Stream([])
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

    def filter_channels_locations(rows, preferred_location_codes:list = ['00', '', '10'],
                                  include_sband=False):
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
            net, sta, loc, cha, _, _, _ = row

            if (cha[0] == 'S' and not include_sband): continue
            if (not is_preferred_component(cha)): continue

            comp = cha[-1]
            itype_dict[net][sta][loc][comp].append([cha, irow])
        # end for

        #################################################
        # Remove multiple location codes
        #################################################
        for nk, nv in itype_dict.items():
            for sk, sv in itype_dict[nk].items():
                if(len(itype_dict[nk][sk]) > 1): # multiple location codes found
                    # first remove location codes not in preferred_location_codes
                    for lk in list(itype_dict[nk][sk].keys()):
                        if(lk not in preferred_location_codes):
                            itype_dict[nk][sk].pop(lk)
                        # end if
                    # end for

                    if(len(itype_dict[nk][sk]) > 1): # still multiple location codes found
                        for lk in list(itype_dict[nk][sk].keys()):
                            # remove all location codes other than the first entry in
                            # preferred_location_codes
                            if(lk in preferred_location_codes[1:]):
                                itype_dict[nk][sk].pop(lk)
                            # end if
                        # end for
                    # end if
                # end if
            # end for
        # end for

        # sanity checks
        for nk, nv in itype_dict.items():
            for sk, sv in itype_dict[nk].items():
                assert len(itype_dict[nk][sk]) <= 1, 'Failed to remove multiple location codes ' \
                                                     'found for station {}.{}:{}'.format(nk, sk,
                                                                                         itype_dict[nk][sk].keys())
            # end for
        # end for

        #################################################
        # Remove unwanted bands
        #################################################
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

    ods = pyasdf.ASDFDataSet(h5_ofn, mode='w', mpi=False)

    oinv = None
    otime = UTCDateTime(edict[0]['arguments']['reference_time_utc_string'])

    st, et = otime - seconds_before, otime + seconds_after
    rows = fds.get_stations(st, et)
    rows = filter_channels_locations(rows) # filter out unwanted channels and band-codes

    data_added_set = set()
    receivers_dict = defaultdict(dict)
    pbar = tqdm(total=len(rows), desc=ename)
    for row in rows:
        net, sta, loc, cha, lon, lat, elev = row

        if (not domain.contains(lon, lat)): continue
        #if (DEBUG and net not in nets): continue

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

                    data_added_set.add(seed_id)
                except Exception as e:
                    print('Failed to add inventory/waveform with error: {}. Moving along..'.format(str(e)))
                # end try
            # end if

            #if (DEBUG): break
        # end if
    # end for
    pbar.close()

    ###############################################################
    # Now download additional data from IRIS
    ###############################################################
    client = Client("IRIS")
    for net in inventory.networks:
        for sta in net.stations:
            netsta = '{}.{}'.format(net.code, sta.code)

            already_added = False
            for seed_id in data_added_set:
                if (netsta in seed_id):
                    already_added = True
                    break
                # end if
            # end for
            if(already_added): continue

            stream = []
            try:
                stream = client.get_waveforms(net.code, sta.code, '*', 'BH?,HH?',
                                              st-1, et+1)
                stream.trim(st, et)
            except:
                pass
            # end try

            if (len(stream)):

                rows = []
                for tr in stream: rows.append([tr.stats.network,
                                               tr.stats.station,
                                               tr.stats.location,
                                               tr.stats.channel,
                                               None, None, None])
                rows = filter_channels_locations(rows)

                for row in rows:
                    tr = stream.select(network=row[0],
                                       station=row[1],
                                       location=row[2],
                                       channel=row[3])
                    tr = tr[0]
                    if((tr.stats.endtime - tr.stats.starttime) < (et - st)): continue
                    print('Adding trace: {}'.format(tr.id))

                    resp = None
                    try:
                        resp = inventory.get_response(tr.id, tr.stats.starttime)
                    except Exception as e:
                        continue
                    # end try

                    if(resp):
                        try:
                            oinv = inventory.select(network=tr.stats.network,
                                                    station=tr.stats.station,
                                                    location=tr.stats.location,
                                                    channel=tr.stats.channel)
                            oinv = remove_comments(oinv)
                            ods.add_stationxml(oinv)

                            ods.add_waveforms(tr, tag)

                            receivers_dict['{}.{}.{}'.format(tr.stats.network,
                                                             tr.stats.station,
                                                             tr.stats.location)] = \
                                {"class_name": "salvus.flow.simple_config.receiver.seismology.SideSetPoint3D",
                                 "salvus_version": "0.12.8",
                                 "arguments": {
                                     "depth_in_m": 0.0,
                                     "radius_of_sphere_in_m": 6371000.0,
                                     "network_code": tr.stats.network,
                                     "location_code": tr.stats.location,
                                     "side_set_name": "r1",
                                     "latitude": sta.latitude,
                                     "longitude": sta.longitude,
                                     "station_code": tr.stats.station,
                                     "fields": [receiver_fields]}}

                            data_added_set.add(tr.id)
                        except Exception as e:
                            print('Failed to add inventory/waveform with error: {}. Moving along..'.format(str(e)))
                        # end try
                    # end if
                # end for
                #if (DEBUG): break
            # end if
        # end for
    # end for

    json.dump(receivers_dict, open(receivers_ofn, 'w+'), indent=4)
    del ods
# end func

def trim_inventory(inv:Inventory, dom:Domain):
    netsta = []

    for net in inv.networks:
        for sta in net.stations:
            if(dom.contains(sta.longitude, sta.latitude)):
                netsta.append([net.code, sta.code])
            # end if
        # end for
    # end for

    newInv = None
    for ns in netsta:
        if(newInv is None):
            newInv = inv.select(network=ns[0], station=ns[1])
        else:
            newInv += inv.select(network=ns[0], station=ns[1])
        # emd if
    # end for

    return newInv
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('asdf-source', required=True,
                type=click.Path(exists=True))
@click.argument('salvus-domain-file', required=True,
                type=click.Path(exists=True))
@click.argument('salvus-events-file', required=True,
                type=click.Path(exists=True))
@click.argument('iris_inventory', required=True,
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
            iris_inventory, output_folder,
            seconds_before, seconds_after, data_name, receiver_fields):
    """
    ASDF_SOURCE: Path to text file containing paths to ASDF files\n
    SALVUS_DOMAIN_FILE: Path to Salvus domain file in json format\n
    SALVUS_EVENTS_FILE: Path to Salvus events file in json format\n
    IRIS_INVENTORY: Path to complete channel-level IRIS inventory, containing instrument responses\n
    OUTPUT_FOLDER: Output folder \n

    """

    fds = None
    dom = None
    events = None
    inv = None

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()

    try:
        fds = FederatedASDFDataSet(asdf_source)
    except Exception as e:
        print(str(e))
        assert 0, 'Failed to instantiate FederatedASDFDataSet with input file: {}'.format(asdf_source)
    # end try

    try:
        if(rank == 0): print('Reading Salvus domain file: {}'.format(salvus_domain_file))
        dom = Domain(salvus_domain_file)
    except Exception as e:
        print(str(e))
        assert 0, 'Failed to instantiate Domain with input file: {}'.format(salvus_domain_file)
    # end try

    try:
        if(rank == 0): print('Reading Salvus events file: {}'.format(salvus_events_file))
        events = json.load(open(salvus_events_file, 'r'))
    except Exception as e:
        print(str(e))
        assert 0, 'Failed to load Salvus events file: {}'.format(salvus_events_file)
    # end try

    tempdir = None
    trimmed_inventory_fn = None
    if(rank == 0):
        try:
            print('Reading inventory file with responses: {}'.format(iris_inventory))
            inv = read_inventory(iris_inventory)
        except Exception as e:
            print(str(e))
            assert 0, 'Failed to read inventory file: {}'.format(iris_inventory)
        # end try

        # trim inventory by Domain
        inv = trim_inventory(inv, dom)

        # save trimmed inventory to a temp folder, so other ranks can load it
        tempdir = tempfile.mkdtemp()

        trimmed_inventory_fn = os.path.join(tempdir, 'trimmed_inventory.xml')
        inv.write(trimmed_inventory_fn, format="STATIONXML")
        print(trimmed_inventory_fn)
    # end if
    comm.barrier()

    # Load trimmed inventory on all ranks
    trimmed_inventory_fn = comm.bcast(trimmed_inventory_fn, root=0)
    inv = read_inventory(trimmed_inventory_fn)

    # extract data for all events
    proc_ekeys = split_list(list(events.keys()), nproc)[rank]
    for ek in tqdm(proc_ekeys, desc='Rank: {}: Events: '.format(rank)):
        extract_data_for_event(fds, dom, {ek:events[ek]}, inv, output_folder, data_name,
                               receiver_fields,
                               seconds_before, seconds_after)
        #if DEBUG: break
    # end for

    if(rank == 0):
        # copy salvus-events to output folder
        shutil.copyfile(salvus_events_file, os.path.join(output_folder, 'EVENTS/event_store.json'))
        shutil.rmtree(tempdir)
    # end if
# end func

if (__name__ == '__main__'):
    process()
# end if
