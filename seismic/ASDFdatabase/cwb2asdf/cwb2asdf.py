#!/usr/env python
"""
Description:
    Reads waveform data from a folder containing demultiplexed mseeds and
    outputs an ASDF file.
References:

CreationDate:   02/04/19
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     22/08/18   RH
    LastUpdate:     2020-04-10 FZ did clean up and test run successfully
"""

import glob
import os, fnmatch, re
import random
from collections import defaultdict

import click
import numpy as np
import pyasdf
from mpi4py import MPI
from obspy import read
from obspy.core.inventory import read_inventory
from obspy.core import Stream
from ordered_set import OrderedSet as set
from tqdm import tqdm
from seismic.misc import split_list
from seismic.misc import recursive_glob
from seismic.ASDFdatabase.utils import remove_comments

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
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
BUFFER_LENGTH = 1000

@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('input-folder', required=True,
                type=click.Path(exists=True))
@click.argument('inventory-folder', required=True,
                type=click.Path(exists=True))
@click.argument('output-file-name', required=True)
@click.option('--file-pattern', type=str, default='*.mseed',
              help="File pattern to be used while looking for data files")
@click.option('--channels-to-extract', type=str, default=None, help="Channels to extract, within quotes and space- "
                                                                   "separated.")
@click.option('--min-length-sec', type=int, default=None, help="Minimum length in seconds")
@click.option('--merge-threshold', type=int, default=None, help="Merge traces if the number of traces fetched for an "
                                                                "interval exceeds this threshold")
@click.option('--ntraces-per-file', type=int, default=3600, help="Maximum number of traces per file; if exceeded, the "
                                                                 "file is ignored.")
@click.option('--dry-run', default=False, is_flag=True, show_default=True,
              help="Dry run only reports stations that were not found in the stationXML files, for "
                   "which waveform data exists")
def process(input_folder, inventory_folder, output_file_name, file_pattern,
            channels_to_extract, min_length_sec, merge_threshold,
            ntraces_per_file, dry_run):
    """
    INPUT_FOLDER: Path to input folder containing miniseed files \n
    INVENTORY_FOLDER: Path to folder containing FDSNStationXML inventories containing
               station-level metadata for all stations \n
    OUTPUT_FILE_NAME: Name of output ASDF file \n
    """

    def _read_inventories(inventory_folder):
        inv_files = recursive_glob(inventory_folder, '*.xml')

        # Read inventories
        inv = None
        for inv_file in inv_files:
            try:
                print('Reading inventory file: {}'.format(inv_file))
                if(inv is None): inv = read_inventory(inv_file)
                else: inv += read_inventory(inv_file)
            except Exception as e:
                print(e)
                assert 0, 'Failed to read inventory file: {}'.format(inv_file)
            # end try
        # end for
        return inv
    # end func

    def _write(ds, ostream, inventory_dict, netsta_set):
        try:
            ds.add_waveforms(ostream, tag='raw_recording')
        except Exception as e:
            print(e)
            print('Failed to append stream:')
            print(ostream)
        # end try

        for item in netsta_set:
            try:
                ds.add_stationxml(inventory_dict[item])
            except Exception as e:
                print(e)
                print('Failed to append inventory:')
                print((inventory_dict[item]))
            # end try
        # end for
    # end func

    # process channels-to-extract
    if(channels_to_extract is not None):
        channels_to_extract = set([item.strip().upper() for item in channels_to_extract.split()])
    # end if

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()

    files = None
    my_files = None
    ustationInv = defaultdict(list)
    if(rank == 0):
        inv = _read_inventories(inventory_folder)

        # generate a list of files
        paths = [i for i in os.listdir(input_folder) if os.path.isfile(os.path.join(input_folder, i))]
        expr = re.compile(fnmatch.translate(file_pattern), re.IGNORECASE)
        files = [os.path.join(input_folder, j) for j in paths if re.match(expr, j)]

        files = np.array(files)
        random.Random(nproc).shuffle(files)
        #print(files); exit(0)
        #files = files[370:380]

        ustations = set()
        networklist = []
        stationlist = []
        filtered_files = []
        for file in tqdm(files, desc='Reading trace headers: '):
            #_, _, net, sta, _ = file.split('.')
            #tokens = os.path.basename(file).split('.')
            #net, sta = tokens[0], tokens[1]

            st = []
            try:
                st = read(file, headonly=True)
            except Exception as e:
                print(e)
                continue
            # end try
            if(len(st) == 0): continue

            net = st[0].meta.network
            sta = st[0].meta.station

            ustations.add('%s.%s' % (net, sta))
            networklist.append(net)
            stationlist.append(sta)
            filtered_files.append(file)
        # end for
        files = np.array(filtered_files)

        networklist = np.array(networklist)
        stationlist = np.array(stationlist)

        idx = np.lexsort((networklist, stationlist))
        files = files[idx]
        my_files = split_list(files, nproc)

        # station inventories
        for i, ustation in enumerate(ustations):
            net, sta = ustation.split('.')
            sinv = inv.select(network=net, station=sta)
            if (not len(sinv.networks)):
                print(('Missing station: %s.%s' % (net, sta)))
                ustationInv[ustation] = None
            else:
                # remove comments from inventory
                ustationInv[ustation] = remove_comments(sinv)
            # end if
        # end for
    # end if

    if(dry_run):
        # nothing more to do
        return
    # end if

    myfiles = comm.scatter(my_files, root=0)
    ustationInv = comm.bcast(ustationInv, root=0)

    # Extract trace-count-lists in parallel
    mytrccountlist = np.zeros(len(myfiles))
    for ifile, file in enumerate(tqdm(myfiles, desc='Rank: {}'.format(rank))):
        try:
            st = read(file, headonly=True)
            mytrccountlist[ifile] = len(st)
            # end if
        except Exception as e:
            print(e)
        # end try
    # end for

    trccountlist = comm.gather(mytrccountlist, root=0)

    if (rank == 0):
        trccountlist = np.array([item for sublist in trccountlist for item in sublist])

        # Some mseed files can be problematic in terms of having way too many traces --
        # e.g. 250k+ traces, each a couple of samples long, for a day mseed file. We
        # need to blacklist them and exclude them from the ASDF file.
        print(('Blacklisted %d files out of %d files; ' \
               'average trace-count %f, std: %f' % (np.sum(trccountlist > ntraces_per_file),
                                                    len(trccountlist),
                                                    np.mean(trccountlist),
                                                    np.std(trccountlist))))
        f = open(str(os.path.splitext(output_file_name)[0] + '.trccount.txt'), 'w+')
        for i in range(len(files)):
            f.write('%s\t%d\n' % (files[i], trccountlist[i]))
        # end for
        f.close()

        if (os.path.exists(output_file_name)): os.remove(output_file_name)
        ds = pyasdf.ASDFDataSet(output_file_name, compression='gzip-3', mpi=False)

        ostream = Stream()
        netsta_set = set()
        for ifile, file in enumerate(tqdm(files, desc='Rank: {}'.format(rank))):
            st = []
            if (trccountlist[ifile] > ntraces_per_file):
                continue
            else:
                try:
                    st = read(file)
                except Exception as e:
                    print(e)
                    continue
                # end try
            # end if

            if (len(st)):
                netsta = st[0].stats.network + '.' + st[0].stats.station

                if (ustationInv[netsta]):
                    # process data only if corresponding metadata exists
                    if (merge_threshold):
                        ntraces = len(st)
                        if (ntraces > merge_threshold):
                            try:
                                st = st.merge(method=1, fill_value='interpolate')
                            except:
                                print('Failed to merge traces. Moving along..')
                                continue
                            # end try
                            print(('Merging stream with %d traces' % (ntraces)))
                        # end if
                    # end if

                    if(len(ostream) < BUFFER_LENGTH):
                        for tr in st:
                            if (channels_to_extract):
                                if (tr.meta.channel not in channels_to_extract): continue
                            # end if

                            if (tr.stats.npts == 0): continue
                            if (min_length_sec):
                                if (tr.stats.npts * tr.stats.delta < min_length_sec): continue
                            # end if

                            ostream += tr
                        # end for

                        netsta_set.add(netsta)
                    else:
                        _write(ds, ostream, ustationInv, netsta_set)
                        ostream = Stream()
                        netsta_set = set()
                    # end if
                # end if
            # end if
        # end for

        _write(ds, ostream, ustationInv, netsta_set)

        print('Closing asdf file..')
        del ds
    # end if
# end func

####################################################################################################################
# Example commandline run:
# python cwb2asdf.py /tmp/Miniseed/ /Datasets/networks_fdsnstationxml/inventory.xml /Datasets/miniseeds_2asdf.h5
####################################################################################################################
if (__name__ == '__main__'):
    process()

