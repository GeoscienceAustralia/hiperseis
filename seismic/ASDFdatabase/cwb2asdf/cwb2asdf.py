#!/usr/env python
"""
Description:
    Reads waveform data from a folder containing mseeds dumped from CWB server and
    dumps out an ASDF file.
References:

CreationDate:   02/04/19
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     22/08/18   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

from mpi4py import MPI
import click
import os
import pyasdf
import io
import glob
import numpy as np
from obspy import read, warnings, Stream, Trace
from obspy.core import UTCDateTime
from obspy.core.inventory import read_inventory
from collections import defaultdict
from tqdm import tqdm
import random

def split_list(lst, npartitions):
    k, m = divmod(len(lst), npartitions)
    return [lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(npartitions)]
# end func

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

@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('input-folder', required=True,
                type=click.Path(exists=True))
@click.argument('inventory', required=True,
                type=click.Path(exists=True))
@click.argument('output-file-name', required=True)
@click.option('--min-length-sec', type=int, default=None, help="Minimum length in seconds")
@click.option('--merge-threshold', type=int, default=None, help="Merge traces if the number of traces fetched for an "
                                                                "interval exceeds this threshold")
@click.option('--ntraces-per-file', type=int, default=3600, help="Maximum number of traces per file; if exceeded, the "
                                                                 "file is ignored.")
def process(input_folder, inventory, output_file_name, min_length_sec, merge_threshold, ntraces_per_file):
    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()

    # Read inventory
    inv = None
    try:
        inv = read_inventory(inventory)
    except Exception as e:
        print e
    # end try

    files = np.array(glob.glob(input_folder+'/*.mseed'))
    random.Random(10).shuffle(files)
    #files = files[:100]

    ustations = set()
    ustationInv = defaultdict(list)
    networklist = []
    stationlist = []
    for file in files:
        _, _, net, sta, _ = file.split('.')
        ustations.add('%s.%s'%(net, sta))
        networklist.append(net)
        stationlist.append(sta)
    # end for

    networklist = np.array(networklist)
    stationlist = np.array(stationlist)

    idx = np.lexsort((networklist, stationlist))
    files = files[idx]

    myfiles = split_list(files, nproc)[rank]

    if (rank == 0):
        for i, ustation in enumerate(ustations):
            net, sta = ustation.split('.')
            sinv = inv.select(network=net, station=sta)
            if(not len(sinv.networks)):
                print net + '.' + sta
                ustationInv[ustation] = None
            else:
                ustationInv[ustation] = sinv
            # end if
        # end for
    # end if
    ustationInv = comm.bcast(ustationInv, root=0)

    # Extract trace-count-lists in parallel
    mytrccountlist = np.zeros(len(myfiles))
    for ifile, file in enumerate(tqdm(myfiles)):
        try:
            st = read(file, headonly=True)
            mytrccountlist[ifile] = len(st)
            # end if
        except Exception as e:
            print e
        # end try
    # end for

    trccountlist = comm.gather(mytrccountlist, root=0)

    if (rank == 0):
        trccountlist = np.array([item for sublist in trccountlist for item in sublist])

        print 'Blacklisted %d files out of %d files; ' \
              'average trace-count %f, std: %f'%(np.sum(trccountlist>ntraces_per_file),
                                                 len(trccountlist),
                                                 np.mean(trccountlist),
                                                 np.std(trccountlist))
        f = open(str(os.path.splitext(output_file_name)[0] + '.trccount.txt'), 'w+')
        for i in range(len(files)):
            f.write('%s\t%d\n'%(files[i], trccountlist[i]))
        # end for
        f.close()

        #exit(0)
        if (os.path.exists(output_file_name)): os.remove(output_file_name)
        ds = pyasdf.ASDFDataSet(output_file_name, compression='gzip-3', mpi=False)

        for ifile, file in enumerate(tqdm(files)):
            st = []
            if(trccountlist[ifile] > ntraces_per_file):
                continue
            else:
                try:
                    st = read(file)
                except Exception as e:
                    print e
                    continue
                # end try
            # end if

            if(len(st)):
                netsta = st[0].stats.network + '.' + st[0].stats.station

                if(ustationInv[netsta]):
                    if(merge_threshold):
                        ntraces = len(st)
                        if(ntraces > merge_threshold):
                            try:
                                st = st.merge(method=1, fill_value='interpolate')
                            except:
                                print 'Failed to merge traces. Moving along..'
                                continue
                            # end try
                            print 'Merging stream with %d traces'%(ntraces)
                        # end if
                    # end if
                # end if

                for tr in st:

                    if(tr.stats.npts == 0): continue
                    if(min_length_sec):
                        if(tr.stats.npts*tr.stats.delta < min_length_sec): continue
                    # end if

                    asdfTag = make_ASDF_tag(tr, "raw_recording").encode('ascii')

                    try:
                        ds.add_waveforms(tr, tag='raw_recording')
                    except Exception as e:
                        print e
                        print 'Failed to append trace:'
                        print tr
                    # end try
                # end for

                try:
                    ds.add_stationxml(ustationInv[netsta])
                except Exception as e:
                    print e
                    print 'Failed to append inventory:'
                    print ustationInv[netsta]
                # end try
            # end if
        # end for
        print 'Closing asdf file..'
        del ds
    # end if
# end func

if (__name__ == '__main__'):
    process()
# end if
