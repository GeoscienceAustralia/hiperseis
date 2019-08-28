#!/usr/env python
"""
Description:
    This script reads in parallel xml files containing events in SC3ML format
    and outputs two txt files, one for p- and the other for s-arrivals, containing
    the same columns as output by sc3_extract_events.py.

References:

CreationDate:   28/08/19
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     28/08/19   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import sys
import click
from obspy import Stream, Trace, UTCDateTime
import subprocess
from tqdm import tqdm
from mpi4py import MPI
from collections import defaultdict
import numpy as np
from random import shuffle
import os, glob, fnmatch, sys
from obspy import read_events, read_inventory
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees
import MySQLdb
import os

def split_list(lst, npartitions):
    k, m = divmod(len(lst), npartitions)
    return [lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(npartitions)]
# end func

def recursive_glob(treeroot, pattern):
    results = []
    for base, dirs, files in os.walk(treeroot):
        goodfiles = fnmatch.filter(files, pattern)
        results.extend(os.path.join(base, f) for f in goodfiles)
    return results
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('data-path', required=True,
                type=str)
@click.argument('inventory-file', required=True,
                type=click.Path(exists=True))
@click.argument('scratch-path', required=True,
                type=str)
@click.argument('output-file-stem', required=True,
                type=click.Path(exists=False))
def process(data_path, inventory_file, scratch_path, output_file_stem):
    """
    DATA_PATH: input-folder
    Inventory: FDSN inventory
    SCRATCH_PATH: scratch-folder
    OUTPUT_FILE_STEM: output file stem
    """

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()

    # Fetch and distribute xml files
    inventory = defaultdict(list)
    files = []
    if (rank==0):
        files = recursive_glob(data_path, "*.xml")
        files = split_list(files, nproc)

        inv = read_inventory(inventory_file)
        for n in inv.networks:
            for s in n.stations:
                inventory['%s.%s'%(n.code, s.code)] = [s.longitude, s.latitude, s.elevation]
            # end for
        # end for
    # end if

    inventory = comm.bcast(inventory, root=0)
    files = comm.scatter(files, root=0)

    pprocfile = open('%s/pproc.%d.txt'%(scratch_path, rank), 'w+')
    sprocfile = open('%s/sproc.%d.txt'%(scratch_path, rank), 'w+')
    for file in tqdm(files):
        cat = read_events(file, format='SC3ML')
        notFound = defaultdict(int)

        linesp = []
        liness = []
        for e in cat.events:
            po = e.preferred_origin()

            # retrieve depth; some preferred origins don't have depth values
            poDepth = po.depth
            if(poDepth == None):
                for o in e.origins:
                    if(o.depth): poDepth = o.depth
                # end for
            # end if

            if(poDepth == None): continue

            for a in po.arrivals:
                try:
                    ncode = a.pick_id.get_referred_object().waveform_id.network_code
                    scode = a.pick_id.get_referred_object().waveform_id.station_code
                    ccode = a.pick_id.get_referred_object().waveform_id.channel_code
                except:
                    continue


                slon = None
                slat = None
                try:
                    slon, slat = inventory[ncode][scode]
                except:
                    notFound['%s.%s'%(ncode, scode)] += 1
                    continue
                # end try

                # get band-index and snr from comments
                pick_attribs = defaultdict(lambda:-999.)
                pick = a.pick_id.get_referred_object()
                for c in pick.comments:
                    if ('text' in c.keys()):
                        for item in c['text'].split(','):
                            k,v = item.split('=')
                            pick_attribs[k.strip()] = float(v)
                        # end for
                    # end if
                # end for

                # get az, baz and distance
                da = gps2dist_azimuth(po.latitude, po.longitude, slat, slon)

                # create row
                if(a.phase not in ['P', 'S']): continue
                if(a.time_residual is None): continue

                line = [e.resource_id, '{:<25s}',
                        po.time.timestamp, '{:f}',
                        e.magnitudes[0].mag if (len(e.magnitudes)) else -999, '{:f}',
                        po.longitude, '{:f}',
                        po.latitude, '{:f}',
                        poDepth/1e3, '{:f}',
                        ncode, '{:<5s}',
                        scode, '{:<5s}',
                        ccode, '{:<5s}',
                        pick.time.timestamp, '{:f}',
                        a.phase, '{:<5s}',
                        slon, '{:f}',
                        slat, '{:f}',
                        da[1], '{:f}',
                        da[2], '{:f}',
                        kilometers2degrees(da[0]/1e3), '{:f}',
                        a.time_residual, '{:f}',
                        pick_attribs['phasepapy_snr'], '{:f}',
                        pick_attribs['quality_measure_cwt'], '{:f}',
                        pick_attribs['dom_freq'], '{:f}',
                        pick_attribs['quality_measure_slope'], '{:f}',
                        int(pick_attribs['band_index']), '{:d}',
                        int(pick_attribs['nsigma']), '{:d}']

                if(a.phase == 'P'): linesp.append(line)
                elif(a.phase == 'S'): liness.append(line)
            # end for
        # end for

        for line in linesp:
            lineout = ' '.join(line[1::2]).format(*line[::2])
            pprocfile.write(lineout + '\n')
        # end for

        for line in liness:
            lineout = ' '.join(line[1::2]).format(*line[::2])
            sprocfile.write(lineout + '\n')
        # end for
        if (len(notFound)): print 'Rank: %d'%(rank), notFound
    # end for

    pprocfile.close()
    sprocfile.close()

    header = '#eventID originTimestamp mag originLon originLat originDepthKm net sta cha pickTimestamp phase stationLon stationLat az baz distance ttResidual snr qualityMeasureCWT domFreq qualityMeasureSlope bandIndex nSigma\n'
    comm.barrier()
    if(rank == 0):
        ofp = open(output_file_stem + '.p.txt', 'w+')
        ofs = open(output_file_stem + '.s.txt', 'w+')

        ofp.write(header)
        ofs.write(header)
        for i in range(nproc):
            pfn = '%s/pproc.%d.txt'%(scratch_path, i)
            sfn = '%s/sproc.%d.txt'%(scratch_path, i)

            lines = open(pfn, 'r').readlines()
            for line in lines:
                ofp.write(line)
            # end for

            lines = open(sfn, 'r').readlines()
            for line in lines:
                ofs.write(line)
            # end for

            if (os.path.exists(pfn)): os.remove(pfn)
            if (os.path.exists(sfn)): os.remove(sfn)
        # end for
        ofp.close()
        ofs.close()
    # end if
# end func

if __name__=="__main__":
    process()
# end if
