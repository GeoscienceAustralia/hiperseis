#!/usr/env python
"""
Description:
    Cleanses data and populates missing metadata in an ASDF file
References:

CreationDate:   24/10/2023
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     03/29/2022   RH
"""

import glob
import os
import random
from collections import defaultdict

import click
import numpy as np
import pyasdf
from obspy.core.inventory import read_inventory
from ordered_set import OrderedSet as set
from tqdm import tqdm

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('input-asdf', required=True,
                type=click.Path(exists=True))
@click.argument('inventory', required=True,
                type=click.Path(exists=True))
@click.argument('output-file-name', required=True)
@click.option('--min-length-sec', type=int, default=250, help="Minimum length in seconds")
def process(input_asdf, inventory, output_file_name, min_length_sec):
    """
    INPUT_ASDF: Path to input ASDF file \n
    INVENTORY: Path to FDSNStationXML inventory containing channel-level metadata for all stations \n
    OUTPUT_FILE_NAME: Name of output ASDF file \n
    """

    # Read inventory
    inv = None
    try:
        print('Reading inventory..')
        inv = read_inventory(inventory)
    except Exception as e:
        print(e)
    # end try

    # read asdf input file
    ids = pyasdf.ASDFDataSet(input_asdf, mode='r')
    ids.single_item_read_limit_in_mb = 20 * 1024 # todo: expose param
    netsta_list = ids.waveforms.list()

    ustations = set()
    ustationInv = defaultdict(list)
    networklist = []
    stationlist = []
    for netsta in netsta_list:
        net, sta = netsta.split('.')
        ustations.add('%s.%s' % (net, sta))
        networklist.append(net)
        stationlist.append(sta)
    # end for

    networklist = np.array(networklist)
    stationlist = np.array(stationlist)

    for i, ustation in enumerate(ustations):
        net, sta = ustation.split('.')
        sinv = inv.select(network=net, station=sta)
        if (not len(sinv.networks)):
            print(('Missing station: %s.%s' % (net, sta)))
            ustationInv[ustation] = None
        else:
            ustationInv[ustation] = sinv
        # end if
    # end for


    if (os.path.exists(output_file_name)): os.remove(output_file_name)
    ods = pyasdf.ASDFDataSet(output_file_name, compression='gzip-3', mpi=False)
    
    print('Reading and sanitizing waveforms..')
    sta_xml_added = set()
    for i, netsta in enumerate(tqdm(netsta_list)):
        net, sta = netsta.split('.')
        
        sta_acc = None
        tags = None
        try:
            sta_acc = ids.waveforms[netsta]
            tags = sta_acc.list()
        except Exception as e:
            print(netsta, str(e))
            continue
        # end try

        for tag in tags:
            if('raw_recording' not in tag): continue
            
            st = []
            try:
                st = sta_acc[tag]
            except Exception as e:
                print(tag, str(e))
            # end try

            if (len(st)):
                for tr in st:

                    if (tr.stats.npts == 0): continue
                    if (min_length_sec):
                        if (tr.stats.npts * tr.stats.delta < min_length_sec): 
                            print('Dropping trace: ', tr)
                            continue
                        # end if
                    # end if

                    try:
                        ods.add_waveforms(tr, tag='raw_recording')
                    except Exception as e:
                        print(e)
                        print('Failed to append trace:')
                        print(tr)
                    # end try
                # end for

                try:
                    if(netsta not in sta_xml_added):
                        ods.add_stationxml(ustationInv[netsta])
                        sta_xml_added.add(netsta)
                    # end if
                except Exception as e:
                    print(e)
                    print('Failed to append inventory:')
                    print((ustationInv[netsta]))
                # end try
            # end if
            #break
        # end for
        #break
    # end for
    
    print('Closing asdf file..')
    del ods
    del ids
# end func

if (__name__ == '__main__'):
    process()


