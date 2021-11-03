"""
Description:
    Exports an SGrid file based on the H5 output of rf_3dmigrate.py

CreationDate:   01/11/21
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     01/11/21   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

from mpi4py import MPI
from pyproj import Geod
from shapely.geometry import Point, LineString, Polygon

import logging
import json

import numpy as np
import click
import os, sys
from seismic.receiver_fn.rf_ccp_util import CCPVolume, CCP_VerticalProfile
from seismic.receiver_fn.rf_util import split_list
from collections import defaultdict

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
log = logging.getLogger('export_sgrid')

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('rf-h5-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output_sgrid_file', type=click.Path(exists=False, dir_okay=False), required=True)
def main(rf_h5_file, output_sgrid_file):
    """Export SGrid file
    """
    log.setLevel(logging.DEBUG)

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()

    geod = Geod(ellps="WGS84")

    # load ccp volume
    vol = CCPVolume(rf_h5_file)

    slon, slat = 128.1026, -25.5007
    elon, elat = 136.001816, -18.502329
    dx = dy = 2
    dz = 0.5
    ly = geod.geometry_length(LineString([Point(slon, slat),
                                          Point(slon, elat)])) / 1e3
    ny = np.int_(np.ceil(ly / dy)) + 1

    lats = np.vstack([np.array([slon, slat]),
                      np.array(geod.npts(slon, slat,
                                         slon, elat,
                                         ny-2)),
                      np.array([slon, elat])])[:,1]

    profiles = np.array([[slon, lat, elon, lat] for lat in lats])
    print(len(profiles))

    profiles = split_list(profiles, nproc)

    if(rank==0):
        print(profiles)
    # end if

    localNodes = np.array([])
    for profile in profiles[rank]:
        slice = CCP_VerticalProfile(vol, profile[0:2], profile[2:], dx=dx, dz=dz)

        if(len(localNodes) == 0):
            localNodes = np.hstack((slice._grid, slice._grid_vals.flatten()[:, None]))
        else:
            currNodes = np.hstack((slice._grid, slice._grid_vals.flatten()[:, None]))
            localNodes = np.concatenate((localNodes, currNodes), axis=0)
    # end for

    sendcount = 0
    if len(localNodes):
        print(rank, localNodes.shape)
        sendcount = np.prod(localNodes.shape)
    # end if

    sendcounts = np.array(comm.gather(sendcount, 0))

    globalNodes = np.array([])
    if(rank == 0):
        print('sendcounts: ', sendcounts)
        globalNodes = np.empty(sum(sendcounts), dtype=np.float64)
    # end if

    comm.Gatherv(sendbuf=localNodes.flatten() if len(localNodes) else np.array([]),
                 recvbuf=(globalNodes, sendcounts), root=0)
    if(rank == 0):
        globalNodes = globalNodes.reshape(-1, 4)

        np.savetxt('./nodes.txt', globalNodes)
    # end if
# end func

if __name__ == "__main__":
    # call main function
    main()
