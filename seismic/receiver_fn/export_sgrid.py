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
from seismic.receiver_fn.rf_util import split_list
from collections import defaultdict
from pyproj import Proj
from seismic.receiver_fn.rf_ccp_util import CCPVolume, rtp2xyz
from collections import defaultdict
import struct

from scipy.spatial import cKDTree

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
    print('rank({}): loaded ccpVolume..'.format(rank))

    proj = Proj('+proj=utm +zone=53 +south +datum=WGS84 +units=m +no_defs')

    #OA2
    s = [129.1026-0.5, -24.5007-0.5]
    e = [135.001816+0.5, -19.502329+0.5]

    #OA1
    #s = [132.8912,  -22.506968]
    #e = [142.001816, -16]

    su = np.array(proj(s[0], s[1]))
    eu = np.array(proj(e[0], e[1]))

    ER = 6371
    MAX_DEPTH=150
    dx = 5
    dy = 5
    dz = 0.5

    xyz = None
    z = None

    nx = np.int(np.fabs(eu[0] - su[0]) / dx / 1e3) + 1
    ny = np.int(np.fabs(eu[1] - su[1]) / dy / 1e3) + 1
    nz = np.int(MAX_DEPTH / dz) + 1

    xs = np.linspace(su[0], eu[0], nx)
    ys = np.linspace(su[1], eu[1], ny)
    zs = np.linspace(ER, ER - MAX_DEPTH, nz)

    ugz, ugy, ugx = np.meshgrid(zs, ys, xs)
    ugy = ugy.transpose(1, 0, 2)
    ugz = ugz.transpose(1, 0, 2)

    if(rank == 0):
        assert np.all(ugx[0,0,:]==xs)
        assert np.all(ugy[0,:,0]==ys)
        assert np.all(ugz[:,0,0]==zs)

        ggx, ggy = proj(ugx, ugy, inverse=True)

        z = ugz.flatten()
        xyz = rtp2xyz(ugz.flatten(), np.radians(90-ggy.flatten()), np.radians(ggx.flatten()))

        xyz = split_list(xyz, nproc)
        z = split_list(z, nproc)
    # end if

    xyz = comm.scatter(xyz, root=0)
    z = comm.scatter(z, root=0)

    p = 3
    pwe = 1
    r = 40
    tree = vol._tree
    data = vol._data

    # iterate over nodes in each swath
    vals = np.zeros(xyz.shape[0])
    for i in np.arange(xyz.shape[0]):

        indices = np.array(tree.query_ball_point(xyz[i, :], r=r))

        if (len(indices) == 0): continue

        # filter out all but nodes within a disc of dz thickness
        indices = np.array(indices)
        indices = indices[np.fabs(data[indices, 3] -
                                  np.fabs((z[i] - ER))) < dz]

        if (len(indices) == 0): continue
        d = np.zeros(len(indices))

        # compute distance of kdtree nodes from current node in swath
        d[:] = np.sqrt(np.sum(np.power(xyz[i] - data[indices, :3], 2), axis=1))

        # filter out nodes outside a cone, defined as current_radius = current depth;
        # this is done to avoid lateral smearing at shallow depths, where piercing
        # points are sparse, except immediately under stations
        indices = indices[d < np.fabs((z[i] - ER))]
        d = d[d < np.fabs((z[i] - ER))]

        if (len(indices) == 0): continue

        # compute IDW weights
        idwIndices = indices
        idw = np.zeros(d.shape)
        idw = 1. / np.power(d, p)

        # compute mean instantaneous phase weight
        pw = np.mean(data[idwIndices, 5] + 1j * data[idwIndices, 6])

        # compute grid values
        v = np.sum(idw * data[idwIndices, 4]) / np.sum(idw)

        vals[i] = v * np.power(np.abs(pw), pwe)
    # end for
    print('rank({}): computed {} entries..'.format(rank, xyz.shape[0]))

    vals = comm.gather(vals, root=0)

    comm.Barrier()
    if(rank == 0):
        vals = np.concatenate(vals)


        nodataval = -9999
        prop_name='ccp_amp'
        fn = 'ccp_amp_ewns.sg'
        #prop_name='oa1_ccp_amp'
        #fn = 'oa1_ccp_amp.sg'

        ascii_data_file = fn.replace('.sg', '')+'__ascii@@'

        headerlines = [r'' + item + '\n' for item in ['GOCAD SGrid 1 ',
                                                              'HEADER {',
                                                              'name:{}'.format(prop_name),
                                                              'ascii:on',
                                                              'double_precision_binary:off',
                                                              '*painted*variable: {}'.format(prop_name),
                                                              '}',
                                                              'GOCAD_ORIGINAL_COORDINATE_SYSTEM',
                                                              'NAME Default',
                                                              'PROJECTION Unknown'
                                                              'DATUM Unknown'                                              
                                                              'AXIS_NAME "X" "Y" "Z"',
                                                              'AXIS_UNIT "m" "m" "m"',
                                                              'ZPOSITIVE Depth',
                                                              'END_ORIGINAL_COORDINATE_SYSTEM',
                                                              'AXIS_N {} {} {} '.format(
                                                                  nx, ny, nz),
                                                              'PROP_ALIGNMENT POINTS',
                                                              'ASCII_DATA_FILE {}'.format(ascii_data_file),
                                                              '',
                                                              '',
                                                              'PROPERTY 1 "{}"'.format(prop_name),
                                                              'PROPERTY_CLASS 1 "{}"'.format(prop_name),
                                                              'PROPERTY_KIND 1 "Amplitude"',
                                                              'PROPERTY_CLASS_HEADER 1 "{}" '.format(
                                                                  str.lower(prop_name)) + '{',
                                                              'low_clip:-0.5',
                                                              'high_clip:0.5',
                                                              'pclip:99',
                                                              'colormap:flag',
                                                              'last_selected_folder:Property',
                                                              '}',
                                                              'PROPERTY_SUBCLASS 1 QUANTITY Float',
                                                              'PROP_ORIGINAL_UNIT 1 arb',
                                                              'PROP_UNIT 1 arb',
                                                              'PROP_NO_DATA_VALUE 1 {}'.format(
                                                                  nodataval),
                                                              'PROP_ESIZE 1 4',
                                                              'END']]

        k, j, i = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx))
        j = j.transpose(1,0,2)
        k = k.transpose(1,0,2)

        od = np.vstack([ugx.flatten(), ugy.flatten(), (ER-ugz.flatten())*1e3,
                        vals, i.flatten(), j.flatten(), k.flatten()]).T
        datahdr = '\n X Y Z {} I J K\n'.format(prop_name)

        with open(fn, 'w') as hdrfile:
            hdrfile.writelines(headerlines)
        # end with

        np.savetxt(ascii_data_file,
            od,
            header=datahdr,
            comments='*',
            fmt=['%10.6f'] *
            4 +
            ['%10i'] *
            3)

        # output stations
        bb = np.array([[s[0], s[1]],
               [e[0], s[1]],
               [e[0], e[1]],
               [s[0], e[1]]])
        bbx, bby = proj(bb[:,0], bb[:,1])
        poly = Polygon(zip(bbx, bby))

        sta_file = fn.replace('.sg', '')+'_stations.csv'
        f=open(sta_file, 'w+')
        for k in vol._meta.keys():
            sta = vol._meta[k]
            x, y = proj(sta[0], sta[1])
            if(poly.contains(Point(x, y))):
                f.write('{}, {}, {}\n'.format(k, x, y))
            # end if
        # end for
        f.close()
    # end if
# end func

if __name__ == "__main__":
    # call main function
    main()
