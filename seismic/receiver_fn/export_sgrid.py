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
from shapely.geometry import Point, MultiPoint, Polygon

import logging
import json

import numpy as np
import click
import os, sys
from seismic.receiver_fn.rf_util import split_list
from pyproj import Proj, Geod
from seismic.receiver_fn.rf_ccp_util import CCPVolume, rtp2xyz

def validate_coordinates(start, end, ccp_vol):
    result = []
    for i, coord in enumerate([start, end]):
        lon = lat = None
        if(coord is None):
            coords = []
            for k in ccp_vol._meta.keys():
                coords.append(ccp_vol._meta[k][:2])
            # end for

            coords = np.array(coords)
            points = MultiPoint(coords)

            if(i==0): lon, lat = points.bounds[:2]
            else: lon, lat = points.bounds[2:]
        elif(len(coord.split())==2):
            try:
                lon, lat = map(float, coord.split())
            except Exception as e:
                print(str(e))
                assert 0, "Invalid input: should 'lon lat' (within quotes, space-separated)"
            # end try
        elif(len(coord.split())==1):
            try:
                lon, lat = ccp_vol._meta[coord][:2]
            except Exception as e:
                print(str(e))
                assert 0, 'Station {} not found. The format should be NET.STA.LOC. Aborting..'.format(coord)
            # end try
        else:
            assert 0, "Invalid coordinate. Must be either NET.STA.LOC or 'lon lat' " \
                      "(within quotes, space-separated)"
        # endif

        result.append(np.array([lon, lat]))
    # end for

    return result
# end func

def get_utm_epsg(ccp_vol):
    coords = []
    for k in ccp_vol._meta.keys():
        coords.append(ccp_vol._meta[k][:2])
    # end for
    coords = np.array(coords)

    lon = np.mean(coords[:,0])
    lat = np.mean(coords[:,1])
    utm_band = str(np.int(((np.floor((lon + 180) / 6) % 60) + 1)))

    epsg_code = None
    if len(utm_band) == 1:
        utm_band = '0' + utm_band
    if lat >= 0:
        epsg_code = '326' + utm_band
    else:
        epsg_code = '327' + utm_band
    # end if

    return int(epsg_code)
# end func

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
log = logging.getLogger('export_sgrid')

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('rf-h5-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output_sgrid_file', type=click.Path(exists=False, dir_okay=False), required=True)
@click.option('--start', default=None, type=str, show_default=True,
              help="Starting station name (NET.STA.LOC) or minimum 'lon lat' (within quotes, space-separated) "
                   "for bounding box. A default value of None triggers starting coordinates to be inferred.")
@click.option('--end', default=None, type=str, show_default=True,
              help="Ending station name (NET.STA.LOC) or maximum 'lon lat' (within quotes, space-separated) "
                   "for bounding box. A default value of None triggers ending coordinates to be inferred.")
@click.option('--epsg-code', type=int, default=-1, show_default=True,
              help="EPSG code for UTM Zone to transform geodetic coordinates to. The default value triggers the "
                   "EPSG code to be inferred.")
@click.option('--extend', type=float, default=50, show_default=True,
              help="The length (km) by which to expand the bounding box, diagonally, "
                   'since CCP coverage extends laterally with depth')
@click.option('--dx', type=float, default=5, show_default=True,
              help='Longitudinal grid-cell dimension (km)')
@click.option('--dy', type=float, default=5, show_default=True,
              help='Latitudinal grid-cell dimension (km)')
@click.option('--dz', type=float, default=0.5, show_default=True,
              help='Depth-step (km)')
@click.option('--max-depth', type=click.FloatRange(0, 150), default=150, show_default=True,
              help='Maximum depth (km)')
@click.option('--cell-radius', type=float, default=40, show_default=True,
              help='CCP amplitudes at each grid element are computed through a combination of '
                   'inverse-distance- and instantaneous-phase-weighting applied to CCP values '
                   'that fall within a disk, defined by cell-radius (km) and dz. Cell-radius '
                   'should be guided by dx, dy, and areal distribution of stations in the CCP '
                   'volume')
@click.option('--idw-exponent', type=float, default=3, show_default=True,
              help='Exponent used in inverse-distance-weighted interpolation, which determines the '
                   'relative contribution of near and far values. A larger exponent diminishes the '
                   'contribution of faraway values')
@click.option('--pw-exponent', type=float, default=1, show_default=True,
              help='Exponent used in instantaneous phase-weighting of CCP amplitudes')
def process(rf_h5_file, output_sgrid_file, start, end, epsg_code, extend, dx, dy, dz, max_depth, cell_radius,
            idw_exponent, pw_exponent):
    """Export SGrid file
    RF_H5_FILE: Migrated H5 volume (output of rf_3dmigrate.py)
    OUTPUT_SGRID_FILE: Name of output file
    """

    log.setLevel(logging.DEBUG)

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()

    geod = Geod(ellps="WGS84")

    # load ccp volume
    vol = CCPVolume(rf_h5_file)

    if(rank == 0): log.info('Loaded ccpVolume..'.format(rank))

    s, e = validate_coordinates(start, end, vol)
    if(rank == 0): log.info('Grid start: {}, end: {}..'.format(s, e))

    # expand bounding box
    if(extend > 0):
        az, baz, _ = geod.inv(s[0], s[1], e[0], e[1])
        s[0], s[1], _ = geod.fwd(s[0], s[1], baz, extend * 1e3)
        e[0], e[1], _ = geod.fwd(e[0], e[1], az, extend * 1e3)

        if(rank == 0): log.info('Expanding bounding box. Grid start: {}, end: {}..'.format(s, e))
    # end if

    # infer epsg code if needed
    if(epsg_code < 0): epsg_code = get_utm_epsg(vol)

    if(rank == 0): log.info('Projecting to EPSG:{}..'.format(epsg_code))
    proj = Proj(epsg_code)

    ER = vol._earth_radius #km
    xyz = None
    z = None

    # initialize UTM grid start and end coordinates
    su = np.array(proj(s[0], s[1]))
    eu = np.array(proj(e[0], e[1]))

    # initialize grid
    nx = np.int(np.fabs(eu[0] - su[0]) / dx / 1e3) + 1
    ny = np.int(np.fabs(eu[1] - su[1]) / dy / 1e3) + 1
    nz = np.int(max_depth / dz) + 1

    xs = np.linspace(su[0], eu[0], nx)
    ys = np.linspace(su[1], eu[1], ny)
    zs = np.linspace(ER, ER - max_depth, nz)

    ugz, ugy, ugx = np.meshgrid(zs, ys, xs)
    ugy = ugy.transpose(1, 0, 2)
    ugz = ugz.transpose(1, 0, 2)

    if(rank == 0):
        assert np.all(ugx[0,0,:]==xs)
        assert np.all(ugy[0,:,0]==ys)
        assert np.all(ugz[:,0,0]==zs)

        log.info('Grid dimensions: nx:{}, ny:{} nz:{}..'.format(nx, ny, nz))

        ggx, ggy = proj(ugx, ugy, inverse=True)

        z = ugz.flatten()
        xyz = rtp2xyz(ugz.flatten(), np.radians(90-ggy.flatten()), np.radians(ggx.flatten()))

        xyz = split_list(xyz, nproc)
        z = split_list(z, nproc)
    # end if

    xyz = comm.scatter(xyz, root=0)
    z = comm.scatter(z, root=0)

    tree = vol._tree
    data = vol._data

    if(rank == 0): log.info('Computing CCP amplitudes..')

    # iterate over local nodes and compute CCP amplitude
    vals = np.zeros(xyz.shape[0])
    for i in np.arange(xyz.shape[0]):

        indices = np.array(tree.query_ball_point(xyz[i, :], r=cell_radius))

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
        idw = 1. / np.power(d, idw_exponent)

        # compute mean instantaneous phase weight
        pw = np.mean(data[idwIndices, 5] + 1j * data[idwIndices, 6])

        # compute grid values
        v = np.sum(idw * data[idwIndices, 4]) / np.sum(idw)

        vals[i] = v * np.power(np.abs(pw), pw_exponent)
    # end for

    vals = comm.gather(vals, root=0)

    comm.Barrier()

    # write SGrid file
    if(rank == 0):
        vals = np.concatenate(vals)

        nodataval = -9999
        prop_name='ccp_amp'
        fn = output_sgrid_file if (output_sgrid_file[-3:] == '.sg') else output_sgrid_file + '.sg'
        ascii_data_file = fn.replace('.sg', '')+'__ascii@@'

        log.info('Writing Sgrid file: {}..'.format(fn))

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

        np.savetxt(ascii_data_file, od,
                   header=datahdr,
                   comments='*',
                   fmt=['%10.6f'] * 4 + ['%10i'] * 3)

        # output stations
        bb = np.array([[s[0], s[1]],
               [e[0], s[1]],
               [e[0], e[1]],
               [s[0], e[1]]])
        bbx, bby = proj(bb[:,0], bb[:,1])
        poly = Polygon(zip(bbx, bby))

        sta_fn = fn.replace('.sg', '')+'_stations.csv'
        f=open(sta_fn, 'w+')
        log.info('Writing station-coordinates file: {}..'.format(sta_fn))
        for k in vol._meta.keys():
            sta = vol._meta[k]
            x, y = proj(sta[0], sta[1])
            if(poly.contains(Point(x, y))):
                f.write('{}, {}, {}\n'.format(k, x, y))
            # end if
        # end for
        f.close()

        log.info('Done..')
    # end if
# end func

if __name__ == "__main__":
    # call main function
    process()
# end if
