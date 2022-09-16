#!/bin/env python
"""
Description:
    Generate an sGrid file from tomographic inversion output

References:

CreationDate:   15/09/22
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     15/05/22   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import numpy as np
import os
from pyproj import Proj
import click
import traceback

DEFAULT_PROJ = 3577  # AU Albers

def generate_sgrid(input_file_name, output_file_name):
    fh = open(input_file_name, 'r')

    ncx = ncy = ncz = dx = dy = 0
    nx = ny = nz = 0
    xs = ys = zs = vals = None

    #===============================================================
    # Read data
    #===============================================================
    try:
        ncx, ncy, ncz, dx, dy, romi, rami, r1, r2, r3 = map(float, fh.readline().split())
        ncx, ncy, ncz = map(int, [ncx, ncy, ncz])

        cz = np.array(list(map(float, fh.readline().split())))

        nx = ncx + 1
        ny = ncy + 1
        nz = ncz + 1
        xs = np.linspace(romi - dx / 2., romi + (ncx - 1) * dx + dx / 2., nx)
        ys = np.linspace(rami - dy / 2., rami + (ncy - 1) * dy + dy / 2., ny)[::-1]
        xs[xs > 180.] -= 360

        zs = np.zeros(nz)
        for i, z in enumerate(cz):
            if (i == 0):
                zs[i + 1] = z * 2
            else:
                zs[i + 1] = zs[i] + (z - zs[i]) * 2
            # end if
        # end for

        vals = []
        for line in fh.readlines():
            vals.append(list(map(float, line.split())))
        # end for
        vals = np.array(vals).flatten()

        assert len(vals) == ncx * ncy * ncz
    except Exception as e:
        print(traceback.format_exc())
        raise RuntimeError('Failed to read file: {}..'.format(input_file_name))
    # end try
    fh.close()

    #===============================================================
    # Generate grids
    #===============================================================
    # prepare grids

    nodataval = -9999
    p = Proj(DEFAULT_PROJ)  # AU Albers

    # initialize Albers grid
    ugz, ugy, ugx = np.meshgrid(zs, ys, xs)
    ugy = ugy.transpose(1, 0, 2)
    ugz = ugz.transpose(1, 0, 2)

    assert np.all(ugx[0, 0, :] == xs)
    assert np.all(ugy[0, :, 0] == ys)
    assert np.all(ugz[:, 0, 0] == zs)

    ugx, ugy = p(ugx, ugy)

    uvals = np.ones(ugz.shape) * nodataval

    uvals[:-1, :-1, :-1] = vals.reshape((ncz, ncy, ncx))

    #===============================================================
    # Output sGrid file
    #===============================================================

    prop_name = 'perturb'
    ascii_data_file = output_file_name.replace('.sg', '')+'__ascii@@'
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
                                                          'PROP_ALIGNMENT CELLS',
                                                          'ASCII_DATA_FILE {}'.format(os.path.basename(
                                                              ascii_data_file)),
                                                          '',
                                                          '',
                                                          'PROPERTY 1 "{}"'.format(prop_name),
                                                          'PROPERTY_CLASS 1 "{}"'.format(prop_name),
                                                          'PROPERTY_KIND 1 "Amplitude"',
                                                          'PROPERTY_CLASS_HEADER 1 "{}" '.format(
                                                              str.lower(prop_name)) + '{',
                                                          #'low_clip:-0.5',
                                                          #'high_clip:0.5',
                                                          #'pclip:99',
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

    od = np.vstack([ugx.flatten(), ugy.flatten(), ugz.flatten()*1e3,
                    uvals.flatten(), i.flatten(), j.flatten(), k.flatten()]).T
    datahdr = '\n X Y Z {} I J K\n'.format(prop_name)

    with open(output_file_name, 'w') as hdrfile:
        hdrfile.writelines(headerlines)
    # end with

    np.savetxt(ascii_data_file, od,
               header=datahdr,
               comments='*',
               fmt=['%10.6f'] * 4 + ['%10i'] * 3)
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('input_file', type=click.Path(exists=True, dir_okay=False),
                required=True)
@click.argument('output_file', type=click.Path(exists=False, dir_okay=False),
                required=True)
def process(input_file, output_file):
    """
    INPUT_FILE: output of tomographic inversion
    OUTPUT_FILE_NAME: name of output file
    """

    ofn = output_file if (output_file[-3:] == '.sg') else output_file + '.sg'

    generate_sgrid(input_file, ofn)
# end func

if __name__=="__main__":
    process()
# end if
