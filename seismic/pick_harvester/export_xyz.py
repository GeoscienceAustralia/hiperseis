#!/bin/env python
"""
Description:
    Generate an xyz file from tomographic inversion output

References:

CreationDate:   01/06/23
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     01/06/23   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import numpy as np
import click
import traceback


def generate_xyz(input_file_name, output_file_name):
    fh = open(input_file_name, 'r')

    ncx = ncy = ncz = dx = dy = 0
    xs = ys = zs = vals = None

    #===============================================================
    # Read data
    #===============================================================
    try:
        ncx, ncy, ncz, dx, dy, romi, rami, r1, r2, r3 = map(float, fh.readline().split())
        ncx, ncy, ncz = map(int, [ncx, ncy, ncz])

        cz = np.array(list(map(float, fh.readline().split())))

        cx = np.linspace(romi + dx/2, romi + ncx * dx - dx/2, ncx)
        cy = np.linspace(rami + dy/2, rami + ncy * dy - dy/2, ncy)[::-1]
        cx[cx > 180.] -= 360

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

    gz, gy, gx = np.meshgrid(cz, cy, cx)
    gy = gy.transpose(1, 0, 2)
    gz = gz.transpose(1, 0, 2)

    assert np.all(gx[0, 0, :] == cx)
    assert np.all(gy[0, :, 0] == cy)
    assert np.all(gz[:, 0, 0] == cz)

    gvals = np.zeros(gz.shape)

    gvals[:, :, :] = vals.reshape((ncz, ncy, ncx))

    np.savetxt(output_file_name,
               np.array([gx.flatten(), gy.flatten(), gz.flatten(), gvals.flatten()]).T,
               delimiter=',', header='lon, lat, depth(km), perturbation')
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
    OUTPUT_FILE_NAME: name of output file, with columns: lon, lat, depth(km), perturbation
    """

    ofn = output_file if (output_file[-4:] == '.csv') else output_file + '.csv'

    generate_xyz(input_file, ofn)
# end func

if __name__=="__main__":
    process()
# end if
