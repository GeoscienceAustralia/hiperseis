"""
Description:
    Exports contents of XDMF to CSV.

References:

CreationDate:   02/05/24
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     02/05/24   RH
"""

import click
import vtk
from vtk.numpy_interface import dataset_adapter as dsa
import numpy as np
from seismic.misc import print_exception
import pyproj
import pandas as pd
import os

def tranforms_coords(x:np.ndarray, y:np.ndarray, z:np.ndarray) -> \
    (np.ndarray, np.ndarray, np.ndarray):

    transformer = pyproj.Transformer.from_crs(
        {"proj": 'geocent', "ellps": 'WGS84', "datum": 'WGS84'},
        {"proj": 'latlong', "ellps": 'WGS84', "datum": 'WGS84'})

    xyz2lonlatalt = lambda x, y, z: np.vstack(transformer.transform(x, y, z,
                                                                    radians=False)).T

    result = xyz2lonlatalt(x, y, z)
    result[:, -1] *= -1 # convert to depth

    return result[:, 0], result[:, 1], result[:, 2] # lons, lats, depths
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('input-file', required=True,
                type=click.Path(exists=True))
@click.argument('output-folder', required=True,
                type=click.Path(exists=True))
def process(input_file, output_folder):
    reader = vtk.vtkXdmfReader()

    try:
        print('Reading input file: {}'.format(input_file))
        reader.SetFileName(input_file)
        reader.UpdateInformation()
    except Exception as e:
        print_exception(e)
    # end func


    pd2cd = vtk.vtkPointDataToCellData()
    pd2cd.SetInputConnection(reader.GetOutputPort())
    pd2cd.PassPointDataOn()

    pd2cd.Update()

    data = dsa.WrapDataObject(pd2cd.GetOutput())

    data.PointData.keys()

    xyz = data.GetPoints()

    rkeys = ['VP', 'VSH', 'VSV']
    valsDict = {}

    for key in rkeys:
        valsDict[key] = data.PointData[key]
    # end for

    # transform coordinates from geocentric xyz to geographic (both in wgs84)
    print('Transforming coordinates..')
    lons, lats, depths = tranforms_coords(xyz[:, 0], xyz[:, 1], xyz[:, 2])

    # write csv
    ofn = os.path.splitext(os.path.basename(input_file))[0] + '.csv'
    ofn = os.path.join(output_folder, ofn)

    df = pd.DataFrame()

    df['lon'] = lons
    df['lat'] = lats
    df['depth_km'] = depths
    for key in rkeys:
        df[key] = valsDict[key]
    # end for

    print('Writing output file: {}'.format(ofn))
    df.to_csv(ofn, header=True, index=False)
# end func


if __name__=="__main__":
    process()
# end if