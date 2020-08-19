"""
Non-workflow script to convert hand-digitized data from H-k stacking at stations
into point dataset of local Moho depth.

Python package requirements:
- xlrd
"""

import os

import click
import numpy as np
import xlrd

from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
from seismic.receiver_fn.sandbox.conversion_helper import NETWORK_CODE_MAPPINGS, SPECIAL_CHARS

@click.command()
@click.option('--fds-file', type=click.Path(dir_okay=False, exists=True), required=True,
              help='Input file for FederatedASDFDataSet containing station coordinates')
@click.argument('infile', type=click.Path(dir_okay=False, exists=True), required=True)
@click.argument('sheet-names', type=str, nargs=-1)
def main(infile, fds_file, sheet_names):
    """
    Process Excel spreadsheet into point dataset based on station codes.

    Example usage:

    python hk_stations2point.py --fds-file /g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt \
        network_hk_data_sample.xlsx

    Output format is csv file containing point data in the form of lon/lat
    coordinates and depth measurement.
    For example:
        # Sta,Lon,Lat,Depth
        I8,133.035951,-19.473353,37.9
        H8,133.006100,-20.003900,45.9
        G8,132.997000,-20.486800,40.8
        F8,132.991205,-20.997177,47.3
        D8,132.989100,-21.506900,30.0
        ...

    Output file name is inferred from input Excel file name with extension changed to '.csv'

    :param infile: Input Excel file containing bespoke digitized Moho depths
    :param fds_file: Index file used to instantiate FederatedASDFDataSet
    :return: None
    """

    with xlrd.open_workbook(infile) as wb:
        if not sheet_names:
            sheet_names = wb.sheet_names()
            print('Processing all sheets:\n', sheet_names)
        else:
            _sheet_names = wb.sheet_names()
            for name in sheet_names:
                assert name in _sheet_names, 'Sheet {} not found in workbook!'.format(name)

    fds = FederatedASDFDataSet(fds_file)
    sta_coords = fds.unique_coordinates
    pts = []
    with xlrd.open_workbook(infile) as wb:
        for sheet_name in sheet_names:
            sheet = wb.sheet_by_name(sheet_name)
            print('Processing sheet {}'.format(sheet_name))
            try:
                network = sheet.cell_value(0, 0)
                network = network.split()[-1]
            except IndexError:
                print('Network code not found in string "{}", exiting'.format(network))
                exit(1)
            network = NETWORK_CODE_MAPPINGS.get(network, network)
            # end try
            for i, row in enumerate(sheet.get_rows()):
                if i == 0:
                    print('Skipping header row:', row)
                    continue
                # end if
                if not row or not row[0].value:
                    break
                # end if
                station = str(row[0].value)
                for sc in SPECIAL_CHARS:
                    station = station.split(sc)[0]
                station = '.'.join([network, station])
                h_val = float(row[1].value)
                if np.isnan(h_val):
                    print(f"Invalid depth value for {station}, skipping")
                    continue
                coords = sta_coords[station]
                if not coords:
                    print(f"Couldn't find coordinates for {station}, skipping")
                pt_data = [station] + coords + [h_val]
                pts.append(pt_data)

    all_data = np.array(pts)
    print('Collected {} samples from {} sheets'.format(all_data.shape[0], len(sheet_names)))
    filebase = os.path.splitext(infile)[0]
    outfile = filebase + '.csv'
    print('Saving point data to file "{}"'.format(outfile))
    np.savetxt(outfile, all_data, fmt=['%s', '%s', '%s', '%s'], delimiter=',',
               header='Sta,Lon,Lat,Depth')
# end func


if __name__ == '__main__':
    main()
# end if
