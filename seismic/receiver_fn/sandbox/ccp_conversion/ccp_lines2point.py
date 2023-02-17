"""
Non-workflow script to convert hand-digitized data from CCP stacking along lines
into point dataset of local Moho depth.

Python package requirements:
- pandas
- xlrd
"""

import os

import click
import numpy as np
import pandas as pd
import xlrd

from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
from seismic.units_utils import KM_PER_DEG
from seismic.receiver_fn.plot_ccp_classic_batch import LEAD_INOUT_DIST_KM
from seismic.receiver_fn.sandbox.conversion_helper import NETWORK_CODE_MAPPINGS, SPECIAL_CHARS


@click.command()
@click.option('--fds-file', type=click.Path(dir_okay=False, exists=True), required=True,
              help='Input file for FederatedASDFDataSet containing station coordinates')
@click.option('--raise-errors', is_flag=True)
@click.argument('infile', type=click.Path(dir_okay=False, exists=True), required=True)
def main(infile, fds_file, raise_errors=False):
    """
    Process Excel spreadsheet into point dataset based on line profiles.

    Example usage:

    python ccp_lines2point.py --fds-file /g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt \
        ccp_line_data_sample.xlsx

    Output format is csv file containing point data in the form of sta, 
    lon/lat coordinates and depth measurement.
    For example:

        # Sta,Lon,Lat,Depth
        I8,134.909765,-17.572545,47.9
        H8,135.017670,-17.570829,47.3
        G8,135.134567,-17.568970,48.9
        F8,135.395337,-17.564823,52.1
        D8,135.494250,-17.563250,52.1
        ...

    Output file name is inferred from input Excel file name with extension changed to '.csv'

    :param infile: Input Excel file containing bespoke digitized Moho depths
    :param fds_file: Index file used to instantiate FederatedASDFDataSet
    :return: None
    """
    fds = FederatedASDFDataSet(fds_file)
    sta_coords = fds.unique_coordinates
    all_network_codes = {s.split('.')[0] for s in sta_coords.keys()}
    vol_data_list = []

    with xlrd.open_workbook(infile) as wb:
        sheet = wb.sheet_by_index(0)
        for i, sheet in enumerate(wb.sheets()):
            network = sheet.cell_value(0, 3)
            network = NETWORK_CODE_MAPPINGS.get(network, network)
            lines_row = sheet.row_values(3)
            lines = [line for line in lines_row if line]
            df = pd.read_excel(infile, sheet_name=i, header=4)
            df.drop(df.columns[0], axis=1, inplace=True)
    
            for i, line in enumerate(lines):
                line = line.strip()
                sta_start, sta_end = line.split(',')
                sta_start = sta_start.strip()
                sta_end = sta_end.strip()
                netsta_start = '.'.join([network, sta_start])
                netsta_end = '.'.join([network, sta_end])
                start = np.array(sta_coords[netsta_start])
                end = np.array(sta_coords[netsta_end])
                if start.size == 0 or end.size == 0:
                    msg = f"Can't get coordinates for {netsta_start} or {netsta_end}"
                    if raise_errors:
                        raise Exception(msg) 
                    else:
                        print(msg)
                        continue
                if not np.any(end != start):
                    msg =  f"Invalid profile line {netsta_start} to {netsta_end}"
                    if raise_errors:
                        raise Exception(msg) 
                    else:
                        print(msg)
                        continue
                dirn = (end - start)
                dirn = dirn / np.linalg.norm(dirn)
                dist_col = df.iloc[:, 3*i + 1]
                dist_col = pd.to_numeric(dist_col, errors='coerce').astype(float)
                valid = dist_col.notna()
                if not np.any(valid):
                    msg = f"No valid values for profile line {netsta_start} to {netsta_end}"
                    if raise_errors:
                        raise Exception(msg) 
                    else:
                        print(msg)
                        continue
                dist = dist_col[valid].values - LEAD_INOUT_DIST_KM
                depth = df.iloc[:, 3*i + 2][valid].values
                stations = df.iloc[:, 3*i][valid].values
                for sc in SPECIAL_CHARS:
                    stations = [s.split(sc)[0] for s in stations]
                stations = np.array(['.'.join([network, s]) for s in stations])
                lonlat = start + np.outer(dist, dirn)/KM_PER_DEG
                vol_data = np.hstack((stations[:, np.newaxis], lonlat, depth[:,np.newaxis]))
                vol_data_list.append(vol_data)

    filebase = os.path.splitext(infile)[0]
    outfile = filebase + '.csv'
    all_data = np.vstack(tuple(v for v in vol_data_list))
    np.savetxt(outfile, all_data, fmt=['%s', '%s', '%s', '%s'], delimiter=',',
               header='Sta,Lon,Lat,Depth')


if __name__ == '__main__':
    main()
