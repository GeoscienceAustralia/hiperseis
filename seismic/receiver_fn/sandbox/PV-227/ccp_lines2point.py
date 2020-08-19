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
from seismic.receiver_fn.plot_ccp_batch import LEAD_INOUT_DIST_KM


@click.command()
@click.option('--fds-file', type=click.Path(dir_okay=False, exists=True), required=True,
              help='Input file for FederatedASDFDataSet containing station coordinates')
@click.argument('infile', type=click.Path(dir_okay=False, exists=True), required=True)
def main(infile, fds_file):
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

    with xlrd.open_workbook(infile) as wb:
        sheet = wb.sheet_by_index(0)
        network = sheet.cell_value(0, 3)
        lines_row = sheet.row_values(3)
        lines = [line for line in lines_row if line]
    # end with

    df = pd.read_excel(infile, header=4)
    df.drop(df.columns[0], axis=1, inplace=True)

    fds = FederatedASDFDataSet(fds_file)
    sta_coords = fds.unique_coordinates
    vol_data_dict = {}
    for i, line in enumerate(lines):
        line = line.strip()
        sta_start, sta_end = line.split(',')
        sta_start = sta_start.strip()
        sta_end = sta_end.strip()
        start = '.'.join([network, sta_start])
        end = '.'.join([network, sta_end])
        start = np.array(sta_coords[start])
        end = np.array(sta_coords[end])
        assert np.any(end != start)
        dirn = (end - start)
        dirn = dirn / np.linalg.norm(dirn)
        dist_col = df.iloc[:, 3*i + 1]
        dist_col = pd.to_numeric(dist_col, errors='coerce').astype(float)
        valid = dist_col.notna()
        if not np.any(valid):
            continue
        dist = dist_col[valid].values - LEAD_INOUT_DIST_KM
        depth = df.iloc[:, 3*i + 2][valid].values
        stations = df.iloc[:, 3*i].values
        stations = [s.strip().split('.')[0] for s in stations]
        lonlat = start + np.outer(dist, dirn)/KM_PER_DEG
        # Difficult to correct for differences in station elevation because
        # FDS does not include it in station coords. Ignore for now.
        vol_data = np.hstack((stations, lonlat, depth[:, np.newaxis]))
        vol_data_dict[line] = vol_data
    # end for

    filebase = os.path.splitext(infile)[0]
    outfile = filebase + '.csv'
    all_data = np.vstack(tuple(v for v in vol_data_dict.values()))
    np.savetxt(outfile, all_data, fmt=['%s', '%.6f', '%.6f', '%.1f'], delimiter=',',
               header='Sta,Lon,Lat,Depth')
# end func


if __name__ == '__main__':
    main()
# end if
