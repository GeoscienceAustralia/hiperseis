"""
Requirements:
- pandas
- xlrd
"""

import os

import numpy as np
import pandas as pd
import xlrd

from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
from seismic.units_utils import KM_PER_DEG
from seismic.receiver_fn.plot_ccp_batch import LEAD_INOUT_DIST_KM


def main(infile, fds_file):

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
        lonlat = start + np.outer(dist, dirn)/KM_PER_DEG
        # Difficult to correct for differences in station elevation because
        # FDS does not include it in station coords. Ignore for now.
        vol_data = np.hstack((lonlat, depth[:, np.newaxis]))
        vol_data_dict[line] = vol_data
    # end for

    filebase = os.path.splitext(infile)[0]
    outfile = filebase + '.csv'
    all_data = np.vstack(tuple(v for v in vol_data_dict.values()))
    np.savetxt(outfile, all_data, fmt=['%.6f', '%.6f', '%.1f'], delimiter=',',
               header='Lon,Lat,Depth')
# end func


if __name__ == '__main__':
    main('ccp_line_data_sample.xlsx',
         '/g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt')
# end if
