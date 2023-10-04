from obspy import Inventory
from obspy.core.inventory import read_inventory
from obspy.core.util import AttribDict
import pandas as pd
import numpy as np
import ast
import sys

def export_clock_corrections(inv: Inventory, ofn: str):
    result = []
    for net in inv.networks:
        for sta in net.stations:
            if(hasattr(sta.extra, 'clock_corrections')):
                loc_corr_dict = ast.literal_eval(sta.extra.clock_corrections.value)
                for loc, corr_list in loc_corr_dict.items():
                    for corr in corr_list:
                        for date_range, corr_val in corr.items():
                            start_time, end_time = date_range.split(' - ')
                            result.append([net.code, sta.code, loc,
                                           start_time, end_time, corr_val])
                        # end for
                    # end for
                # end for
            # end if
        # end for
    # end for
    ds = pd.DataFrame(result, columns=['Network', 'Station', 'Location',
                                       'Start-time', 'End-time',
                                       'Correction_in_seconds'])
    ds = ds.drop_duplicates()
    with open(ofn, 'w') as f:
        f.write('# GPS clock-corrections grouped by network, station and location\n')
    # end with
    ds.to_csv(ofn, mode='a', index=False)
# end func

def export_orientation_corrections(inv: Inventory, ofn:str):
    result = []
    for net in inv.networks:
        for sta in net.stations:
            if(hasattr(sta.extra, 'rf_orientation_corrections')):
                loc_corr_dict = ast.literal_eval(sta.extra.rf_orientation_corrections.value)
                for loc, corr in loc_corr_dict.items():
                    for date_range, corr_val in corr.items():
                        start_time, end_time = date_range.split(' - ')
                        result.append([net.code, sta.code, loc, 'RF',
                                       start_time, end_time, corr_val, ''])
                    # end for
                # end for
            # end if
        # end for
    # end for

    for net in inv.networks:
        for sta in net.stations:
            if(hasattr(sta.extra, 'swp_orientation_corrections')):
                loc_corr_dict = ast.literal_eval(sta.extra.swp_orientation_corrections.value)
                for loc, corr in loc_corr_dict.items():
                    for date_range, corr in corr.items():
                        start_time, end_time = date_range.split(' - ')
                        corr_val, corr_uncert = corr
                        result.append([net.code, sta.code, loc, 'SWP',
                                       start_time, end_time,
                                       corr_val, corr_uncert])
                    # end for
                # end for
            # end if
        # end for
    # end for

    ds = pd.DataFrame(result, columns=['Network', 'Station', 'Location', 'Method',
                                       'Start-time', 'End-time',
                                       'Azimuth_correction_in_degrees',
                                       'Uncertainty±'])
    ds = ds.drop_duplicates()
    with open(ofn, 'w') as f:
        f.write('# Orientation corrections are derived from two separate methods: (i) Receiver Function (RF) '
                '(ii) Surface-wave Polarization (SWP). Only the latter method provides uncertainty estimates.\n')
    # end with
    ds.to_csv(ofn, mode='a', index=False)
# end func

if __name__=="__main__":
    if(len(sys.argv) < 3):
        print('Usage: python export_corrections.py <inventory_with_corrections.xml> <output_file_name_base>')
        exit(0)
    # end if

    inv = None
    try:
        inv = read_inventory(sys.argv[1])
    except Exception as e:
        print('Failed to read inventory with error: {}. Aborting.. '.format(str(e)))
        exit(0)
    # end try

    clock_corrs_fn = sys.argv[2] + '.clock_corrections.csv'
    orientation_corrs_fn = sys.argv[2] + '.orientation_corrections.csv'

    export_clock_corrections(inv, clock_corrs_fn)
    export_orientation_corrections(inv, orientation_corrs_fn)
# end if
