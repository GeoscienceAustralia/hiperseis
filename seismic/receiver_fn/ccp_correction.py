"""
Corrects CCP data by H-k data or other 1D inversion data.
"""
import json
import click
import os
from collections import defaultdict
import logging

import numpy as np


def correct(ccp, corr, outfile=None):
    """
    Corrects CCP data using 1D inversion data. The median of all a 
    station's samples are taken for both the CCP and 1D data.
     
    The difference is taken between the CCP median and the 1D inversion 
    data for the same station. This difference becomes the
    correction value, which is applied to all measurements for that
    station in the CCP data.

    Parameters
    ----------
    ccp_data : MethodDataset
        MethodDataset containg CCP data to correct.
    corr_data : str or bytes or os.PathLike
        MethodDataset containing data to correct with. Network and
        station codes need be of the same format and spelling between 
        the two files, otherwise correction won't work.
    outfile : str or bytes or os.PathLike
        The output file to save the corrected CCP CSV data file to.
        If not provided, the corrected data is saved to the working 
        directory as '{ccp_data}_corrected.csv'.
    """
    if ccp.net is None and ccp.sta is None:
        raise ValueError("Can't perform correction without network and/or station names")
    elif ccp.net is None and ccp.sta is not None:
        all_sta = set(ccp.sta)
    elif ccp.net is not None and ccp.sta is None:
        all_sta = set(ccp.net)
    else:
        all_sta = set(['.'.join((net, sta)) for net, sta in zip(ccp.net, ccp.sta)])

                #'.'.join((ccp.net, ccp.sta)))
    
    for sta in all_sta:
        ccp_med = np.median(ccp.val[ccp.sta == sta])
        corr_med = np.median(corr.val[corr.sta == sta])
        if np.isnan(corr_med):
            print(f"Not enough data to compute correction for {sta}")
            continue
        corr_value = corr_med - ccp_med
        ccp.val[ccp.sta == sta] += corr_value

    if outfile is None:
        outfile = os.path.join(os.getcwd(), os.path.splitext(os.path.basename(ccp))[0])
        outfile += '_corrected.csv'

    data = [ccp.lon, ccp.lat, ccp.val, ccp.sw]
    header = 'Lon,Lat,Depth,Weight'
    fmt = ['%s', '%s', '%s', '%s']
    if ccp.sta is not None:
        data.insert(0, ccp.sta)
        fmt.append('%s')
        header = 'Sta,' + header
    if ccp.net is not None:
        data.insert(0, ccp.net)
        fmt.append('%s')
        header = 'Net,' + header
    data = np.array(data).T

    with open(outfile, 'w') as fw:
        fw.write('# START\n')
        np.savetxt(fw, data, fmt=fmt, delimiter=',', header=header)

    print(f"Complete! Corrected data saved to '{outfile}'")
