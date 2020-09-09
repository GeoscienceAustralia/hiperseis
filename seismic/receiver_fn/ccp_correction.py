"""
Corrects CCP data by H-k data or other 1D inversion data.
"""
import json
import click
import os
from collections import defaultdict
import logging

import numpy as np


def correct(ccp_data, corr_data, outfile=None):
    """
    Corrects CCP data using 1D inversion data. The median of all a 
    station's samples are taken for both the CCP and 1D data.
     
    The difference is taken between the CCP median and the 1D inversion 
    data for the same station. This difference becomes the
    correction value, which is applied to all measurements for that
    station in the CCP data.

    Parameters
    ----------
    ccp_data : str or bytes or os.PathLike
        The CSV file containg CCP data to correct. Needs to be of 
        format `STA,LON,LAT,DEPTH,SAMPLE_WEIGHT`.
    corr_data : str or bytes or os.PathLike
        The CSV file containing 1D data to correct by. Needs to be of
        format `STA,LON,LAT,DEPTH,SAMPLE_WEIGHT`. Station codes `STA`
        should be of the same format and spelling between the two files.
    outfile : str or bytes or os.PathLike
        The output file to save the corrected CCP CSV data file to.
        If not provided, the corrected data is saved to the working 
        directory as '{ccp_data}_corrected.csv'.
    """
    print(f"Correcting {ccp_data} by {corr_data}")
    def _load(fname):
        return np.genfromtxt(fname, delimiter=',', dtype=None, encoding=None,
                             names=['sta', 'lon', 'lat', 'depth', 'weight'])
    
    ccp = _load(ccp_data)
    corr = _load(corr_data)

    all_sta = set(ccp['sta'])

    for sta in all_sta:
        ccp_med = np.median(ccp['depth'][ccp['sta'] == sta])
        corr_med = np.median(corr['depth'][corr['sta'] == sta])
        if np.isnan(corr_med):
            print(f"Not enough data to compute correction for {sta}")
            continue
        corr_value = corr_med - ccp_med
        ccp['depth'][ccp['sta'] == sta] += corr_value

    if outfile is None:
        outfile = os.path.join(os.getcwd(), os.path.splitext(os.path.basename(ccp_data))[0])
        outfile += '_corrected.csv'

    with open(outfile, 'w') as fw:
        np.savetxt(fw, ccp, fmt=['%s', '%s', '%s', '%s', '%s'], delimiter=',',
                   header='Sta,Lon,Lat,Depth,Weight')

    print(f"Complete! Corrected data saved to '{outfile}'")


@click.command()
@click.option('--outfile', type=click.Path(dir_okay=False, exists=False),
              required=False)
@click.argument('ccp-data', type=click.Path(dir_okay=False, exists=True),
                required=True)
@click.argument('correction-data', type=click.Path(dir_okay=False, exists=True),
                required=True)
def main(ccp_data, correction_data, outfile=None):
    correct(ccp_data, correction_data, outfile)


if __name__ == '__main__':
    main()
