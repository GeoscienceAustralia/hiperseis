#!/usr/bin/env python

import os

import numpy as np
import click
import rf

import seismic.receiver_fn.rf_util as rf_util
from seismic.ASDFdatabase import FederatedASDFDataSet
from seismic.receiver_fn.plot_ccp import run


def run_batch(transect_file, rf_waveform_file, fed_db_file, amplitude_filter=False, similarity_filter=False,
              stack_scale=0.4, width=30.0, spacing=2.0, max_depth=200.0,
              channels='R', output_folder='', colormap='seismic'):
    """Run CCP generation in batch mode along a series of transects.

    :param transect_file: File containing specification of network and station locations of ends of transects
    :type transect_file: str or Path
    :param rf_waveform_file: HDF5 file of QA'd receiver functions for the network matching the transect file
    :type rf_waveform_file: str or Path
    :param fed_db_file: Name of file with which to initialize FederatedASDFDataBase
    :type fed_db_file: str or Path
    :param amplitude_filter: Whether to use amplitude-based filtering of waveforms beform plotting.
    :type amplitude_filter: bool
    :param similarity_filter: Whether to use RF waveform similarity filtering of waveforms beform plotting.
    :type similarity_filter: bool
    :param stack_scale: Max value to represent on color scale of CCP plot
    :type stack_scale: float
    :param width: Width of transect (km)
    :type width: float
    :param spacing: Discretization size (km) for RF ray sampling
    :type spacing: float
    :param max_depth: Maximum depth of slice below the transect line (km)
    :type max_depth: float
    :param channels: String of comma-separated component IDs to source for the RF amplitude
    :type channels: str, comma separated
    :return: None
    """

    print("Reading HDF5 file...")
    rf_stream = rf.read_rf(rf_waveform_file, 'H5')

    rf_type = 'ZRT'
    if amplitude_filter:
        # Label and filter quality
        rf_util.label_rf_quality_simple_amplitude(rf_type, rf_stream)
        rf_stream = rf.RFStream([tr for tr in rf_stream if tr.stats.predicted_quality == 'a'])
    # end if
    if similarity_filter and len(rf_stream) >= 3:
        rf_stream = rf_util.filter_crosscorr_coeff(rf_stream)
    # end if

    db = FederatedASDFDataSet.FederatedASDFDataSet(fed_db_file)
    sta_coords = db.unique_coordinates

    if output_folder and not os.path.isdir(output_folder):
        assert not os.path.isfile(output_folder)
        os.makedirs(output_folder, exist_ok=True)
    # end if

    with open(transect_file, 'r') as f:
        net = f.readline().strip()
        for transect in f.readlines():
            if not transect.strip():
                continue
            sta_start, sta_end = transect.split(',')
            sta_start = sta_start.strip()
            sta_end = sta_end.strip()
            start = '.'.join([net, sta_start])
            end = '.'.join([net, sta_end])
            start = np.array(sta_coords[start])
            end = np.array(sta_coords[end])
            # Offset ends slightly to make sure we don't lose end stations due to truncation error.
            # Note: for simplicity this treats lat/lon like cartesian coords, but this is approximate
            # and will break down near poles, for long transects, or if transect crosses the antimeridian.
            dirn = (end - start)
            dirn = dirn/np.linalg.norm(dirn)
            start -= 25*dirn/rf_util.KM_PER_DEG
            end += 25*dirn/rf_util.KM_PER_DEG
            start_latlon = (start[1], start[0])
            end_latlon = (end[1], end[0])

            outfile = '{}-ZRT-R_CCP_stack_{}-{}_{}km_spacing.png'.format(net, sta_start, sta_end, spacing)
            title = 'Network {} CCP R-stacking (profile {}-{})'.format(net, sta_start, sta_end)

            outfile = os.path.join(output_folder, outfile)
            run(rf_stream, outfile, start_latlon, end_latlon, width, spacing, max_depth, channels,
                stacked_scale=stack_scale, title=title, colormap=colormap)

    # end for
    # end with
# end func


@click.command()
@click.option('--rf-file', type=click.Path(exists=True, dir_okay=False), required=True,
              help='HDF5 file containing receiver functions')
@click.option('--waveform-database', type=click.Path(exists=True, dir_okay=False), required=True,
              help='Location of waveform source database used to generate FederatedASDFDataSet. '
                   'Provides station location metadata. e.g. "/g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt".')
@click.option('--stack-scale', type=float, default=0.4, show_default=True,
              help='Max value to represent on color scale of CCP plot')
@click.option('--width', type=float, default=40.0, show_default=True,
              help='Width of transect (km)')
@click.option('--depth', type=float, default=200.0, show_default=True,
              help='Depth of slice below the transect line (km)')
@click.option('--spacing', type=float, default=2.0, show_default=True,
              help='Discretization size (km) of grid beneath transect line')
@click.option('--channel', type=click.Choice(['R', 'T', 'Q']), default='R', show_default=True,
              help='Channel to use for stacking, e.g. R')
@click.option('--colormap', type=str, default='seismic', show_default=True,
              help='Colormap to use. Must be recognized by matplotlib. Suggest seismic, coolwarm or jet.')
@click.option('--apply-amplitude-filter', is_flag=True, default=False, show_default=True,
              help='Apply RF amplitude filtering to the RFs.')
@click.option('--apply-similarity-filter', is_flag=True, default=False, show_default=True,
              help='Apply RF similarity filtering to the RFs.')
@click.argument('transect-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output-folder', type=click.Path(dir_okay=True, file_okay=False), required=False, default='')
def main(transect_file, output_folder, rf_file, waveform_database, stack_scale, width, depth, spacing, channel,
         colormap, apply_amplitude_filter, apply_similarity_filter):
    """
    Batch mode generation of CCP stacks.

    transect_file: File containing network code (first line) followed by transect line start/end points in
    the form of pairs on comma-separated station codes, two per line

    output_folder: Folder in which to place output files.
    """
    run_batch(transect_file, rf_file, waveform_database, stack_scale=stack_scale, width=width, spacing=spacing,
              max_depth=depth, channels=channel, output_folder=output_folder, colormap=colormap,
              amplitude_filter=apply_amplitude_filter, similarity_filter=apply_similarity_filter)
# end main


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
# end if
