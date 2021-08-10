#!/usr/bin/env python
# coding: utf-8
"""Produce PDF report of network stations showing RF waveforms
"""

from mpi4py import MPI
import os
import re
import logging
import itertools

import numpy as np
import click

import rf
import rf.imaging

import tqdm.auto as tqdm

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import pandas as pd
import seismic.receiver_fn.rf_util as rf_util
import seismic.receiver_fn.rf_corrections as rf_corrections
import seismic.receiver_fn.rf_plot_utils as rf_plot_utils
import seismic.receiver_fn.rf_stacking as rf_stacking

# pylint: disable=invalid-name, logging-format-interpolation, too-many-arguments, too-many-statements, too-many-locals
from seismic.receiver_fn.rf_plot_utils import pdf_merge

logging.basicConfig()

paper_size_A4 = (8.27, 11.69)  # inches

DEFAULT_HK_SOLN_LABEL = 'global'

def _get_aspect(ax):
    """Compute aspect ratio of given axes data."""
    from operator import sub
    # Total figure size
    figW, figH = ax.get_figure().get_size_inches()
    # Axis size on figure
    _, _, w, h = ax.get_position().bounds
    # Ratio of display units
    disp_ratio = (figH * h) / (figW * w)
    # Ratio of data units
    # Negative over negative because of the order of subtraction
    data_ratio = sub(*ax.get_ylim()) / sub(*ax.get_xlim())

    return disp_ratio / data_ratio
# end func


def _rf_layout_A4(fig):
    """Layout plots for A4 paper size"""
    # Fix aspect ratio of stack plot
    ax = fig.axes[3]
    ax.set_aspect(_get_aspect(ax))

    # Fix aspect ratio of traces plot
    ax = fig.axes[0]
    ax.set_aspect(_get_aspect(ax))

    # Set to A4 paper size
    fig.set_size_inches(*paper_size_A4)

    # Adjust position of stack plot to fixed location on page
    ax3 = fig.axes[3]
    ax3_pos = ax3.get_position()
    ax3_top = 0.95
    ax3.set_position([ax3_pos.x0, ax3_top - ax3_pos.height, ax3_pos.width, ax3_pos.height])
    ax3.set_anchor('NW')
    ax3_pos = ax3.get_position()

    # Adjust position of traces plot to pack at top of page with against stack plot
    ax0 = fig.axes[0]
    ax0_pos = ax0.get_position()
    ax0.set_position([ax0_pos.x0, ax3_pos.y0 - 0.01 - ax0_pos.height, ax0_pos.width, ax0_pos.height])
    ax0.set_anchor('NW')
    ax0_pos = ax0.get_position()

    # Adjust y-position of distance/azimuth marker plots to match traces plot
    ax1 = fig.axes[1]
    ax2 = fig.axes[2]
    ax1_pos = ax1.get_position()
    ax2_pos = ax2.get_position()
    ax1.set_position([ax1_pos.x0, ax0_pos.y0, ax1_pos.width, ax0_pos.height])
    ax2.set_position([ax2_pos.x0, ax0_pos.y0, ax2_pos.width, ax0_pos.height])
# end func


def _produce_hk_stacking(channel_data, weighting=rf_stacking.DEFAULT_WEIGHTS,
                         labelling=DEFAULT_HK_SOLN_LABEL, depth_colour_range=(20, 70)):
    """Helper function to produce H-k stacking figure."""

    k_grid, h_grid, hk_stack = rf_stacking.compute_hk_stack(channel_data,
                                                            h_range=rf_stacking.DEFAULT_H_RANGE,
                                                            k_range=rf_stacking.DEFAULT_k_RANGE,
                                                            weights=weighting,
                                                            root_order=2)

    sta = channel_data[0].stats.station
    loc = channel_data[0].stats.location
    channel = channel_data[0].stats.channel
    num = len(channel_data)
    title = '.'.join([sta, loc, channel])

    fig = rf_plot_utils.plot_hk_stack(k_grid, h_grid, hk_stack, title=title, num=num,
                                      depth_colour_range=depth_colour_range,
                                      stack_ylabel='Moho depth')

    # Stamp weightings onto plot
    xl = plt.xlim()
    yl = plt.ylim()
    txt_x = xl[0] + 0.95*(xl[1] - xl[0])
    txt_y = yl[0] + 0.90*(yl[1] - yl[0])
    plt.text(txt_x, txt_y, "w={}".format(weighting), horizontalalignment='right', color="#ffffff",
             fontsize=12, rasterized=True)

    # Find and label location of maximum
    soln = []
    if labelling == 'global':
        h_max, k_max = rf_stacking.find_global_hk_maximum(k_grid, h_grid, hk_stack)
        soln = [(h_max, k_max)]
        log.info("Numerical solution (H, k) = ({:.3f}, {:.3f})".format(*soln[0]))
    elif labelling == 'local':
        soln = rf_stacking.find_local_hk_maxima(k_grid, h_grid, hk_stack)
        log.info("Numerical solutions (H, k) = {}".format(soln))
    # end if

    # Plot the local maxima
    for i, (h, k) in enumerate(soln):
        _plot_hk_solution_point(plt.gca(), k, h, i+1)
    # end for

    return fig, soln
# end func

def _produce_sediment_hk_stacking(channel_data, H_c, k_c, labelling=DEFAULT_HK_SOLN_LABEL):
    """Helper function to produce H-k stacking figure."""

    k_grid, h_grid, hk_stack, weighting = rf_stacking.compute_sediment_hk_stack(channel_data,
                                                                     H_c=H_c, k_c=k_c,
                                                                     h_range=rf_stacking.DEFAULT_SED_H_RANGE,
                                                                     k_range=rf_stacking.DEFAULT_SED_k_RANGE,
                                                                     root_order=9)

    sta = channel_data[0].stats.station
    channel = channel_data[0].stats.channel
    num = len(channel_data)
    title = sta + '.{} (Sediment thickness)'.format(channel)

    fig = rf_plot_utils.plot_hk_stack(k_grid, h_grid, hk_stack, title=title, num=num,
                                      depth_colour_range=(np.min(h_grid), np.max(h_grid)),
                                      stack_ylabel='Sediment thickness')

    # Stamp weightings onto plot
    xl = plt.xlim()
    yl = plt.ylim()
    txt_x = xl[0] + 0.95*(xl[1] - xl[0])
    txt_y = yl[0] + 0.90*(yl[1] - yl[0])
    plt.text(txt_x, txt_y, "w={}".format(weighting), horizontalalignment='right', color="#ffffff",
             fontsize=12, rasterized=True)

    # Find and label location of maximum
    soln = []
    if labelling == 'global':
        h_max, k_max = rf_stacking.find_global_hk_maximum(k_grid, h_grid, hk_stack)
        soln = [(h_max, k_max)]
        log.info("Numerical solution (H, k) = ({:.3f}, {:.3f})".format(*soln[0]))
    elif labelling == 'local':
        soln = rf_stacking.find_local_hk_maxima(k_grid, h_grid, hk_stack)
        log.info("Numerical solutions (H, k) = {}".format(soln))
    # end if

    # Plot the local maxima
    for i, (h, k) in enumerate(soln):
        _plot_hk_solution_point(plt.gca(), k, h, i+1)
    # end for

    return fig, soln
# end func

def _plot_hk_solution_point(axes, k, h, idx):
    xl = axes.get_xlim()
    yl = axes.get_ylim()

    axes.text(k, h, "  %d"%(idx),
              color="#000000", fontsize=10, horizontalalignment='right',
              clip_on=False, rasterized=True)

    axes.scatter(k, h, marker='+', c="#000000", s=20)

    x = np.mean(np.array(xl))
    axes.text(x, h, "H{}={:.3f}, k{}={:.3f}".format(idx, h, idx, k),
              color="#000000", fontsize=12, horizontalalignment='left',
              clip_on=False)
# end func

@click.command()
@click.argument('input-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output-file', type=click.Path(dir_okay=False), required=True)
@click.option('--network-list', default='*', help='A space-separated list of networks (within quotes) to process.', type=str,
              show_default=True)
@click.option('--station-list', default='*', help='A space-separated list of stations (within quotes) to process.', type=str,
              show_default=True)
@click.option('--event-mask-folder', type=click.Path(dir_okay=True, exists=True, file_okay=False),
              help='Folder containing event masks to use to filter traces. Such masks are generated '
                   'using rf_handpick_tool')
@click.option('--apply-amplitude-filter', is_flag=True, default=False, show_default=True,
              help='Apply RF amplitude filtering to the RFs. The default filtering logic includes: '
                   'Signal SNR >= 2.0 '
                   'RMS amplitude of signal < 0.2 '
                   'Maximum amplitude of signal < 1.0')
@click.option('--apply-similarity-filter', is_flag=True, default=False, show_default=True,
              help='Apply RF similarity filtering to the RFs.')
@click.option('--min-slope-ratio', type=float, default=-1, show_default=True,
              help='Apply filtering to the RFs based on the "slope_ratio" metric that indicates robustness'
                   'of P-arrival. Typically, a minimum slope-ratio of 5 is able to pick out strong arrivals. '
                   'The default value of -1 does not apply this filter')
@click.option('--hk-weights', type=(float, float, float), default=(0.5, 0.4, 0.1), show_default=True,
              help='Weightings per arrival multiple for H-k stacking')
@click.option('--hk-solution-labels', type=click.Choice(['global', 'local', 'none']), default=DEFAULT_HK_SOLN_LABEL,
              show_default=True, help='Method of labeling automatically selected solutions on H-k stack plots. '
              'global: find and label global maximum, local: find and label up to 3 local maxima after '
              'clustering, none: do not label any solutions on H-k stack plot.')
@click.option('--depth-colour-range', type=(float, float), default=(20, 70), show_default=True,
              help='The range of depth values from which to choose the maximum hk_stack value for plotting '
                   'purposes. Note that this parameter has no effect on the computation of the hk_stack.')
@click.option('--hk-hpf-freq', type=float, default=None, show_default=True,
              help='If present, cutoff frequency for high pass filter to use prior to generating H-k stacking plot.')
def main(input_file, output_file, network_list='*', station_list='*', event_mask_folder='',
         apply_amplitude_filter=False, apply_similarity_filter=False, min_slope_ratio=-1,
         hk_weights=rf_stacking.DEFAULT_WEIGHTS, hk_solution_labels=DEFAULT_HK_SOLN_LABEL,
         depth_colour_range=(20, 70), hk_hpf_freq=None):
    # docstring redundant since CLI options are already documented.

    log.setLevel(logging.INFO)

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()
    proc_hdfkeys = None

    if(rank == 0):
        # retrieve all available hdf_keys
        proc_hdfkeys = rf_util.get_hdf_keys(input_file)

        # trim stations to be processed based on the user-provided network- and station-list
        proc_hdfkeys = rf_util.trim_hdf_keys(proc_hdfkeys, network_list, station_list)

        # split work-load over all procs
        proc_hdfkeys = rf_util.split_list(proc_hdfkeys, nproc)
    # end if

    # broadcast workload to all procs
    proc_hdfkeys = comm.bcast(proc_hdfkeys, root=0)

    pbar = tqdm.tqdm(total=len(proc_hdfkeys[rank]))

    pdf_names = []
    hk_soln = dict()
    sediment_hk_soln = dict()
    station_coords = dict()
    sediment_station_coords = dict()

    for proc_hdfkey in proc_hdfkeys[rank]:
        data_all = rf.read_rf(input_file, format='h5', group='waveforms/%s'%(proc_hdfkey))
        # Convert to hierarchical dictionary format
        data_dict = rf_util.rf_to_dict(data_all)

        nsl = proc_hdfkey # network-station-location
        pbar.update()
        pbar.set_description("Rank {}: {}".format(rank, nsl))

        event_mask_dict = None
        if event_mask_folder and os.path.isdir(event_mask_folder):
            log.info("Applying event mask from folder {}".format(event_mask_folder))
            mask_files = os.listdir(event_mask_folder)
            mask_files = [f for f in mask_files if os.path.isfile(os.path.join(event_mask_folder, f))]
            pattern = r"([A-Za-z0-9\.]{5,})_event_mask\.txt"
            pattern = re.compile(pattern)
            event_mask_dict = dict()
            for f in mask_files:
                match_result = pattern.match(f)
                if not match_result:
                    continue
                code = match_result[1]
                with open(os.path.join(event_mask_folder, f), 'r') as _f:
                    events = _f.readlines()
                    events = set([e.strip() for e in events])
                    event_mask_dict[code] = events
                # end with
            # end for
        # end if

        if event_mask_dict:
            log.info("Loaded {} event masks".format(len(event_mask_dict)))
        # end if

        # Plot all data to PDF file
        fixed_stack_height_inches = 0.8
        y_pad_inches = 1.6
        total_trace_height_inches = paper_size_A4[1] - fixed_stack_height_inches - y_pad_inches
        max_trace_height = 0.2

        log.setLevel(logging.WARNING)

        curr_output_file, _ = os.path.splitext(output_file)
        curr_output_file += '.{}.pdf'.format(nsl)

        with PdfPages(curr_output_file) as pdf:
            # Would like to use Tex, but lack desktop PC privileges to update packages to what is required
            plt.rc('text', usetex=False)
            network = data_dict.network
            rf_type = data_dict.rotation
            for st in sorted(data_dict.keys()):
                station_db = data_dict[st]

                # Choose RF channel
                channel = rf_util.choose_rf_source_channel(rf_type, station_db)
                channel_data = station_db[channel]
                if not channel_data:
                    continue
                # end if
                full_code = '.'.join([network, st, channel])

                t_channel = list(channel)
                t_channel[-1] = 'T'
                t_channel = ''.join(t_channel)

                rf_stream = rf.RFStream(channel_data).sort(['back_azimuth'])
                if event_mask_dict and full_code in event_mask_dict:
                    # Select events from external source
                    event_mask = event_mask_dict[full_code]
                    rf_stream = rf.RFStream(
                        [tr for tr in rf_stream if tr.stats.event_id in event_mask]).sort(['back_azimuth'])
                # end if
                if apply_amplitude_filter:
                    # Label and filter quality
                    rf_util.label_rf_quality_simple_amplitude(rf_type, rf_stream)
                    rf_stream = rf.RFStream(
                        [tr for tr in rf_stream if tr.stats.predicted_quality == 'a']).sort(['back_azimuth'])
                # end if

                if(min_slope_ratio>0):
                    rf_stream = rf.RFStream([tr for tr in rf_stream \
                                             if tr.stats.slope_ratio > min_slope_ratio]).sort(['back_azimuth'])
                # end if

                if apply_similarity_filter and len(rf_stream) >= 3:
                    rf_stream = rf_util.filter_crosscorr_coeff(rf_stream)
                # end if

                if not rf_stream:
                    continue

                # Find matching T-component data
                events = [tr.stats.event_id for tr in rf_stream]
                transverse_data = station_db[t_channel]
                t_stream = rf.RFStream(
                    [tr for tr in transverse_data if tr.stats.event_id in events]).sort(['back_azimuth'])

                # Plot pinwheel of primary and transverse components
                fig = rf_plot_utils.plot_rf_wheel([rf_stream, t_stream], fontscaling=0.8)
                fig.set_size_inches(*paper_size_A4)
                plt.tight_layout()
                plt.subplots_adjust(hspace=0.15, top=0.95, bottom=0.15)
                ax = fig.gca()
                fig.text(-0.32, -0.32, "\n".join(rf_stream[0].stats.processing), fontsize=6,
                         transform=ax.transAxes, rasterized=True)
                pdf.savefig(dpi=300, orientation='portrait')
                plt.close()

                num_traces = len(rf_stream)
                assert len(t_stream) == num_traces or not t_stream

                # Filter rf_stream if needed
                if(hk_hpf_freq and hk_hpf_freq>0):
                    rf_stream.filter(type='highpass', freq=hk_hpf_freq,
                                     corners=1, zerophase=True)
                # end if

                # Plot RF stack of primary component
                trace_ht = min(total_trace_height_inches/num_traces, max_trace_height)
                fig = rf_plot_utils.plot_rf_stack(rf_stream, trace_height=trace_ht, stack_height=fixed_stack_height_inches,
                                                  fig_width=paper_size_A4[0])
                fig.suptitle("Channel {}".format(rf_stream[0].stats.channel))
                # Customize layout to pack to top of page while preserving RF plots aspect ratios
                _rf_layout_A4(fig)
                # Save to new page in file
                pdf.savefig(dpi=300, orientation='portrait')
                plt.close()

                reverberations_removed = False
                # Apply reverberation filter if needed
                if(rf_corrections.has_reverberations(rf_stream)):
                    rf_stream = rf_corrections.apply_reverberation_filter(rf_stream)
                    reverberations_removed = True
                # end if

                if(reverberations_removed):
                    # Plot reverberation-filtered RF stack
                    trace_ht = min(total_trace_height_inches/num_traces, max_trace_height)
                    fig = rf_plot_utils.plot_rf_stack(rf_stream, trace_height=trace_ht, stack_height=fixed_stack_height_inches,
                                                      fig_width=paper_size_A4[0])
                    fig.suptitle("Channel {} (Reverberations removed)".format(rf_stream[0].stats.channel))
                    # Customize layout to pack to top of page while preserving RF plots aspect ratios
                    _rf_layout_A4(fig)
                    # Save to new page in file
                    pdf.savefig(dpi=300, orientation='portrait')
                    plt.close()
                # end if

                # Plot RF stack of transverse component
                if t_stream:
                    fig = rf_plot_utils.plot_rf_stack(t_stream, trace_height=trace_ht,
                                                      stack_height=fixed_stack_height_inches,
                                                      fig_width=paper_size_A4[0])
                    fig.suptitle("Channel {}".format(t_stream[0].stats.channel))
                    # Customize layout to pack to top of page while preserving RF plots aspect ratios
                    _rf_layout_A4(fig)
                    # Save to new page in file
                    pdf.savefig(dpi=300, orientation='portrait')
                    plt.close()
                # end if

                # Plot H-k stack using primary RF component
                fig, maxima = _produce_hk_stacking(rf_stream, weighting=hk_weights,
                                                   labelling=hk_solution_labels,
                                                   depth_colour_range=depth_colour_range)
                hk_soln[nsl] = maxima
                station_coords[nsl] = (channel_data[0].stats.station_latitude, channel_data[0].stats.station_longitude)

                paper_landscape = (paper_size_A4[1], paper_size_A4[0])
                fig.set_size_inches(*paper_landscape)
                # plt.tight_layout()
                # plt.subplots_adjust(hspace=0.15, top=0.95, bottom=0.15)
                pdf.savefig(dpi=300, orientation='landscape')
                plt.close()

                if reverberations_removed:
                    H_c = hk_soln[nsl][0][0]
                    k_c = hk_soln[nsl][0][1]
                    fig, maxima = _produce_sediment_hk_stacking(rf_stream, H_c=H_c, k_c=k_c)
                    sediment_hk_soln[nsl] = maxima

                    sediment_station_coords[nsl] = (channel_data[0].stats.station_latitude,
                                                    channel_data[0].stats.station_longitude)
                    fig.set_size_inches(*paper_landscape)
                    pdf.savefig(dpi=300, orientation='landscape')
                    plt.close()
                # end if

                plt.close('all')
            # end for
        # end with
        pdf_names.append(curr_output_file)
    # end for
    pbar.close()

    # gather hk_soln, sediment_hk_soln and associated coordinates on rank 0
    hk_soln = comm.gather(hk_soln, root=0)
    sediment_hk_soln = comm.gather(sediment_hk_soln, root=0)
    station_coords = comm.gather(station_coords, root=0)
    sediment_station_coords = comm.gather(sediment_station_coords, root=0)
    pdf_names = comm.gather(pdf_names, root=0)

    comm.Barrier()
    if(rank==0):
        def flatten_dict_list(dict_list):
            return {k: v for d in dict_list for k, v in d.items()}
        #end func

        hk_soln = flatten_dict_list(hk_soln)
        sediment_hk_soln = flatten_dict_list(sediment_hk_soln)
        station_coords = flatten_dict_list(station_coords)
        sediment_station_coords = flatten_dict_list(sediment_station_coords)
        # write solutions to csv files
        for result_hk, result_coords, fname in zip([hk_soln, sediment_hk_soln],
                                                     [station_coords, sediment_station_coords],
                                                     [os.path.splitext(output_file)[0] + '.csv',
                                                      os.path.splitext(output_file)[0] + '.sed.csv']):
            assert len(result_hk) == len(result_coords)

            if(len(result_hk) == 0): continue

            # Sort H-k solutions by depth from low to high
            update_dict = {}
            for nsl, hks in result_hk.items():
                sorted_hks = sorted([tuple(hk) for hk in hks])
                update_dict[nsl] = np.array(list(result_coords[nsl]) + [i for hk in sorted_hks for i in hk])
            # end for
            result_hk.update(update_dict)

            df = pd.DataFrame.from_dict(result_hk, orient='index')
            colnames = [('H{}'.format(i), 'k{}'.format(i)) for i in range((len(df.columns) - 2)//2)]
            colnames = ['Latitude', 'Longitude'] + list(itertools.chain.from_iterable(colnames))
            df.columns = colnames
            df.index.name = 'Station'
            df.to_csv(fname)
        # end for

        # gather pdf-names, flatten list and merge pdfs
        pdf_names = [item for items in pdf_names for item in items]
        pdf_merge(pdf_names, output_file)
    # end if
# end main

if __name__ == "__main__":
    log = logging.getLogger(__name__)
    main()  # pylint: disable=no-value-for-parameter
# end if
