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
from seismic.stream_io import get_obspyh5_index
import uuid
from shutil import rmtree
from collections import defaultdict

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


def _produce_hk_stacking(channel_data, Vp=rf_stacking.DEFAULT_Vp,
                         weighting=rf_stacking.DEFAULT_WEIGHTS,
                         semblance_weighted=True, labelling=DEFAULT_HK_SOLN_LABEL,
                         depth_colour_range=(20, 70)):
    """Helper function to produce H-k stacking figure."""

    k_grid, h_grid, hk_stack = rf_stacking.compute_hk_stack(channel_data,
                                                            Vp=Vp,
                                                            h_range=rf_stacking.DEFAULT_H_RANGE,
                                                            k_range=rf_stacking.DEFAULT_k_RANGE,
                                                            weights=weighting,
                                                            root_order=2,
                                                            semblance_weighted=semblance_weighted)

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

    # ==========================================================================
    # Plot the local maxima
    # ==========================================================================
    def compute_arrival_times(p, h, k):
        Vp_inv = 1. / Vp

        Vs_inv = k * Vp_inv
        term1 = np.sqrt(Vs_inv ** 2 - p ** 2)
        term2 = np.sqrt(Vp_inv ** 2 - p ** 2)

        t1 = h * (term1 - term2)
        t2 = h * (term1 + term2)
        t3 = h * 2 * term1

        return t1, t2, t3
    # end func
    for i, (h, k) in enumerate(soln):
        phase_time_dict = {}
        p_array = np.array([trc.stats.slowness / rf_stacking.DEG2KM for trc in channel_data])
        p_mean = np.mean(p_array)
        p_std = np.std(p_array)
        min_times = compute_arrival_times(p_mean - p_std, h, k)
        max_times = compute_arrival_times(p_mean + p_std, h, k)
        phase_time_dict['t1'] = sorted(np.array([min_times[0], max_times[0]]))
        phase_time_dict['t2'] = sorted(np.array([min_times[1], max_times[1]]))
        phase_time_dict['t3'] = sorted(np.array([min_times[2], max_times[2]]))
        _plot_hk_solution_point(plt.gca(), k, h, i, phase_time_dict)
    # end for

    return fig, soln
# end func

def _produce_sediment_hk_stacking(channel_data, H_c, k_c, Vp=rf_stacking.DEFAULT_Vp,
                                  labelling=DEFAULT_HK_SOLN_LABEL):
    """Helper function to produce H-k stacking figure."""

    k_grid, h_grid, hk_stack, weighting = rf_stacking.compute_sediment_hk_stack(channel_data,
                                                                     H_c=H_c, k_c=k_c, Vp=Vp,
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
        _plot_hk_solution_point(plt.gca(), k, h, i)
    # end for

    return fig, soln
# end func

def _plot_hk_solution_point(axes, k, h, idx, phase_time_dict=None):
    xl = axes.get_xlim()
    yl = axes.get_ylim()

    axes.text(k, h, "  %d"%(idx),
              color="#000000", fontsize=10, horizontalalignment='right',
              clip_on=False, rasterized=True)

    axes.scatter(k, h, marker='+', c="#000000", s=20)

    x = np.mean(np.array(xl))
    x1 = x + (xl[1]-x)*0.55
    axes.text(x, h, "H{}={:.3f}, k{}={:.3f}".format(idx, h, idx, k),
              color="#000000", fontsize=12, horizontalalignment='left',
              clip_on=False)
    if(phase_time_dict):
        axes.text(x1, h+0.6, "Ps: [{:.1f} - {:.1f}] s".format(*phase_time_dict['t1']),
                  color="#000000", fontsize=4, horizontalalignment='left',
                  clip_on=False)
        axes.text(x1, h, "PpPs: [{:.1f} - {:.1f}] s".format(*phase_time_dict['t2']),
                  color="#000000", fontsize=4, horizontalalignment='left',
                  clip_on=False)
        axes.text(x1, h-0.6, "PpSs+PsPs: [{:.1f} - {:.1f}] s".format(*phase_time_dict['t3']),
                  color="#000000", fontsize=4, horizontalalignment='left',
                  clip_on=False)
    # end if
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
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
              help='Discard RFs that do not meet the specified minimum "slope_ratio" metric that indicates '
                   'robustness of P-arrival. Typically, a minimum slope-ratio of 5 is able to pick out '
                   'strong arrivals. The default value of -1 does not apply this filter')
@click.option('--force-dereverberate', default=None, show_default=False,
              help='A space-separated list of stations (within quotes) for which the de-reverberation filter '
                   'is to be triggered by force.')
@click.option('--filter-by-distance', default=None, type=(str, float, float), show_default=False, multiple=True,
              help='RFs for particular stations can be restricted to a given range of angular distance. Parameters '
                   'are specified as a space-separated triplet: station-name min-dist max-dist. Note that this '
                   'parameter can be repeated for multiple stations')
@click.option('--hk-vp', type=float, default=rf_stacking.DEFAULT_Vp, show_default=True,
              help='Value to use for Vp (crustal P-wave velocity [km/s]) for H-k stacking')
@click.option('--hk-weights', type=(float, float, float), default=(0.5, 0.4, 0.1), show_default=True,
              help='Weightings per arrival multiple for H-k stacking')
@click.option('--hk-solution-labels', type=click.Choice(['global', 'local', 'none']), default=DEFAULT_HK_SOLN_LABEL,
              show_default=True, help='Method of labeling automatically selected solutions on H-k stack plots. '
              'global: find and label global maximum, local: find and label up to 3 local maxima after '
              'clustering, none: do not label any solutions on H-k stack plot')
@click.option('--depth-colour-range', type=(float, float), default=(20, 70), show_default=True,
              help='The range of depth values from which to choose the maximum hk_stack value for plotting '
                   'purposes. Note that this parameter has no effect on the computation of the hk_stack')
@click.option('--hk-hpf-freq', type=float, default=None, show_default=True,
              help='If present, cutoff frequency for high pass filter to use prior to generating H-k stacking plot')
@click.option('--disable-semblance-weighting', is_flag=True, default=False, show_default=True,
              help='Disables default semblance-weighting applied to H-k stacks')
def main(input_file, output_file, network_list='*', station_list='*', event_mask_folder='',
         apply_amplitude_filter=False, apply_similarity_filter=False, min_slope_ratio=-1,
         force_dereverberate=None, filter_by_distance=None,
         hk_vp=rf_stacking.DEFAULT_Vp, hk_weights=rf_stacking.DEFAULT_WEIGHTS,
         hk_solution_labels=DEFAULT_HK_SOLN_LABEL, depth_colour_range=(20, 70),
         hk_hpf_freq=None, disable_semblance_weighting=False):

    """
    INPUT_FILE : Input RFs in H5 format\n
                 (output of generate_rf.py or rf_quality_filter.py)\n
    OUTPUT_FILE : Output pdf file name
    """


    log.setLevel(logging.INFO)

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()
    proc_hdfkeys = None

    tempdir = None
    fd_sta_list = []
    filter_by_distance_dict = defaultdict(list)
    if(rank == 0):
        # retrieve all available hdf_keys
        proc_hdfkeys = get_obspyh5_index(input_file, seeds_only=True)

        ###############################################################################
        # Parse special parameters e.g. for invoking forced de-reverberation or for
        # limiting RFs by distance for particular stations.
        ###############################################################################
        def match_stations(sta, sta_list):
            matches_found = []
            for tsta in sta_list:
                if sta in tsta:
                    matches_found.append(tsta)
                # end if
            # end for
            return matches_found
        # end func
        # parse --force-dereverberate
        if(force_dereverberate is not None):
            fd_sta_list = re.findall('\S+', force_dereverberate)
            assert len(fd_sta_list), 'Invalid station list for triggering forced dereverberation. Aborting..'
            for fd_sta in fd_sta_list:
                matches_found = match_stations(fd_sta, proc_hdfkeys)
                if(len(matches_found) == 0): assert 0, 'Station {} not found. Aborting..'.format(fd_sta)
            # end for
        # end if
        # parse --filter-by-distance
        if(filter_by_distance is not None):
            for row in filter_by_distance:
                fd_sta, min_dist, max_dist = row
                matches_found = match_stations(fd_sta, proc_hdfkeys)
                if(len(matches_found) == 0): assert 0, 'Station {} not found. Aborting..'.format(fd_sta)
                if(min_dist > max_dist): assert 0, 'Invalid parameters for station {}: min_dist > max_dist ' \
                                                   'Aborting..'.format(fd_sta)
                filter_by_distance_dict[fd_sta] = [min_dist, max_dist]
            # end for
        # end if

        # trim stations to be processed based on the user-provided network- and station-list
        proc_hdfkeys = rf_util.trim_hdf_keys(proc_hdfkeys, network_list, station_list)

        # split work-load over all procs
        proc_hdfkeys = rf_util.split_list(proc_hdfkeys, nproc)

        tempdir = os.path.join(os.path.dirname(output_file), str(uuid.uuid4()))
        os.makedirs(tempdir, exist_ok=True)
    # end if
    tempdir = comm.bcast(tempdir, root=0)

    # broadcast workload to all procs
    proc_hdfkeys = comm.bcast(proc_hdfkeys, root=0)
    fd_sta_list = comm.bcast(fd_sta_list, root=0)
    filter_by_distance_dict = comm.bcast(filter_by_distance_dict, root=0)

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

        curr_output_file = os.path.join(tempdir, '{}.pdf'.format(nsl))

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
                transverse_data = station_db[t_channel]

                rf_stream = rf.RFStream(channel_data).sort(['back_azimuth'])
                if event_mask_dict and full_code in event_mask_dict:
                    # Select events from external source
                    event_mask = event_mask_dict[full_code]
                    rf_stream = rf.RFStream(
                        [tr for tr in rf_stream if tr.stats.event_id in event_mask]).sort(['back_azimuth'])
                # end if

                ###############################################################################
                # Restrict RFs by distance if specified
                ###############################################################################
                if st in filter_by_distance_dict.keys():
                    min_dist, max_dist = filter_by_distance_dict[st]

                    rf_stream = rf_util.filter_by_distance(rf_stream, min_dist, max_dist)
                    if (len(rf_stream) == 0):
                        log.warning("Filter-by-distance has removed all traces for {}.".format(nsl))
                    # end if
                # end if

                ###############################################################################
                # Apply amplitude filter
                ###############################################################################
                if apply_amplitude_filter:
                    # Label and filter quality
                    rf_util.label_rf_quality_simple_amplitude(rf_type, rf_stream)
                    rf_stream = rf.RFStream(
                        [tr for tr in rf_stream if tr.stats.predicted_quality == 'a']).sort(['back_azimuth'])
                    if(len(rf_stream) == 0):
                        log.warning("Amplitude filter has removed all traces for {}. "
                                    "Ensure rf_quality_filter was run beforehand..".format(nsl))
                    # end if
                # end if

                ###############################################################################
                # Apply slope-ratio filter
                ###############################################################################
                if(min_slope_ratio>0):
                    rf_stream = rf.RFStream([tr for tr in rf_stream \
                                             if tr.stats.slope_ratio > min_slope_ratio]).sort(['back_azimuth'])
                # end if

                ###############################################################################
                # Apply similarity filter
                ###############################################################################
                if apply_similarity_filter and len(rf_stream) >= 3:
                    rf_stream = rf_util.filter_crosscorr_coeff(rf_stream)
                # end if

                if not rf_stream:
                    continue

                ###############################################################################
                # Filter rf_stream if needed
                ###############################################################################
                if(hk_hpf_freq and hk_hpf_freq>0):
                    rf_stream.filter(type='highpass', freq=hk_hpf_freq,
                                     corners=1, zerophase=True)
                # end if

                rf_stream_raw = rf_stream.copy()
                ###############################################################################
                # Apply reverberation filter if needed
                ###############################################################################
                reverberations_removed = False
                # Apply reverberation filter if needed (or forced)
                if(rf_corrections.has_reverberations(rf_stream) or \
                   st in fd_sta_list):
                    rf_stream = rf_corrections.apply_reverberation_filter(rf_stream)
                    reverberations_removed = True
                # end if

                ###############################################################################
                # RF amplitudes should not exceed 1.0 and should peak around onset time --
                # otherwise, such traces are deemed problematic and excluded while computing
                # H-k stacks.
                ###############################################################################
                before = len(rf_stream)
                rf_stream = rf_util.filter_invalid_radial_component(rf_stream)
                after = len(rf_stream)

                if(before > after):
                    print("{}: {}/{} RF traces with amplitudes > 1.0 or troughs around onset time dropped "
                          "before computing H-k stack ..".format(full_code,
                                                                 before - after,
                                                                 before))
                # end if

                if not rf_stream:
                    continue

                ############################################
                # Collate streams to plot
                ############################################
                events = [tr.stats.event_id for tr in rf_stream]
                num_traces = len(rf_stream)

                rf_stream_raw = rf.RFStream(
                    [tr for tr in rf_stream_raw if tr.stats.event_id in events]).sort(['back_azimuth'])
                t_stream = rf.RFStream(
                    [tr for tr in transverse_data if tr.stats.event_id in events]).sort(['back_azimuth'])
                assert len(t_stream) == num_traces or not t_stream

                ############################################
                # Plot psd
                ############################################
                fig, ax = plt.subplots()
                fig.set_size_inches(paper_size_A4[1], paper_size_A4[0])
                fig.suptitle('.'.join([nsl, channel]))
                ax.set_rasterized(True)

                rf_plot_utils.plot_rf_psd(rf_stream, ax, min_slope_ratio=min_slope_ratio)
                pdf.savefig(dpi=300, orientation='landscape')
                plt.close()

                ############################################
                # Plot pinwheel of primary and transverse components
                ############################################
                fig = rf_plot_utils.plot_rf_wheel([rf_stream_raw, t_stream], fontscaling=0.8)
                fig.set_size_inches(*paper_size_A4)
                plt.tight_layout()
                plt.subplots_adjust(hspace=0.15, top=0.95, bottom=0.15)
                ax = fig.gca()
                fig.text(-0.32, -0.32, "\n".join(rf_stream[0].stats.processing), fontsize=6,
                         transform=ax.transAxes, rasterized=True)
                pdf.savefig(dpi=300, orientation='portrait')
                plt.close()

                ############################################
                # Plot RF waveform stacks
                ############################################
                trace_ht = min(total_trace_height_inches/num_traces, max_trace_height)
                # Plot raw RF stack of primary component
                fig = rf_plot_utils.plot_rf_stack(rf_stream_raw, trace_height=trace_ht,
                                                  stack_height=fixed_stack_height_inches,
                                                  fig_width=paper_size_A4[0])
                fig.suptitle("Channel {}".format(rf_stream[0].stats.channel))
                # Customize layout to pack to top of page while preserving RF plots aspect ratios
                _rf_layout_A4(fig)
                # Save to new page in file
                pdf.savefig(dpi=300, orientation='portrait')
                plt.close()

                if(reverberations_removed):
                    # Plot reverberation-filtered RF stack
                    fig = rf_plot_utils.plot_rf_stack(rf_stream, trace_height=trace_ht,
                                                      stack_height=fixed_stack_height_inches,
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

                ############################################
                # Plot H-k stack using primary RF component
                ############################################
                fig, maxima = _produce_hk_stacking(rf_stream, Vp=hk_vp, weighting=hk_weights,
                                                   semblance_weighted=(not disable_semblance_weighting),
                                                   labelling=hk_solution_labels,
                                                   depth_colour_range=depth_colour_range)
                hk_soln[nsl] = maxima
                station_coords[nsl] = (channel_data[0].stats.station_longitude, channel_data[0].stats.station_latitude)

                paper_landscape = (paper_size_A4[1], paper_size_A4[0])
                fig.set_size_inches(*paper_landscape)
                # plt.tight_layout()
                # plt.subplots_adjust(hspace=0.15, top=0.95, bottom=0.15)
                pdf.savefig(dpi=300, orientation='landscape')
                plt.close()

                if reverberations_removed:
                    if(len(hk_soln[nsl])):
                        H_c = hk_soln[nsl][0][0]
                        k_c = hk_soln[nsl][0][1]
                        fig, maxima = _produce_sediment_hk_stacking(rf_stream, H_c=H_c, k_c=k_c, Vp=hk_vp)
                        sediment_hk_soln[nsl] = maxima

                        sediment_station_coords[nsl] = (channel_data[0].stats.station_longitude,
                                                        channel_data[0].stats.station_latitude)
                        fig.set_size_inches(*paper_landscape)
                        pdf.savefig(dpi=300, orientation='landscape')
                        plt.close()
                    else:
                        log.warning("Sediment H-K stacking for {} failed. Moving along..".format(nsl))
                    # end if
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
            colnames = ['Longitude', 'Latitude'] + list(itertools.chain.from_iterable(colnames))
            df.columns = colnames
            df.index.name = 'Station'
            df.to_csv(fname)
        # end for

        # gather pdf-names, flatten list and merge pdfs
        pdf_names = [item for items in pdf_names for item in items]
        pdf_merge(pdf_names, output_file)
    # end if

    if(rank == 0):
        rmtree(tempdir)

        print("Finishing...")
        print("bulk_rf_report SUCCESS!")
    # end if
# end main

if __name__ == "__main__":
    log = logging.getLogger(__name__)
    main()  # pylint: disable=no-value-for-parameter
# end if
