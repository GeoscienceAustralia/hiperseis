#!/usr/bin/env python
# coding: utf-8
"""Produce PDF report of network stations showing RF waveforms
"""

import os
import re
import logging

import numpy as np
import click
import rf
import rf.imaging

import tqdm.auto as tqdm

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

import seismic.receiver_fn.rf_util as rf_util
import seismic.receiver_fn.rf_plot_utils as rf_plot_utils
import seismic.receiver_fn.rf_stacking as rf_stacking

# pylint: disable=invalid-name, logging-format-interpolation

logging.basicConfig()

paper_size_A4 = (8.27, 11.69)  # inches

DEFAULT_HK_WEIGHTS = (0.5, 0.5, 0.0)


def _get_aspect(ax):
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


def _produce_hk_stacking(channel_data, V_p=6.4, weighting=(0.5, 0.5, 0.0)):

    k_grid, h_grid, hk_stack = rf_stacking.compute_hk_stack(channel_data,
                                                            h_range=np.linspace(20.0, 70.0, 501),
                                                            root_order=2, V_p=V_p)

    layer_maxes = [np.max(hk_stack[i, :, :]) for i in range(3)]
    log.info("Max amplitude by layer: {}".format(layer_maxes))

    # Sum the phases
    hk_stack_sum = rf_stacking.compute_weighted_stack(hk_stack, weighting)

    # Raise the final sum over phases to power >1 to increase contrast
    hk_stack_sum = rf_util.signed_nth_power(hk_stack_sum, 2)
    hk_stack_sum = hk_stack_sum/np.nanmax(hk_stack_sum[:])

    sta = channel_data[0].stats.station
    channel = channel_data[0].stats.channel
    num = len(channel_data)
    fig = rf_plot_utils.plot_hk_stack(k_grid, h_grid, hk_stack_sum, title=sta + '.{}'.format(channel), num=num)

    # Stamp weightings onto plot
    xl = plt.xlim()
    yl = plt.ylim()
    txt_x = xl[0] + 0.95*(xl[1] - xl[0])
    txt_y = yl[0] + 0.90*(yl[1] - yl[0])
    plt.text(txt_x, txt_y, "w={}".format(weighting), horizontalalignment='right', color="#ffffff", fontsize=12)

    # Find and label location of maximum
    h_max, k_max = rf_stacking.find_global_hk_maximum(k_grid, h_grid, hk_stack_sum)
    log.info("Numerical solution (H, k) = ({:.3f}, {:.3f})".format(h_max, k_max))
    plt.scatter(k_max, h_max, marker='+', c="#000000", s=20)
    if k_max >= 1.7:
        plt.text(k_max - 0.01, h_max + 1, "Solution H = {:.3f}, k = {:.3f}".format(h_max, k_max),
                 color="#ffffff", fontsize=16, horizontalalignment='right')
    else:
        plt.text(k_max + 0.01, h_max + 1, "Solution H = {:.3f}, k = {:.3f}".format(h_max, k_max),
                 color="#ffffff", fontsize=16)
    # end if

    return fig
# end func


@click.command()
@click.argument('input-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output-file', type=click.Path(dir_okay=False), required=True)
@click.option('--event-mask-folder', type=click.Path(dir_okay=True, exists=True, file_okay=False),
              help='Folder containing event masks to use to filter traces. Such masks are generated '
                   'using rf_handpick_tool')
@click.option('--apply-amplitude-filter', is_flag=True, default=False, show_default=True,
              help='Apply RF amplitude filtering to the RFs.')
@click.option('--apply-similarity-filter', is_flag=True, default=False, show_default=True,
              help='Apply RF similarity filtering to the RFs.')
@click.option('--hk-weights', type=(float, float, float), default=DEFAULT_HK_WEIGHTS, show_default=True,
              help='Weightings per arrival multiple for H-k stacking')
def main(input_file, output_file, event_mask_folder='', apply_amplitude_filter=False,
         apply_similarity_filter=False, hk_weights=DEFAULT_HK_WEIGHTS):

    # Read source file
    data_all = rf_util.read_h5_rf(input_file)

    # Convert to hierarchical dictionary format
    data_dict = rf_util.rf_to_dict(data_all)

    event_mask_dict = None
    if event_mask_folder and os.path.isdir(event_mask_folder):
        mask_files = os.listdir(event_mask_folder)
        mask_files = [f for f in mask_files if os.path.isfile(os.path.join(event_mask_folder, f))]
        # print(mask_files)
        pattern = r"([A-Za-z0-9\.]{5,})_event_mask\.txt"
        pattern = re.compile(pattern)
        event_mask_dict = dict()
        for f in mask_files:
            match_result = pattern.match(f)
            if not match_result:
                continue
            code = match_result[1]
            # print(code)
            with open(os.path.join(event_mask_folder, f), 'r') as f:
                events = f.readlines()
                events = set([e.strip() for e in events])
                event_mask_dict[code] = events
            # end with
        # end for
    # end if

    if event_mask_dict:
        print("Loaded {} event masks".format(len(event_mask_dict)))
    # end if


    # Plot all data to PDF file
    fixed_stack_height_inches = 0.8
    y_pad_inches = 1.6
    total_trace_height_inches = paper_size_A4[1] - fixed_stack_height_inches - y_pad_inches
    max_trace_height = 0.2

    with PdfPages(output_file) as pdf:
        # Would like to use Tex, but lack desktop PC privileges to update packages to what is required
        plt.rc('text', usetex=False)
        pbar = tqdm.tqdm(total=len(data_dict))
        network = data_dict.network
        rf_type = data_dict.rotation
        for st, station_db in data_dict:
            pbar.update()
            pbar.set_description("{}.{}".format(network, st))

            # Choose RF channel
            channel = rf_util.choose_rf_source_channel(rf_type, station_db)
            channel_data = station_db[channel]
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
            if apply_similarity_filter:
                rf_stream = rf_util.filter_crosscorr_coeff(rf_stream)
            # end if

            if not rf_stream:
                continue

            # Find matching T-component data
            events = [tr.stats.event_id for tr in rf_stream]
            transverse_data = station_db[t_channel]
            t_stream = rf.RFStream(
                [tr for tr in transverse_data if tr.stats.event_id in events]).sort(['back_azimuth'])
            if not t_stream:
                continue

            # Plot pinwheel of primary and transverse components
            fig = rf_plot_utils.plot_rf_wheel([rf_stream, t_stream], fontscaling=0.8)
            fig.set_size_inches(*paper_size_A4)
            plt.tight_layout()
            plt.subplots_adjust(hspace=0.15, top=0.95, bottom=0.15)
            ax = fig.gca()
            fig.text(-0.32, -0.32, "\n".join(rf_stream[0].stats.processing), fontsize=6, transform=ax.transAxes)
            pdf.savefig(dpi=300, papertype='a4', orientation='portrait')
            plt.close()

            num_traces = len(rf_stream)
            assert len(t_stream) == num_traces

            # Plot RF stack of primary component
            trace_ht = min(total_trace_height_inches/num_traces, max_trace_height)
            fig = rf_plot_utils.plot_rf_stack(rf_stream, trace_height=trace_ht, stack_height=fixed_stack_height_inches,
                                              fig_width=paper_size_A4[0])
            fig.suptitle("Channel {}".format(rf_stream[0].stats.channel))
            # Customize layout to pack to top of page while preserving RF plots aspect ratios
            _rf_layout_A4(fig)
            # Save to new page in file
            pdf.savefig(dpi=300, papertype='a4', orientation='portrait')
            plt.close()

            # Plot RF stack of transverse component
            fig = rf_plot_utils.plot_rf_stack(t_stream, trace_height=trace_ht, stack_height=fixed_stack_height_inches,
                                              fig_width=paper_size_A4[0])
            fig.suptitle("Channel {}".format(t_stream[0].stats.channel))
            # Customize layout to pack to top of page while preserving RF plots aspect ratios
            _rf_layout_A4(fig)
            # Save to new page in file
            pdf.savefig(dpi=300, papertype='a4', orientation='portrait')
            plt.close()

            # Plot H-k stack using primary RF component
            fig = _produce_hk_stacking(rf_stream, weighting=hk_weights)
            paper_landscape = (paper_size_A4[1], paper_size_A4[0])
            fig.set_size_inches(*paper_landscape)
            # plt.tight_layout()
            # plt.subplots_adjust(hspace=0.15, top=0.95, bottom=0.15)
            pdf.savefig(dpi=300, papertype='a4', orientation='landscape')
            plt.close()

        # end for
        pbar.close()
    # end with

# end main


if __name__ == "__main__":
    log = logging.getLogger(__name__)
    main()  # pylint: disable=no-value-for-parameter
# end if
