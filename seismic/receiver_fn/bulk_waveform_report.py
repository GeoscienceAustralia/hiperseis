#!/usr/bin/env python
# coding: utf-8
"""Produce PDF report of network stations showing waveform characteristics, e.g. power spectral density, etc.
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
from seismic.network_event_dataset import NetworkEventDataset
from scipy import stats
from scipy import signal

logging.basicConfig()

paper_size_A4 = (8.27, 11.69)  # inches

def plot_psd(ned, ax, min_slope_ratio=-1, component='Z'):
    tw = [-10, 25]
    all_trace_lens = []
    traces = []
    for sta, ev, st in ned:
        z = st.select(component=component)[0]

        ot = z.stats.onset
        z = z.copy().slice(ot + tw[0], ot + tw[1])
        if (len(z)):
            all_trace_lens.append(len(z))
            traces.append(z)
        # end if
    # end for

    all_trace_lens = np.array(all_trace_lens)
    most_common_len, _ = stats.mode(all_trace_lens, axis=None)

    # plot psd
    psd_array = []
    slope_ratios = []
    fbins = None
    for length, trace in zip(all_trace_lens, traces):
        if (length == most_common_len):
            fbins, psd = signal.welch(trace.data, fs=trace.stats.sampling_rate,
                                      detrend='linear')

            psd_array.append(psd)
            slope_ratios.append(trace.stats.slope_ratio)

            c = 'm' if trace.stats.slope_ratio >= min_slope_ratio else 'c'
            zo = 10 if trace.stats.slope_ratio >= min_slope_ratio else 1

            ax.loglog(fbins, psd, alpha=0.05, c=c, zorder=zo)
        # end if
        # break
    # end for

    if (len(psd_array)):
        psd_array = np.array(psd_array)
        slope_ratios = np.array(slope_ratios)

        psd_good = np.nanmean(psd_array[slope_ratios >= min_slope_ratio, :], axis=0)
        psd_bad = np.nanmean(psd_array[slope_ratios < min_slope_ratio, :], axis=0)

        ax.loglog(fbins, psd_good, alpha=1, c='m', lw=2, label='Mean PSD (Slope_Ratio >= {})'.format(min_slope_ratio))
        ax.loglog(fbins, psd_bad, alpha=1, c='c', lw=2, label='Mean PSD (Slope_Ratio < {})'.format(min_slope_ratio))
        ax.set_xlabel('Freq. [Hz]')
        ax.set_ylabel('Power Spectral Density [arb. units]')
        ax.text(x=0.9, y=0.85, s='{} Traces'.format(len(psd_array)), transform=ax.transAxes)
        ax.legend()
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
@click.option('--min-slope-ratio', type=float, default=-1, show_default=True,
              help='Filter waveforms based on the "slope_ratio" metric that indicates robustness'
                   'of P-arrival. Typically, a minimum slope-ratio of 5 is able to pick out strong arrivals. '
                   'The default value of -1 does not apply this filter')
def main(input_file, output_file, network_list='*', station_list='*', min_slope_ratio=-1):

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
    if(rank == 0):
        # retrieve all available hdf_keys
        proc_hdfkeys = get_obspyh5_index(input_file, seeds_only=True)

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

    pbar = tqdm.tqdm(total=len(proc_hdfkeys[rank]))

    pdf_names = []

    for proc_hdfkey in proc_hdfkeys[rank]:
        nsl = proc_hdfkey # network-station-location
        pbar.update()
        pbar.set_description("Rank {}: {}".format(rank, nsl))

        curr_output_file = os.path.join(tempdir, '{}.pdf'.format(nsl))

        net, sta, loc = nsl.split('.')
        # note that ned contains a single station
        ned = NetworkEventDataset(input_file, network=net, station=sta, location=loc)

        with PdfPages(curr_output_file) as pdf:
            fig, ax = plt.subplots()
            fig.set_size_inches(paper_size_A4[1], paper_size_A4[0])
            fig.suptitle(nsl)

            ax.set_rasterized(True)

            if(len(ned)): plot_psd(ned, ax, min_slope_ratio=min_slope_ratio)
            pdf.savefig(dpi=300, orientation='landscape')
            plt.close('all')
        # end with
        pdf_names.append(curr_output_file)
    # end for
    pbar.close()

    # gather on rank 0
    pdf_names = comm.gather(pdf_names, root=0)

    comm.Barrier()
    if(rank==0):
        # gather pdf-names, flatten list and merge pdfs
        pdf_names = [item for items in pdf_names for item in items]
        pdf_merge(pdf_names, output_file)

        os.removedirs(tempdir)

        print("Finishing...")
        print("bulk_waveform_report SUCCESS!")
    # end if
# end main

if __name__ == "__main__":
    log = logging.getLogger(__name__)
    main()  # pylint: disable=no-value-for-parameter
# end if
