#!/usr/bin/env python
# coding: utf-8
"""Bulk plotting helper functions based on NetworkEventDataset
"""

import os
import math

import click
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from tqdm.auto import tqdm

from seismic.network_event_dataset import NetworkEventDataset
from seismic.stream_quality_filter import curate_stream3c


def plot_ned_seismograms(ned, output_file, channel_order='ZNE'):
    """
    Plot seismograms in NetworkEventDataset to PDF file.
    If dataset is very large, this may take a long time to run.

    :param ned: NetworkEventDataset containing waveforms.
    :param output_file: Output file name
    """

    if os.path.splitext(output_file)[1].lower() != '.pdf':
        output_file += '.pdf'
    # end if

    # Per page
    ncols = 2
    nrows = 8
    pagesize = (8.3, 8.3*1.414)  # A4
    net = ned.network

    bbstyle = dict(boxstyle="round", fc="w", alpha=0.5, linewidth=0.5)
    annot_fontsize = 5
    axes_fontsize = 6

    # logger = logging.getLogger(__name__)

    out_basename, ext = os.path.splitext(output_file)
    output_file = out_basename + ext

    with PdfPages(output_file) as pdf:
        pb = tqdm(ned.by_station())
        for sta, db_evid in pb:
            seedid = '.'.join([net, sta])
            pb.set_description(seedid)
            pb.write('Processing ' + seedid)
            num_events = len(db_evid)
            num_per_page = nrows*ncols
            ev_ids = list(db_evid.keys())
            num_pages = math.ceil(num_events/num_per_page)
            for pagenum in range(num_pages):
                f = plt.figure(constrained_layout=False, figsize=pagesize)
                f.suptitle('{}.{} event seismograms (pg {}/{})'.format(net, sta, pagenum + 1, num_pages), y=0.98,
                           fontsize=11, va='top')
                pgspec = gridspec.GridSpec(ncols=ncols, nrows=nrows, figure=f, left=0.1, right=0.9, bottom=0.05,
                                           top=0.95, hspace=0.3, wspace=0.2)
                for i, evid in enumerate(ev_ids[pagenum*num_per_page:(pagenum + 1)*num_per_page]):
                    stream = db_evid[evid]
                    # Find max range so we can use a common scale
                    max_half_range = 0.0
                    for tr in stream:
                        max_half_range = max(max_half_range, 0.5*(np.nanmax(tr.data) - np.nanmin(tr.data)))
                    # end for
                    max_half_range *= 1.1

                    # Do plots
                    gs = pgspec[i%nrows, i//nrows].subgridspec(3, 1, hspace=0.0)
                    for tr in stream:
                        cha = tr.stats.channel[-1]
                        idx = channel_order.index(cha)
                        axn = f.add_subplot(gs[idx])
                        t = tr.times() - (tr.stats.onset - tr.stats.starttime)
                        tr_mean = np.nanmean(tr.data)
                        # Rasterize so that size is fairly constant, irrespective
                        # of the number of sample points per trace
                        axn.plot(t, tr.data, 'k', rasterized=True, linewidth=0.6, alpha=0.8)
                        axn.set_ylim(tr_mean - max_half_range, tr_mean + max_half_range)
                        axn.tick_params(labelsize=axes_fontsize)
                        tag = '.'.join([tr.stats.network, tr.stats.station, tr.stats.location, tr.stats.channel])
                        axn.text(0.02, 0.90, tag, fontsize=annot_fontsize, ha='left', va='top', transform=axn.transAxes,
                                bbox=bbstyle)
                        if idx == 0:
                            id_label = 'Event: {}'.format(evid)
                            axn.text(0.98, 0.90, id_label, fontsize=annot_fontsize, ha='right', va='top',
                                    transform=axn.transAxes, bbox=bbstyle)
                        # end if
                        if idx != 2:
                            axn.set_xticks([])
                            axn.xaxis.set_visible(False)
                            axn.xaxis.label.set_visible(False)
                        else:
                            axn.xaxis.label.set_size(axes_fontsize)
                            axn.text(0.5, 0.02, 'Onset: {}'.format(str(tr.stats.onset)), fontsize=annot_fontsize,
                                    ha='center', va='bottom', transform=axn.transAxes)
                        # end if
                        axn.yaxis.label.set_size(axes_fontsize)
                        # f.add_subplot(axn)
                    # end for
                # end for

                pdf.savefig(dpi=300, papertype='a4', orientation='portrait')
                plt.close()
            # end for
        # end for
        pb.close()
    # end with
# end func


@click.command()
@click.argument('input-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output-file', type=click.Path(exists=False, dir_okay=False), required=True)
@click.option('--channel-order', type=click.Choice(['ZNE', 'ZRT', 'Z12', 'LQT']),
              required=False, default='ZNE', show_default=True)
def main(input_file, output_file, channel_order):
    ned = NetworkEventDataset(input_file)
    if channel_order == 'ZNE':
        ned.apply(lambda stream: stream.rotate('->ZNE'))
    elif channel_order == 'ZRT':
        # Apply curation to streams prior to rotation. Rotation can't be done on streams
        # that fail this curation.
        ned.curate(lambda _, evid, stream: curate_stream3c(evid, stream))
        ned.apply(lambda stream: stream.rotate('NE->RT'))
    elif channel_order == 'LQT':
        ned.curate(lambda _, evid, stream: curate_stream3c(evid, stream))
        ned.apply(lambda stream: stream.rotate('ZNE->LQT'))
    # end if
    plot_ned_seismograms(ned, output_file, channel_order)
# end func


if __name__ == "__main__":
    main()
# end if
