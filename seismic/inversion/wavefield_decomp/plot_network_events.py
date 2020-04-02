#!/usr/bin/env python
# coding: utf-8
"""Bulk plotting helper functions based on NetworkEventDataset
"""

# pylint: disable-all

import math

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
import obspy

from seismic.inversion.wavefield_decomp.network_event_dataset import NetworkEventDataset

net = 'OA'
src_file = r"D:\temp\ga_work\OA_BT28_event_waveforms_for_rf_20170911T000036-20181128T230620_rev8.h5"

ned = NetworkEventDataset(src_file)

# ---------- BEGIN FILTERING ----------
# Rotate to ZRT coordinates

# # Apply curation to streams prior to rotation
# ned.curate(lambda _, evid, stream: curate_stream3c(evid, stream))

# Rotate to ZRT coordinates
ned.apply(lambda stream: stream.rotate('NE->RT'))

# Trim to shorter time window
time_window = (-20.0, 50.0)
ned.apply(lambda stream: stream.trim(stream[0].stats.onset + time_window[0], stream[0].stats.onset + time_window[1]))

# # Detrend the traces
# ned.apply(lambda stream: stream.detrend('linear'))

# ---------- END FILTERING ----------

# Per page
ncols = 2
nrows = 8
pagesize = (8.3, 8.3*1.414)  # A4
# channel_order = {'Z': 0, 'R': 1, 'T': 2}
channel_order = 'ZRT'

bbstyle = dict(boxstyle="round", fc="w", alpha=0.5, linewidth=0.5)
annot_fontsize = 5
axes_fontsize = 6

output_file = '{}_event_seismograms.pdf'.format(net)

with PdfPages(output_file) as pdf:
    for sta, db_evid in ned.by_station():
        num_events = len(db_evid)
        num_per_page = nrows*ncols
        ev_ids = list(db_evid.keys())
        num_pages = math.ceil(num_events/num_per_page)
        for pagenum in range(num_pages):
            f = plt.figure(constrained_layout=False, figsize=pagesize)
            f.suptitle('{}.{} event seismograms (pg {}/{})'.format(net, sta, pagenum + 1, num_pages), y=0.98, fontsize=11, va='top')
            pgspec = gridspec.GridSpec(ncols=ncols, nrows=nrows, figure=f, left=0.1, right=0.9, bottom=0.05, top=0.95, hspace=0.3, wspace=0.2)
            for i, evid in enumerate(ev_ids[pagenum*num_per_page:(pagenum + 1)*num_per_page]):
                stream = db_evid[evid]
                gs = pgspec[i%nrows, i//nrows].subgridspec(3, 1, hspace=0.0)

                # Find max range so we can use a common scale
                max_half_range = 0.0
                for tr in stream:
                    max_half_range = max(max_half_range, 0.5*(np.nanmax(tr.data) - np.nanmin(tr.data)))
                # end for
                max_half_range *= 1.1

                # Do plots
                for tr in stream:
                    cha = tr.stats.channel[-1]
                    idx = channel_order.index(cha)
                    axn = f.add_subplot(gs[idx])
                    t = tr.times() - (tr.stats.onset - tr.stats.starttime)
                    tr_mean = np.nanmean(tr.data)
                    axn.plot(t, tr.data, 'k')
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
# end with
