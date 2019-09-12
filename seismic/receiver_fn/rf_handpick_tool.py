#!/usr/bin/env python

import logging

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.widgets import RectangleSelector
import matplotlib

import tkinter as tk
from tkinter import filedialog

import rf

import seismic.receiver_fn.rf_util as rf_util
import seismic.receiver_fn.rf_plot_utils as rf_plot_utils


logging.basicConfig()
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

root = tk.Tk()
root.withdraw()

# infile = filedialog.askopenfilename(initialdir=".", title="Select RF file", filetypes=(("h5 files", "*.h5"),))
infile = r"C:\software\hiperseis\seismic\receiver_fn\DATA\OA_rfs_20170911T000036-20181128T230620_ZRT_td_rev9_qual.h5"
log.info("Loading %s", infile)

data_all = rf_util.read_h5_rf(infile, network='OA', station='BS24', loc='0M')
data_dict = rf_util.rf_to_dict(data_all)

# stations = sorted(list(data_dict.keys()))
# stations = ['BV21', 'BV22']
stations = ['BS24']

# rf_type = 'ZRT'  # set manually
rf_type = data_all[0].stats.rotation

def on_select(e_down, e_up, select_mask):
    print("on_select ({}, {}) to ({}, {})".format(e_down.xdata, e_down.ydata, e_up.xdata, e_up.ydata))
    min_y = np.round(np.min([e_down.ydata, e_up.ydata])).astype(int)
    max_y = np.round(np.max([e_down.ydata, e_up.ydata])).astype(int)
    toggle_range = (max(0, min_y - 1), min(len(select_mask), max_y))
    select_mask[toggle_range[0]:toggle_range[1]] = ~select_mask[toggle_range[0]:toggle_range[1]]
    print("Selection:")
    if np.any(select_mask):
        print(np.nonzero(select_mask)[0].tolist())
    else:
        print("none")

for st in stations:
    station_db = data_dict[st]

    # Choose RF channel
    channel = rf_util.choose_rf_source_channel(rf_type, station_db)
    channel_data = station_db[channel]

    # Label and filter quality
    rf_util.label_rf_quality_simple_amplitude(rf_type, channel_data)
    rf_stream = rf.RFStream([tr for tr in channel_data if tr.stats.predicted_quality == 'a']).sort(['back_azimuth'])
    if not rf_stream:
        log.info("No data survived filtering for %s, skipping", st)
        continue

    # Plot RF stack of primary component
    fig = rf_plot_utils.plot_rf_stack(rf_stream)
    fig.set_size_inches(8, 9)
    fig.suptitle("Channel {}".format(rf_stream[0].stats.channel))
    ax0 = fig.axes[0]

    mask = np.array([False]*len(rf_stream))
    rect_select = RectangleSelector(ax0, lambda e0, e1: on_select(e0, e1, mask), useblit=True)
    # cid = fig.canvas.mpl_connect('pick_event', onpick)
    plt.show()

    # fig.canvas.mpl_disconnect(cid)
# end for
