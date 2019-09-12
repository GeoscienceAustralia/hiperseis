#!/usr/bin/env python

import logging

import matplotlib.pyplot as plt

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

infile = filedialog.askopenfilename(initialdir=".", title="Select RF file", filetypes=(("h5 files", "*.h5"),))
log.info("Loading %s", infile)

data_all = rf_util.read_h5_rf(infile, network='OA', station='BV21', loc='0M')
data_dict = rf_util.rf_to_dict(data_all)

# stations = sorted(list(data_dict.keys()))
# stations = ['BV21', 'BV22']
stations = ['BV21']

# rf_type = 'ZRT'  # set manually
rf_type = data_all[0].stats.rotation

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
    fig.set_size_inches(8.27, 11.69)
    fig.suptitle("Channel {}".format(rf_stream[0].stats.channel))
    plt.show()
# end for
