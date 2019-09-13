#!/usr/bin/env python

import logging
import tkinter as tk
from tkinter import filedialog

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
import matplotlib

import rf

import seismic.receiver_fn.rf_util as rf_util
import seismic.receiver_fn.rf_plot_utils as rf_plot_utils

# pylint: disable=invalid-name, logging-format-interpolation

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
    """Event handler for RectangleSelector

    :param e_down: [description]
    :type e_down: [type]
    :param e_up: [description]
    :type e_up: [type]
    :param select_mask: [description]
    :type select_mask: [type]
    """
    log.debug("on_select ({}, {}) to ({}, {})".format(e_down.xdata, e_down.ydata, e_up.xdata, e_up.ydata))
    min_y = np.round(np.min([e_down.ydata, e_up.ydata])).astype(int)
    max_y = np.round(np.max([e_down.ydata, e_up.ydata])).astype(int)
    toggle_range = (max(0, min_y - 1), min(len(select_mask), max_y))
    select_mask[toggle_range[0]:toggle_range[1]] = ~select_mask[toggle_range[0]:toggle_range[1]]

    # Echo selection
    log.info("Selection:")
    if np.any(select_mask):
        log.info(np.nonzero(select_mask)[0].tolist())
    else:
        log.info("none")
    # end if

# end func

def on_release(event, select_mask, background, rect_selector):
    """Event handler when mouse button released.

    :param event: [description]
    :type event: [type]
    :param select_mask: [description]
    :type select_mask: [type]
    :param background: [description]
    :type background: [type]
    :param rect_selector: [description]
    :type rect_selector: [type]
    """

    # Display selection
    event.canvas.restore_region(background)
    lines = [c for c in event.inaxes.get_children() if isinstance(c, matplotlib.lines.Line2D)]
    for i, selected in enumerate(select_mask):
        line = lines[i]
        if selected:
            line.set_color('red')
            event.inaxes.draw_artist(line)
        else:
            line.set_color('black')
        # end if
    # end for
    event.canvas.blit(event.inaxes.bbox)
    # Make sure rect_selector widget now updates it's background to the new colors
    rect_selector.update_background(event)

# end func


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
    # Make sure we draw once first before capturing blit background
    fig.canvas.draw()
    # Disallow resizing to avoid having to handle blitting with resized window.
    win = fig.canvas.window()
    win.setFixedSize(win.size())
    blit_background = fig.canvas.copy_from_bbox(ax0.bbox)

    mask = np.array([False]*len(rf_stream))
    rect_select = RectangleSelector(ax0, lambda e0, e1: on_select(e0, e1, mask), useblit=True,
                                    rectprops=dict(fill=False, edgecolor='red'))
    cid = fig.canvas.mpl_connect('button_release_event', lambda e: on_release(e, mask, blit_background, rect_select))
    plt.show()

    fig.canvas.mpl_disconnect(cid)
    rect_select = None
# end for
