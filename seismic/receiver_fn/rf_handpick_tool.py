#!/usr/bin/env python

import logging

import matplotlib.pyplot as plt
import matplotlib.patches as patches
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

data_all = rf_util.read_h5_rf(infile, network='OA', station='BV21', loc='0M')
data_dict = rf_util.rf_to_dict(data_all)

# stations = sorted(list(data_dict.keys()))
# stations = ['BV21', 'BV22']
stations = ['BV21']

# rf_type = 'ZRT'  # set manually
rf_type = data_all[0].stats.rotation

def onpick(event):
    print("in onpick")
    print(event)
    print(dir(event))
# end func

active_rect = None
down_loc = None
def on_button_down(event, ax):
    global active_rect, down_loc
    print("down at ({}, {})".format(event.xdata, event.ydata))
    down_loc = (event.xdata, event.ydata)
    active_rect = patches.Rectangle(down_loc, 0, 0, linewidth=1, edgecolor='r', facecolor='none')
    ax.add_patch(active_rect)

def on_button_up(event):
    global active_rect, down_loc
    print("up at ({}, {})".format(event.xdata, event.ydata))
    active_rect.set_visible(False)
    active_rect.remove()
    del active_rect
    down_loc = None
    event.canvas.draw()

def on_move(event):
    global active_rect, down_loc
    if down_loc is not None:
        active_rect.set_width(max(event.xdata - down_loc[0], 0))
        active_rect.set_height(max(event.ydata - down_loc[1], 0))
        active_rect.figure.canvas.draw()

# def enter_axes(event):
#     print('enter_axes', event.inaxes)
#     event.inaxes.patch.set_facecolor('yellow')
#     event.canvas.draw()

# def leave_axes(event):
#     print('leave_axes', event.inaxes)
#     event.inaxes.patch.set_facecolor('white')
#     event.canvas.draw()

# def enter_figure(event):
#     print('enter_figure', event.canvas.figure)
#     event.canvas.figure.patch.set_facecolor('red')
#     event.canvas.draw()

# def leave_figure(event):
#     print('leave_figure', event.canvas.figure)
#     event.canvas.figure.patch.set_facecolor('grey')
#     event.canvas.draw()

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
    ax0 = fig.axes[0]
    # fig.picker = True
    # ax0.picker = True
    # # line_plots = [c for c in ax0.get_children() if isinstance(c, matplotlib.lines.Line2D)]
    # # for line_artist in line_plots:
    # #     line_artist.picker = 10.0
    # for w in ax0.get_children():
    #     w.picker = 10.0

    # cid = fig.canvas.mpl_connect('pick_event', onpick)
    fig.canvas.mpl_connect('button_press_event', lambda e: on_button_down(e, ax0))
    fig.canvas.mpl_connect('button_release_event', on_button_up)
    fig.canvas.mpl_connect('motion_notify_event', on_move)
    # fig.canvas.mpl_connect('figure_enter_event', enter_figure)
    # fig.canvas.mpl_connect('figure_leave_event', leave_figure)
    # fig.canvas.mpl_connect('axes_enter_event', enter_axes)
    # fig.canvas.mpl_connect('axes_leave_event', leave_axes)
    plt.show()

    # fig.canvas.mpl_disconnect(cid)
# end for
