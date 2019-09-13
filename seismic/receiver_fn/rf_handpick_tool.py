#!/usr/bin/env python
"""Simple tool for hand picking functions from RF plots. Selected plots event IDs are exported
to text file for later use as RF filter in other tools and workflows.
"""
import os
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


def on_select(e_down, e_up, select_mask):
    """Event handler for RectangleSelector

    :param e_down: Button down event for start of rectangle
    :type e_down: matplotlib.backend_bases.MouseEvent
    :param e_up: Button up event for end of rectangle
    :type e_up: matplotlib.backend_bases.MouseEvent
    :param select_mask: Boolean mask of selected state of each receiver function in the plot
    :type select_mask: numpy.array(bool)
    """
    assert e_up.inaxes == e_down.inaxes, "Up event in different axes, xy data may be incorrect!"
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


def on_release(event, target_axes, select_mask, background, rect_selector):
    """Event handler when mouse button released.

    :param target_axes: [TBD]
    :type target_axes: [TBD]
    :param event: Button up event for end of rectangle area selection
    :type event: matplotlib.backend_bases.MouseEvent
    :param select_mask: Boolean mask of selected state of each receiver function in the plot
    :type select_mask: numpy.array(bool)
    :param background: Background raster from initial render of RFs
    :type background: Return type from fig.canvas.copy_from_bbox
    :param rect_selector: Widget for selecting rectangular region to toggle RF selection
    :type rect_selector: matplotlib.widgets.RectangleSelector
    """

    # We hardwire this to specific target axes, since the mouse release event could
    # occur in a different axes.
    ax = target_axes

    # Display selection
    event.canvas.restore_region(background)
    lines = [c for c in ax.get_children() if isinstance(c, matplotlib.lines.Line2D)]
    for i, selected in enumerate(select_mask):
        line = lines[i]
        if selected:
            line.set_color('red')
            ax.draw_artist(line)
        else:
            line.set_color('black')
        # end if
    # end for
    event.canvas.blit(ax.bbox)
    # Make sure rect_selector widget now updates it's background to the new colors
    rect_selector.update_background(event)
# end func


def main():
    """Main entry function for RF picking tool.
    """
    infile = filedialog.askopenfilename(initialdir=".", title="Select RF file", filetypes=(("h5 files", "*.h5"),))
    output_folder = filedialog.askdirectory(initialdir=os.path.split(infile)[0], title='Select output folder')
    if not os.path.isdir(output_folder):
        log.info("Creating output folder {}".format(output_folder))
        os.makedirs(output_folder, exist_ok=True)
    # end if
    log.info("Output files will be emitted to {}".format(output_folder))

    log.info("Loading %s", infile)
    data_all = rf_util.read_h5_rf(infile)
    data_dict = rf_util.rf_to_dict(data_all)

    stations = sorted(list(data_dict.keys()))

    # Assuming same rotation type for all RFs. This is consistent with the existing workflow.
    rf_type = data_all[0].stats.rotation

    for st in stations:
        station_db = data_dict[st]

        # Choose RF channel
        channel = rf_util.choose_rf_source_channel(rf_type, station_db)
        channel_data = station_db[channel]
        # Check assumption
        for tr in channel_data:
            assert tr.stats.rotation == rf_type, 'Mismatching RF rotation type'

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
        cid = fig.canvas.mpl_connect('button_release_event',
                                     lambda e: on_release(e, ax0, mask, blit_background, rect_select))
        plt.show()

        fig.canvas.mpl_disconnect(cid)
        rect_select = None

        selected_event_ids = [tr.stats.event_id for i, tr in enumerate(rf_stream) if mask[i]]
        log.info("{} streams selected".format(len(selected_event_ids)))
        log.info("Selected event ids:")
        log.info(selected_event_ids)

        network = rf_stream[0].stats.network
        outfile = os.path.join(output_folder, '.'.join([network, st, channel]) + '_event_mask.txt')
        log.info("Writing mask to file {}".format(outfile))
        if os.path.exists(outfile):
            log.warning("Overwriting existing file {} !".format(outfile))
        with open(outfile, 'w') as f:
            f.write('\n'.join(selected_event_ids))
    # end for
# end main

# -------------------------------------
if __name__ == "__main__":
    logging.basicConfig()
    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)

    root = tk.Tk()
    root.withdraw()

    main()
