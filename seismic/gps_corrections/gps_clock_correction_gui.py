#!/bin/env python
"""
GUI interface to generating clock corrections from x-corr results.
"""

import os

try:
    import tkinter as tk
except ImportError:
    import Tkinter as tk

import tkFileDialog

import matplotlib
matplotlib.use("TkAgg")
# import matplotlib.backends.tkagg as tkagg
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
# from matplotlib.figure import Figure

from seismic.ASDFdatabase import FederatedASDFDataSet
from seismic.xcorqc.xcorr_station_clock_analysis import (plot_xcorr_file_clock_analysis,
                                                         read_correlator_config,
                                                         batch_process_xcorr)


dataset = "/g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt"


class GpsClockCorrectionApp(tk.Frame):

    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        self.nc_file = tk.StringVar(self)
        self.pack()
        self._last_dir = None
        self._createStep0Widgets()
        self.fds = FederatedASDFDataSet.FederatedASDFDataSet(dataset)
        self.current_step = 0
        self.station_code = ''
        self.time_window = tk.IntVar(self, value=300)
        self.snr_threshold = tk.IntVar(self, value=6)
        self.xcorr_fig = None
        self.fig_canv = None

    def _createStep0Widgets(self):
        self.TOP_PANE_0 = tk.PanedWindow(self, orient=tk.VERTICAL)
        self.TOP_PANE_0.pack(fill=tk.BOTH, expand=1)

        self.UPPER_PANE_0 = tk.PanedWindow(self.TOP_PANE_0, orient=tk.HORIZONTAL)
        self.TOP_PANE_0.add(self.UPPER_PANE_0)

        self.STN_CODE_LABEL = tk.Label(self.TOP_PANE_0, text="Station code: ")
        self.TOP_PANE_0.add(self.STN_CODE_LABEL)

        self.LOWER_PANE_0 = tk.PanedWindow(self.TOP_PANE_0, orient=tk.HORIZONTAL)
        self.TOP_PANE_0.add(self.LOWER_PANE_0)

        self.NC_FILE_LABEL = tk.Label(self.UPPER_PANE_0, text="Cross-correlation file:")
        self.NC_FILE_LABEL.pack(side=tk.LEFT)

        self.NC_FILE_ENTRY = tk.Entry(self.UPPER_PANE_0)
        self.NC_FILE_ENTRY['width'] = 64
        self.NC_FILE_ENTRY.pack(padx=5, side=tk.LEFT)
        self.NC_FILE_ENTRY['textvariable'] = self.nc_file

        self.OPEN = tk.Button(self.LOWER_PANE_0)
        self.OPEN['text'] = "Open..."
        self.OPEN['command'] = self._openNcFile
        self.OPEN.pack(side=tk.LEFT)

        self.NEXT = tk.Button(self.LOWER_PANE_0)
        self.NEXT['text'] = "Next..."
        self.NEXT['state'] = tk.DISABLED
        self.NEXT['command'] = self._gotoStep1
        self.NEXT.pack(side=tk.LEFT)
        self.nc_file.trace_variable('w', self._updateNextButton)

        self.QUIT = tk.Button(self.LOWER_PANE_0)
        self.QUIT['text'] = "Quit"
        self.QUIT['command'] = self.quit
        self.QUIT.pack(side=tk.RIGHT)

    def _extractCodeFromFilename(self, src_file):
        _, basename = os.path.split(src_file)
        _, file_type = os.path.splitext(src_file)
        name_parts = basename.split('.')
        netcode = name_parts[0]
        statcode = name_parts[1]
        full_code = '.'.join([netcode, statcode])
        return full_code

    def _updateNextButton(self, _1, _2, _3):
        nc_file_value = self.nc_file.get()
        nc_file_valid = bool(nc_file_value) and os.path.exists(nc_file_value)
        self.NEXT['state'] = tk.NORMAL if nc_file_valid else tk.DISABLED
        if nc_file_valid:
            self.station_code= self._extractCodeFromFilename(nc_file_value)
        else:
            self.station_code = ''
        self.STN_CODE_LABEL['text'] = "Station code: " + self.station_code

    def _openNcFile(self):
        initial_dir = '/g/data/ha3/Passive/SHARED_DATA/GPS_Clock/xcorr' if self._last_dir is None else self._last_dir
        file_name = tkFileDialog.askopenfilename(initialdir=initial_dir,
                                                 title='Select .nc file to analyze',
                                                 filetypes=(("nc files", "*.nc"), ("all files", "*")))
        if file_name:
            self._last_dir = os.path.split(file_name)[0]
            self.nc_file.set(file_name)

    def _gotoStep1(self):
        if self.current_step == 0:
            # for child in self.TOP_PANE_0.winfo_children():
            #     self.child.destroy()
            self.TOP_PANE_0.destroy()
            self.TOP_PANE_0 = None
            self._createStep1Widgets()

    def _createStep1Widgets(self):
        self.current_step = 1
        self.TOP_PANE_1 = tk.PanedWindow(self, orient=tk.HORIZONTAL)
        self.TOP_PANE_1.pack(fill=tk.X, expand=0)

        self.LEFT_PANE_1 = tk.PanedWindow(self.TOP_PANE_1, orient=tk.VERTICAL, width=500)
        self.LEFT_PANE_1.pack(fill=tk.X, expand=0)
        self.TOP_PANE_1.add(self.LEFT_PANE_1)

        self.STN_CODE_LABEL = tk.Label(self.LEFT_PANE_1, text="Station code: " + self.station_code)
        self.STN_CODE_LABEL.pack(side=tk.LEFT)
        self.LEFT_PANE_1.add(self.STN_CODE_LABEL)

        self.TIME_WINDOW_PANE = tk.PanedWindow(self.LEFT_PANE_1, orient=tk.VERTICAL, borderwidth=2, relief=tk.RAISED)
        self.TIME_WINDOW_PANE.pack(fill=tk.X, padx=2, pady=2)
        self.LEFT_PANE_1.add(self.TIME_WINDOW_PANE)

        self.TIME_WINDOW_LABEL = tk.Label(self.TIME_WINDOW_PANE, text="Time window (s)")
        self.TIME_WINDOW_LABEL.pack(side=tk.LEFT)
        self.TIME_WINDOW_PANE.add(self.TIME_WINDOW_LABEL)
        self.TIME_WINDOW_ENTRY = tk.Scale(self.TIME_WINDOW_PANE, from_=30, to=1800, resolution=30,
                                          orient=tk.HORIZONTAL, width=200)
        self.TIME_WINDOW_ENTRY.set(self.time_window.get())
        self.TIME_WINDOW_ENTRY['variable'] = self.time_window
        self.TIME_WINDOW_ENTRY['command'] = self._enableRefresh
        self.TIME_WINDOW_ENTRY.pack(fill=tk.X)
        self.TIME_WINDOW_PANE.add(self.TIME_WINDOW_ENTRY)

        self.SNR_PANE = tk.PanedWindow(self.LEFT_PANE_1, orient=tk.VERTICAL, borderwidth=2, relief=tk.RAISED)
        self.SNR_PANE.pack(fill=tk.X, padx=2, pady=2)
        self.LEFT_PANE_1.add(self.SNR_PANE)

        self.SNR_LABEL = tk.Label(self.SNR_PANE, text="SNR")
        self.SNR_LABEL.pack(side=tk.LEFT)
        self.SNR_PANE.add(self.SNR_LABEL)
        self.SNR_ENTRY = tk.Scale(self.SNR_PANE, from_=0, to=100, resolution=1, orient=tk.HORIZONTAL, width=200)
        self.SNR_ENTRY.set(self.snr_threshold.get())
        self.SNR_ENTRY['variable'] = self.snr_threshold
        self.SNR_ENTRY['command'] = self._enableRefresh
        self.SNR_ENTRY.pack(fill=tk.X)
        self.SNR_PANE.add(self.SNR_ENTRY)

        self.RIGHT_FIGURE_CANVAS_1 = tk.Canvas(self.TOP_PANE_1)
        self.RIGHT_FIGURE_CANVAS_1.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.TOP_PANE_1.add(self.RIGHT_FIGURE_CANVAS_1)

        self.REFRESH = tk.Button(self.LEFT_PANE_1)
        self.REFRESH['text'] = "Refresh"
        self.REFRESH['state'] = tk.DISABLED
        self.REFRESH['command'] = self._updateStep1Canvas
        self.REFRESH.pack(side=tk.TOP, fill=tk.X)
        self.LEFT_PANE_1.add(self.REFRESH)

        self.NEXT = tk.Button(self.LEFT_PANE_1)
        self.NEXT['text'] = "Next..."
        self.NEXT['command'] = self._saveStep1AndContinue
        self.NEXT.pack(side=tk.TOP, fill=tk.X)
        self.LEFT_PANE_1.add(self.NEXT)

        self.QUIT = tk.Button(self.LEFT_PANE_1)
        self.QUIT['text'] = "Quit"
        self.QUIT['command'] = self.quit
        self.QUIT.pack(side=tk.BOTTOM, fill=tk.X)
        self.LEFT_PANE_1.add(self.QUIT)

        self.xcorr_settings, self.xcorr_title_tag = read_correlator_config(self.nc_file.get())
        self._updateStep1Canvas()

    def _enableRefresh(self, _new_val):
        self.REFRESH['state'] = tk.NORMAL

    def _updateStep1Canvas(self):
        self.REFRESH['state'] = tk.DISABLED
        self.NEXT['state'] = tk.DISABLED
        self.TOP_PANE_1.remove(self.RIGHT_FIGURE_CANVAS_1)
        if self.xcorr_fig is not None:
            self.xcorr_fig.clear()
            del self.xcorr_fig
        if self.fig_canv is not None:
            self.fig_canv.get_tk_widget().destroy()
        self.RIGHT_FIGURE_CANVAS_1.delete(tk.ALL)
        self.xcorr_fig = plot_xcorr_file_clock_analysis(self.nc_file.get(), self.fds, self.time_window.get(),
                                                        self.snr_threshold.get(), show=False,
                                                        title_tag=self.xcorr_title_tag, settings=self.xcorr_settings)
        self.fig_canv = FigureCanvasTkAgg(self.xcorr_fig, master=self.RIGHT_FIGURE_CANVAS_1)
        self.fig_canv.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.RIGHT_FIGURE_CANVAS_1.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.TOP_PANE_1.add(self.RIGHT_FIGURE_CANVAS_1)
        self.NEXT['state'] = tk.NORMAL
        self.update()

    def _saveStep1AndContinue(self):
        if self.current_step == 1:
            self.NEXT['state'] = tk.DISABLED
            batch_process_xcorr([self.nc_file.get()], self.fds, self.time_window.get(), self.snr_threshold.get(),
                                save_plots=True, force_save=True)
            self._gotoStep2()

    def _gotoStep2(self):
        if self.current_step == 1:
            # for child in self.TOP_PANE_0.winfo_children():
            #     self.child.destroy()
            if self.xcorr_fig is not None:
                self.xcorr_fig.clear()
                del self.xcorr_fig
                self.xcorr_fig = None
            if self.fig_canv is not None:
                self.fig_canv.get_tk_widget().destroy()
                self.fig_canv = None
            self.RIGHT_FIGURE_CANVAS_1.delete(tk.ALL)
            self.TOP_PANE_1.destroy()
            self.TOP_PANE_1 = None
            self._createStep2Widgets()

    def _createStep2Widgets(self):
        self.current_step = 2

#end class


tk_root = tk.Tk()
app = GpsClockCorrectionApp(master=tk_root)
app.master.title("GPS Clock Correction Workflow")
app.master.minsize(320, 160)
app.master.columnconfigure(0, weight=1)
app.master.rowconfigure(0, weight=1)
app.mainloop()
# tk_root.destroy()
