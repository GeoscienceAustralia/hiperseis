#!/bin/env python
"""
GUI interface to generating clock corrections from x-corr results.
"""

import os
import copy

try:
    import Tkinter as tk
except ImportError:
    import tkinter as tk

try:
    import tkFileDialog as filedialog
    import tkMessageBox as messagebox
except ImportError:
    from tkinter import filedialog
    from tkinter import messagebox

import matplotlib
matplotlib.use("TkAgg")
# import matplotlib.backends.tkagg as tkagg
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
# from matplotlib.figure import Figure

import numpy as np
import pandas as pd
import obspy

from seismic.ASDFdatabase import FederatedASDFDataSet
from seismic.xcorqc.xcorr_station_clock_analysis import (plot_xcorr_file_clock_analysis,
                                                         read_correlator_config,
                                                         batch_process_xcorr)


dataset = "/g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt"


class GpsClockCorrectionApp(tk.Frame):

    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        self.nc_file = tk.StringVar(self)
        self.output_folder = tk.StringVar(self)
        self.fds = FederatedASDFDataSet.FederatedASDFDataSet(dataset)
        self.current_step = 0
        self.station_code = ''
        self.time_window = tk.IntVar(self, value=300)
        self.snr_threshold = tk.IntVar(self, value=6)
        self.pearson_cutoff_factor = 0.5
        self.xcorr_ca = None
        self.xcorr_fig = None
        self.fig_canv = None
        self.cluster_fig = None
        self.cluster_fig_canv = None
        self.cluster_coeff0 = tk.DoubleVar(self, value=1.0)
        self.cluster_coeff1 = tk.DoubleVar(self, value=1.0)
        self.cluster_coeff2 = tk.DoubleVar(self, value=0.0)
        self.cluster_ids = None
        self.selected_cluster_ids = None
        self.degrees = None
        self.resampling_period_days = tk.DoubleVar(self, value=1.0)
        self.regression_fig = None
        self.regression_fig_canv = None
        self.resampling_fig = None
        self.resampling_fig_canv = None
        self.display_dpi = 50
        self.pack()

        self._createStep0Widgets()

    def busy(self):
        tk_root.config(cursor="watch")
        tk_root.update()

    def not_busy(self):
        tk_root.config(cursor="")

    def _createStep0Widgets(self):
        self.ROOT_FRAME_0 = tk.Frame(self)
        self.ROOT_FRAME_0.pack(fill=tk.BOTH, expand=1)

        self.UPPER_FRAME_0 = tk.LabelFrame(self.ROOT_FRAME_0, text="File selection", borderwidth=2)
        self.UPPER_FRAME_0.pack(anchor=tk.N, side=tk.TOP, fill=tk.X, padx=2, pady=2)

        self.STN_CODE_LABEL = tk.Label(self.ROOT_FRAME_0, text="Station code: ", font=8)
        self.STN_CODE_LABEL.pack(side=tk.TOP, fill=tk.X, padx=2, pady=8)

        self.LOWER_FRAME_0 = tk.Frame(self.ROOT_FRAME_0)
        self.LOWER_FRAME_0.pack(anchor=tk.S, side=tk.BOTTOM, fill=tk.X, padx=2, pady=2)

        file_entry_frame = tk.Frame(self.UPPER_FRAME_0)
        file_entry_frame.pack(anchor=tk.E, side=tk.TOP, fill=tk.X)
        self.OPEN = tk.Button(file_entry_frame, width=10)
        self.OPEN['text'] = "Open..."
        self.OPEN['command'] = self._openNcFile
        self.OPEN.pack(anchor=tk.W, side=tk.LEFT, padx=2, pady=2)

        self.NC_FILE_LABEL = tk.Label(file_entry_frame, text="Cross-correlation file:")
        self.NC_FILE_LABEL.pack(anchor=tk.E, side=tk.LEFT, padx=(8,0), pady=2)

        self.NC_FILE_ENTRY = tk.Entry(file_entry_frame, exportselection=False, width=64)
        self.NC_FILE_ENTRY.pack(anchor=tk.E, padx=(2,4), pady=2, side=tk.TOP, fill=tk.X)
        self.NC_FILE_ENTRY['textvariable'] = self.nc_file

        output_folder_frame = tk.Frame(self.UPPER_FRAME_0)
        output_folder_frame.pack(anchor=tk.E, side=tk.TOP, fill=tk.X)
        self.SELECT = tk.Button(output_folder_frame, width=10)
        self.SELECT['text'] = "Choose..."
        self.SELECT['command'] = self._chooseOutputFolder
        self.SELECT.pack(anchor=tk.W, side=tk.LEFT, padx=2, pady=2)

        self.OUTPUT_FOLDER_LABEL = tk.Label(output_folder_frame, text="Output folder:")
        self.OUTPUT_FOLDER_LABEL.pack(anchor=tk.E, side=tk.LEFT, padx=(8,0), pady=2)

        self.OUTPUT_FOLDER_ENTRY = tk.Entry(output_folder_frame, exportselection=False, width=64)
        self.OUTPUT_FOLDER_ENTRY.pack(anchor=tk.E, padx=(2,4), pady=2, side=tk.TOP, fill=tk.X)
        self.OUTPUT_FOLDER_ENTRY['textvariable'] = self.output_folder

        self.NEXT = tk.Button(self.LOWER_FRAME_0)
        self.NEXT['text'] = "Continue..."
        self.NEXT['state'] = tk.DISABLED
        self.NEXT['command'] = self._gotoStep1
        self.NEXT.pack(anchor=tk.SW, side=tk.LEFT)
        self.nc_file.trace_variable('w', self._updateNextButton)
        self.output_folder.trace_variable('w', self._updateNextButton)

        self.QUIT = tk.Button(self.LOWER_FRAME_0)
        self.QUIT['text'] = "Quit"
        self.QUIT['command'] = self._quitApp
        self.QUIT.pack(anchor=tk.SE, side=tk.RIGHT)

    def _quitApp(self):
        self._destroyFigures()
        self.quit()

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
        nc_file_valid = bool(nc_file_value) and os.path.isfile(nc_file_value)
        if nc_file_valid:
            self.station_code= self._extractCodeFromFilename(nc_file_value)
        else:
            self.station_code = ''
        self.STN_CODE_LABEL['text'] = "Station code: " + self.station_code
        output_folder_value = self.output_folder.get()
        output_folder_valid = bool(output_folder_value) and os.path.isdir(output_folder_value)
        self.NEXT['state'] = tk.NORMAL if (nc_file_valid and output_folder_valid) else tk.DISABLED

    def _openNcFile(self):
        initial_dir = self.nc_file.get()
        if not initial_dir:
            initial_dir = '/g/data/ha3/Passive/SHARED_DATA/GPS_Clock/xcorr'
        file_name = filedialog.askopenfilename(initialdir=initial_dir,
                                               title='Select .nc file to analyze',
                                               filetypes=(("nc files", "*.nc"), ("all files", "*")))
        if file_name:
            self.nc_file.set(file_name)

    def _chooseOutputFolder(self):
        initial_dir = self.output_folder.get()
        if not initial_dir:
            initial_dir = '/g/data/ha3/Passive/SHARED_DATA/GPS_Clock/corrections'
        folder = filedialog.askdirectory(initialdir=initial_dir, title='Select output folder', mustexist=True)
        if folder:
            self.output_folder.set(folder)

    def _gotoStep1(self):
        if self.current_step == 0:
            # for child in self.TOP_PANE_0.winfo_children():
            #     self.child.destroy()
            self.ROOT_FRAME_0.destroy()
            self.ROOT_FRAME_0 = None
            self._createStep1Widgets()

    def _createStep1Widgets(self):
        self.busy()
        self.current_step = 1
        self.ROOT_FRAME_1 = tk.Frame(self)
        self.ROOT_FRAME_1.pack(fill=tk.BOTH, expand=1)

        self.LEFT_FRAME_1 = tk.LabelFrame(self.ROOT_FRAME_1, width=800)
        self.LEFT_FRAME_1.pack(anchor=tk.NW, side=tk.LEFT, fill=tk.X, padx=2, pady=2)

        self.STN_CODE_LABEL = tk.Label(self.LEFT_FRAME_1, text="Station code: " + self.station_code, font=8)
        self.STN_CODE_LABEL.pack(anchor=tk.W, side=tk.TOP)

        self.TIME_WINDOW_FRAME = tk.LabelFrame(self.LEFT_FRAME_1, borderwidth=2)
        self.TIME_WINDOW_FRAME.pack(anchor=tk.W, side=tk.TOP, fill=tk.X, padx=2, pady=2)

        self.TIME_WINDOW_LABEL = tk.Label(self.TIME_WINDOW_FRAME, text="Time window (s)")
        self.TIME_WINDOW_LABEL.pack(anchor=tk.NW, side=tk.TOP)
        self.TIME_WINDOW_ENTRY = tk.Scale(self.TIME_WINDOW_FRAME, from_=30, to=1800, resolution=10,
                                          orient=tk.HORIZONTAL, length=500)
        self.TIME_WINDOW_ENTRY.set(self.time_window.get())
        self.TIME_WINDOW_ENTRY['variable'] = self.time_window
        self.TIME_WINDOW_ENTRY['command'] = self._enableRefresh
        self.TIME_WINDOW_ENTRY.pack(anchor=tk.SW, side=tk.BOTTOM, fill=tk.X)

        self.SNR_FRAME = tk.LabelFrame(self.LEFT_FRAME_1, borderwidth=2)
        self.SNR_FRAME.pack(anchor=tk.W, side=tk.TOP, fill=tk.X, padx=2, pady=2)

        self.SNR_LABEL = tk.Label(self.SNR_FRAME, text="SNR")
        self.SNR_LABEL.pack(anchor=tk.NW, side=tk.TOP)
        self.SNR_ENTRY = tk.Scale(self.SNR_FRAME, from_=0, to=100, resolution=1, orient=tk.HORIZONTAL, length=500)
        self.SNR_ENTRY.set(self.snr_threshold.get())
        self.SNR_ENTRY['variable'] = self.snr_threshold
        self.SNR_ENTRY['command'] = self._enableRefresh
        self.SNR_ENTRY.pack(anchor=tk.SW, side=tk.BOTTOM, fill=tk.X)

        self.RIGHT_FIGURE_CANVAS_1 = tk.Canvas(self.ROOT_FRAME_1)
        self.RIGHT_FIGURE_CANVAS_1.pack(anchor=tk.NE, side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=4, pady=4)

        self.REFRESH = tk.Button(self.LEFT_FRAME_1)
        self.REFRESH['text'] = "Refresh"
        self.REFRESH['state'] = tk.DISABLED
        self.REFRESH['command'] = self._updateStep1Canvas
        self.REFRESH.pack(anchor=tk.NW, side=tk.TOP, fill=tk.X, padx=2, pady=2)

        self.NEXT = tk.Button(self.LEFT_FRAME_1)
        self.NEXT['text'] = "Continue..."
        self.NEXT['command'] = self._gotoStep2
        self.NEXT.pack(anchor=tk.W, side=tk.LEFT, padx=2, pady=(16,2))

        self.QUIT = tk.Button(self.LEFT_FRAME_1)
        self.QUIT['text'] = "Quit"
        self.QUIT['command'] = self._quitApp
        self.QUIT.pack(anchor=tk.E, side=tk.RIGHT, padx=2, pady=(16,2))

        self.xcorr_settings, self.xcorr_title_tag = read_correlator_config(self.nc_file.get())
        self._updateStep1Canvas()
        self.not_busy()

    def _enableRefresh(self, _new_val):
        self.REFRESH['state'] = tk.NORMAL

    def _destroyFigures(self):
        if self.xcorr_fig is not None:
            self.xcorr_fig.clear()
            self.xcorr_fig = None
        if self.fig_canv is not None:
            self.fig_canv.get_tk_widget().destroy()
            self.fig_canv = None
        if self.cluster_fig is not None:
            self.cluster_fig.clear()
            self.cluster_fig = None
        if self.cluster_fig_canv is not None:
            self.cluster_fig_canv.get_tk_widget().destroy()
            self.cluster_fig_canv = None
        if self.regression_fig is not None:
            self.regression_fig.clear()
            self.regression_fig = None
        if self.regression_fig_canv is not None:
            self.regression_fig_canv.get_tk_widget().destroy()
            self.regression_fig_canv = None
        if self.resampling_fig is not None:
            self.resampling_fig.clear()
            self.resampling_fig = None
        if self.resampling_fig_canv is not None:
            self.resampling_fig_canv.get_tk_widget().destroy()
            self.resampling_fig_canv = None

    def _updateStep1Canvas(self):
        self.REFRESH['state'] = tk.DISABLED
        self.NEXT['state'] = tk.DISABLED
        self._destroyFigures()
        self.RIGHT_FIGURE_CANVAS_1.delete(tk.ALL)
        self.xcorr_ca, self.xcorr_fig = \
            plot_xcorr_file_clock_analysis(self.nc_file.get(), self.fds, self.time_window.get(),
                                           self.snr_threshold.get(), self.pearson_cutoff_factor, show=False,
                                           title_tag=self.xcorr_title_tag, settings=self.xcorr_settings)
        self.fig_canv = FigureCanvasTkAgg(self.xcorr_fig, master=self.RIGHT_FIGURE_CANVAS_1)
        self.fig_canv.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.RIGHT_FIGURE_CANVAS_1.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.NEXT['state'] = tk.NORMAL
        self.update()

    def _gotoStep2(self):
        if self.current_step == 1:
            self.NEXT['state'] = tk.DISABLED

            self._destroyFigures()

            self.RIGHT_FIGURE_CANVAS_1.delete(tk.ALL)
            self.ROOT_FRAME_1.destroy()
            self.ROOT_FRAME_1 = None
            self.update()

            info_label = tk.Label(self, text="Saving plot...")
            info_label.pack()
            self.update()

            # Generate PNG file for the .nc file using the current settings.
            self.busy()
            batch_process_xcorr([self.nc_file.get()], self.fds, self.time_window.get(), self.snr_threshold.get(),
                                pearson_cutoff_factor=self.pearson_cutoff_factor, save_plots=True, force_save=True)
            self.not_busy()

            info_label.destroy()

            self._createStep2Widgets()

    def _createStep2Widgets(self):
        self.current_step = 2

        assert self.xcorr_ca is not None

        self.ROOT_FRAME_2 = tk.Frame(self)
        self.ROOT_FRAME_2.pack(fill=tk.BOTH, expand=1)

        self.UPPER_FRAME_2 = tk.LabelFrame(self.ROOT_FRAME_2, text="Station code: " + self.station_code)
        self.UPPER_FRAME_2.pack(anchor=tk.N, side=tk.TOP)

        self.COEFF_WIDGETS_FRAME_2 = tk.Frame(self.UPPER_FRAME_2)
        self.COEFF_WIDGETS_FRAME_2.pack(anchor=tk.NW, side=tk.TOP, padx=2, pady=2)

        self.COEFF_WIDGETS_LABEL = tk.Label(self.COEFF_WIDGETS_FRAME_2, font=8,
                                            text="Tune coefficients to optimize clustering and filtering of outliers")
        self.COEFF_WIDGETS_LABEL.pack(anchor=tk.W, side=tk.TOP, pady=2)

        self.COEFF0_LABEL = tk.LabelFrame(self.COEFF_WIDGETS_FRAME_2, text="Coefficient 0 (x-separation)")
        self.COEFF0_LABEL.pack(anchor=tk.W, side=tk.LEFT, padx=2)
        self.COEFF0_SPINBOX = tk.Spinbox(self.COEFF0_LABEL, from_=0.0, to=100.0, increment=0.1,
                                         textvariable=self.cluster_coeff0)
        self.COEFF0_SPINBOX.pack(padx=2, pady=2)

        self.COEFF1_LABEL = tk.LabelFrame(self.COEFF_WIDGETS_FRAME_2, text="Coefficient 1 (y-separation)")
        self.COEFF1_LABEL.pack(anchor=tk.W, side=tk.LEFT, padx=2)
        self.COEFF1_SPINBOX = tk.Spinbox(self.COEFF1_LABEL, from_=0.0, to=100.0, increment=0.1,
                                         textvariable=self.cluster_coeff1)
        self.COEFF1_SPINBOX.pack(padx=2, pady=2)

        self.COEFF2_LABEL = tk.LabelFrame(self.COEFF_WIDGETS_FRAME_2, text="Coefficient 2 (slope)")
        self.COEFF2_LABEL.pack(anchor=tk.W, side=tk.LEFT, padx=2)
        self.COEFF2_SPINBOX = tk.Spinbox(self.COEFF2_LABEL, from_=0.0, to=100.0, increment=0.1,
                                         textvariable=self.cluster_coeff2)
        self.COEFF2_SPINBOX.pack(padx=2, pady=2)

        self.CONTROL_BUTTONS_FRAME_2 = tk.Frame(self.UPPER_FRAME_2)
        self.CONTROL_BUTTONS_FRAME_2.pack(anchor=tk.NW, side=tk.TOP, padx=2, pady=2, fill=tk.X)

        self.NEXT = tk.Button(self.CONTROL_BUTTONS_FRAME_2)
        self.NEXT['text'] = "Continue..."
        self.NEXT['command'] = self._gotoStep3
        self.NEXT.pack(anchor=tk.NW, side=tk.LEFT)

        self.QUIT = tk.Button(self.CONTROL_BUTTONS_FRAME_2)
        self.QUIT['text'] = "Quit"
        self.QUIT['command'] = self._quitApp
        self.QUIT.pack(anchor=tk.SE, side=tk.RIGHT)

        self.LOWER_FRAME_2 = tk.LabelFrame(self.ROOT_FRAME_2, text="Clustering Result")
        self.LOWER_FRAME_2.pack(anchor=tk.N, side=tk.TOP)

        self.CLUSTER_FIGURE_CANVAS_2 = tk.Canvas(self.LOWER_FRAME_2)
        self.CLUSTER_FIGURE_CANVAS_2.pack(anchor=tk.N, side=tk.TOP, fill=tk.BOTH, expand=True, padx=4, pady=4)

        self._refreshClusteringCanvas()

        self.update()

        self.COEFF0_SPINBOX['command'] = self._refreshClusteringCanvas
        self.COEFF1_SPINBOX['command'] = self._refreshClusteringCanvas
        self.COEFF2_SPINBOX['command'] = self._refreshClusteringCanvas

    def _refreshClusteringCanvas(self):
        if self.cluster_fig_canv:
            self.cluster_fig_canv.get_tk_widget().destroy()
        self.CLUSTER_FIGURE_CANVAS_2.delete(tk.ALL)
        self._redrawClusteringFigure()
        self.cluster_fig_canv = FigureCanvasTkAgg(self.cluster_fig, master=self.CLUSTER_FIGURE_CANVAS_2)
        self.cluster_fig_canv.get_tk_widget().pack(anchor=tk.N, side=tk.TOP, fill=tk.BOTH, expand=True)
        self.CLUSTER_FIGURE_CANVAS_2.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.update()

    def _redrawClusteringFigure(self):
        cluster_coeffs = (self.cluster_coeff0.get(), self.cluster_coeff1.get(), self.cluster_coeff2.get())
        _, self.cluster_ids = self.xcorr_ca.do_clustering(cluster_coeffs)

        if self.cluster_fig is not None:
            self.cluster_fig.clear()
        else:
            self.cluster_fig = plt.figure(figsize=(16,9))

        self.xcorr_ca.plot_clusters(self.cluster_fig.gca(), self.cluster_ids, cluster_coeffs, self.station_code)
        self.cluster_fig.tight_layout()
        self.cluster_fig.autofmt_xdate()

    def _gotoStep3(self):
        if self.current_step == 2:
            self.NEXT['state'] = tk.DISABLED

            cluster_fig_file = os.path.join(self.output_folder.get(), self.station_code + "_clustering_profile.png")
            self.cluster_fig.savefig(cluster_fig_file, dpi=300)

            assert self.cluster_ids is not None
            num_clusters = len(set(self.cluster_ids[self.cluster_ids != -1]))
            self.degrees = dict(zip(range(num_clusters), [1]*num_clusters))

            self._destroyFigures()

            self.CLUSTER_FIGURE_CANVAS_2.delete(tk.ALL)
            self.ROOT_FRAME_2.destroy()
            self.ROOT_FRAME_2 = None
            self.update()

            self._createStep3Widgets()

    def _createStep3Widgets(self):
        self.current_step = 3

        assert self.degrees is not None

        self.ROOT_FRAME_3 = tk.Frame(self)
        self.ROOT_FRAME_3.pack(fill=tk.BOTH, expand=1)

        self.LEFT_FRAME_3 = tk.LabelFrame(self.ROOT_FRAME_3, text="Station code: " + self.station_code)
        self.LEFT_FRAME_3.pack(anchor=tk.NW, side=tk.LEFT, fill=tk.X, padx=2, pady=2)

        self.DEGREE_CONTROLS = []
        for i, d in self.degrees.items():
            dc = SplineDegreeWidget(i, d, self.LEFT_FRAME_3, self._refreshSplineCanvas)
            dc.pack(anchor=tk.NW, side=tk.TOP, fill=tk.X, padx=4)
            self.DEGREE_CONTROLS.append(dc)

        self.RESAMPLING_PERIOD_LABEL = tk.Label(self.LEFT_FRAME_3, text="Enter resampling period (days):")
        self.RESAMPLING_PERIOD_LABEL.pack(anchor=tk.W, padx=2, pady=2, side=tk.TOP)

        self.RESAMPLING_PERIOD_ENTRY = tk.Entry(self.LEFT_FRAME_3, textvariable=self.resampling_period_days,
                                                exportselection=False)
        self.RESAMPLING_PERIOD_ENTRY['width'] = 8
        self.RESAMPLING_PERIOD_ENTRY.pack(anchor=tk.E, padx=5, pady=2, side=tk.TOP)
        self.resampling_period_days.trace_variable('w', self._resampleRateChanged)

        self.EXPORT = tk.Button(self.LEFT_FRAME_3)
        self.EXPORT['text'] = "Export..."
        self.EXPORT['command'] = self._exportCorrections
        self.EXPORT.pack(anchor=tk.SW, side=tk.LEFT, padx=2, pady=(16, 2))

        self.QUIT = tk.Button(self.LEFT_FRAME_3)
        self.QUIT['text'] = "Quit"
        self.QUIT['command'] = self._quitApp
        self.QUIT.pack(anchor=tk.SE, side=tk.RIGHT, padx=2, pady=(16, 2))

        self.RIGHT_FRAME_3 = tk.Frame(self.ROOT_FRAME_3)
        self.RIGHT_FRAME_3.pack(anchor=tk.NE, side=tk.RIGHT, fill=tk.BOTH)

        self.RIGHT_UPPER_FRAME_3 = tk.LabelFrame(self.RIGHT_FRAME_3, text="Regression curves overlay")
        self.RIGHT_UPPER_FRAME_3.pack(anchor=tk.NE, side=tk.TOP, fill=tk.BOTH, padx=2, pady=2)
        self.CURVE_FIGURE_CANVAS_3 = tk.Canvas(self.RIGHT_UPPER_FRAME_3)
        self.CURVE_FIGURE_CANVAS_3.pack(anchor=tk.N, side=tk.TOP, fill=tk.BOTH, expand=True, padx=2, pady=2)

        self.RIGHT_LOWER_FRAME_3 = tk.LabelFrame(self.RIGHT_FRAME_3, text="Resampled regression curves")
        self.RIGHT_LOWER_FRAME_3.pack(anchor=tk.NE, side=tk.TOP, fill=tk.BOTH, padx=2, pady=2)
        self.RESAMPLE_FIGURE_CANVAS_3 = tk.Canvas(self.RIGHT_LOWER_FRAME_3)
        self.RESAMPLE_FIGURE_CANVAS_3.pack(anchor=tk.N, side=tk.TOP, fill=tk.BOTH, expand=True, padx=2, pady=2)

        self.regression_fig_lims = None
        self.resampling_fig_lims = None
        self._refreshSplineCanvas(None)

    def _resampleRateChanged(self, _1, _2, _3):
        # if float(self.resampling_period_days.get()) <= 0:
        #     self.resampling_period_days.set(1.0)
        self._refreshResamplingCanvas()

    def _refreshRegressionCanvas(self):
        if self.regression_fig_canv:
            self.regression_fig_canv.get_tk_widget().destroy()
        self.CURVE_FIGURE_CANVAS_3.delete(tk.ALL)
        self._redrawRegressionFigure()
        if self.regression_fig_lims is not None:
            self.regression_fig.gca().set_xlim(self.regression_fig_lims[0])
            self.regression_fig.gca().set_ylim(self.regression_fig_lims[1])
        else:
            ax = self.regression_fig.gca()
            self.regression_fig_lims = (ax.get_xlim(), ax.get_ylim())
        self.regression_fig_canv = FigureCanvasTkAgg(self.regression_fig, master=self.CURVE_FIGURE_CANVAS_3)
        self.regression_fig_canv.get_tk_widget().pack(anchor=tk.N, side=tk.TOP, fill=tk.BOTH, expand=True)
        self.CURVE_FIGURE_CANVAS_3.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.update()

    def _redrawRegressionFigure(self):
        # Harvest settings from widgets, then run regression
        self.degrees = {}
        self.selected_cluster_ids = copy.copy(self.cluster_ids)
        for i, c in enumerate(self.DEGREE_CONTROLS):
            if c.is_enabled:
                self.degrees[i] = c.spline_degree
            else:
                self.selected_cluster_ids[self.selected_cluster_ids == i] = -1
        self.regression_curves = self.xcorr_ca.do_spline_regression(self.selected_cluster_ids, self.degrees)

        if self.regression_fig is not None:
            self.regression_fig.clear()
        else:
            self.regression_fig = plt.figure(figsize=(16, 9), dpi=self.display_dpi)
            # self.regression_fig = plt.figure()

        self.xcorr_ca.plot_regressors(self.regression_fig.gca(), self.selected_cluster_ids, self.regression_curves,
                                      self.station_code)
        self.regression_fig.tight_layout()
        self.regression_fig.autofmt_xdate()
        # TODO: Add auto-saving to file

    def _refreshResamplingCanvas(self):
        if self.resampling_fig_canv:
            self.resampling_fig_canv.get_tk_widget().destroy()
        self.RESAMPLE_FIGURE_CANVAS_3.delete(tk.ALL)
        self._redrawResamplingFigure()
        if self.resampling_fig_lims is not None:
            self.resampling_fig.gca().set_xlim(self.resampling_fig_lims[0])
            self.resampling_fig.gca().set_ylim(self.resampling_fig_lims[1])
        else:
            ax = self.resampling_fig.gca()
            self.resampling_fig_lims = (ax.get_xlim(), ax.get_ylim())
        self.resampling_fig_canv = FigureCanvasTkAgg(self.resampling_fig, master=self.RESAMPLE_FIGURE_CANVAS_3)
        self.resampling_fig_canv.get_tk_widget().pack(anchor=tk.N, side=tk.TOP, fill=tk.BOTH, expand=True)
        self.RESAMPLE_FIGURE_CANVAS_3.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.update()

    def _redrawResamplingFigure(self):
        # Harvest settings from widgets, then run resampling
        sec_per_day = 24*3600.0
        sample_period = self.resampling_period_days.get()*sec_per_day
        self.regular_corrections = self.xcorr_ca.do_spline_resampling(self.selected_cluster_ids, self.regression_curves,
                                                                      sample_period)

        if self.resampling_fig is not None:
            self.resampling_fig.clear()
        else:
            self.resampling_fig = plt.figure(figsize=(16, 9), dpi=self.display_dpi)
            # self.resampling_fig = plt.figure()

        self.xcorr_ca.plot_resampled_clusters(self.resampling_fig.gca(), self.selected_cluster_ids,
                                              self.regular_corrections, self.station_code)
        self.resampling_fig.tight_layout()
        self.resampling_fig.autofmt_xdate()
        # TODO: Add auto-saving to file

    def _refreshSplineCanvas(self, _new_val_unused):
        self._refreshRegressionCanvas()
        self._refreshResamplingCanvas()

    def _exportCorrections(self):
        # Save plots to PNG files.
        try:
            self.busy()
            regression_fig_file = os.path.join(self.output_folder.get(), self.station_code + "_regression_profile.png")
            self.regression_fig.savefig(regression_fig_file, dpi=300)
            resampling_fig_file = os.path.join(self.output_folder.get(),
                                               self.station_code + "_clock_correction_profile.png")
            self.resampling_fig.savefig(resampling_fig_file, dpi=300)

            # Save resampled data to csv.
            data_blocks = []
            for c in self.regular_corrections.values():
                # BEWARE: The 'corrections' array sign is negated there, since the correction
                # we have computed up to this point is actually the clock *error*. Subtraction
                # of an error is the same as addition of a correction of opposite sign.
                data_blocks.append(pd.DataFrame(np.column_stack([c['times'], -c['corrections']]),
                                                columns=['timestamp', 'clock_correction']))
            df = pd.concat(data_blocks)

            df['date'] = df['timestamp'].apply(obspy.UTCDateTime).apply(lambda x: x.date)
            net, sta = self.station_code.split('.')
            df['net'] = net
            df['sta'] = sta
            df = df[['net', 'sta', 'date', 'clock_correction']]

            output_file = os.path.join(self.output_folder.get(), self.station_code + "_clock_correction.csv")
            df.to_csv(output_file, index=False)

            self.not_busy()
            messagebox.showinfo("Export success", "Saved following files:\n{}\n{}\n{}".format(
                regression_fig_file, resampling_fig_file, output_file), width=100)
        except Exception as e:
            self.not_busy()
            messagebox.showerror("Export error", "Error occurred during exported - results might not be saved!\n"
                                                 "Error:\n{}".format(str(e)))
            return

        self.ROOT_FRAME_3.destroy()
        self.ROOT_FRAME_3 = None
        self._quitApp()

#end class


class SplineDegreeWidget(tk.LabelFrame):
    def __init__(self, index, initial_value, master=None, command=None):
        assert initial_value >= 1 and initial_value <= 5
        tk.LabelFrame.__init__(self, master, text="Spline {}".format(index))
        self.enabled = tk.BooleanVar(self, True)
        self.degree = tk.IntVar(self, initial_value)
        self.ENABLE_TOGGLE = tk.Checkbutton(self, variable=self.enabled)
        self.ENABLE_TOGGLE['command'] = self._toggled
        self.ENABLE_TOGGLE.pack(anchor=tk.E, side=tk.LEFT)
        self.DEGREE_LABEL = tk.Label(self, text="Spline order:")
        self.DEGREE_LABEL.pack(anchor=tk.E, side=tk.LEFT, padx=2)
        self.DEGREE_CHOOSER = tk.OptionMenu(self, self.degree, 1, 2, 3, 4, 5, command=command)
        self.DEGREE_CHOOSER.pack(anchor=tk.E, side=tk.LEFT, padx=2, fill=tk.X, expand=1)
        self.command_handler = command

    def _toggled(self):
        if self.enabled.get():
            self.DEGREE_LABEL['state'] = tk.NORMAL
            self.DEGREE_CHOOSER['state'] = tk.NORMAL
            if self.command_handler is not None:
                self.command_handler(True)
        else:
            self.DEGREE_LABEL['state'] = tk.DISABLED
            self.DEGREE_CHOOSER['state'] = tk.DISABLED
            if self.command_handler is not None:
                self.command_handler(False)

    @property
    def is_enabled(self):
        return self.enabled.get()

    @property
    def spline_degree(self):
        return self.degree.get()


tk_root = tk.Tk()
app = GpsClockCorrectionApp(master=tk_root)
app.master.title("GPS Clock Correction Workflow")
app.master.minsize(320, 120)
app.master.columnconfigure(0, weight=1)
app.master.rowconfigure(0, weight=1)
app.mainloop()
# tk_root.destroy()
