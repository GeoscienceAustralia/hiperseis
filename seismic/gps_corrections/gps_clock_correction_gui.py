#!/bin/env python
"""
GUI interface to generating clock corrections from x-corr results.
"""

import os
import copy
import stat

# pylint: disable=invalid-name, unresolved-import

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


class GpsClockCorrectionApp(tk.Frame):  # pragma: no cover

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

        self._create_step0_widgets()

    def busy(self):
        tk_root.config(cursor="watch")
        tk_root.update()

    def not_busy(self):
        tk_root.config(cursor="")

    def _create_step0_widgets(self):
        self.ROOT_FRAME_0 = tk.Frame(self)
        self.ROOT_FRAME_0.pack(fill=tk.BOTH, expand=1)

        upper_frame_0 = tk.LabelFrame(self.ROOT_FRAME_0, text="File selection", borderwidth=2)
        upper_frame_0.pack(anchor=tk.N, side=tk.TOP, fill=tk.X, padx=2, pady=2)

        self.stn_code_label = tk.Label(self.ROOT_FRAME_0, text="Station code: ", font=8)
        self.stn_code_label.pack(side=tk.TOP, fill=tk.X, padx=2, pady=8)

        lower_frame_0 = tk.Frame(self.ROOT_FRAME_0)
        lower_frame_0.pack(anchor=tk.S, side=tk.BOTTOM, fill=tk.X, padx=2, pady=2)

        file_entry_frame = tk.Frame(upper_frame_0)
        file_entry_frame.pack(anchor=tk.E, side=tk.TOP, fill=tk.X)
        open_button = tk.Button(file_entry_frame, width=10, text="Open...", command=self._open_nc_file)
        open_button.pack(anchor=tk.W, side=tk.LEFT, padx=2, pady=2)

        nc_file_label = tk.Label(file_entry_frame, text="Cross-correlation file:")
        nc_file_label.pack(anchor=tk.E, side=tk.LEFT, padx=(8,0), pady=2)

        nc_file_entry = tk.Entry(file_entry_frame, exportselection=False, width=64)
        nc_file_entry.pack(anchor=tk.E, padx=(2,4), pady=2, side=tk.TOP, fill=tk.X)
        nc_file_entry['textvariable'] = self.nc_file

        output_folder_frame = tk.Frame(upper_frame_0)
        output_folder_frame.pack(anchor=tk.E, side=tk.TOP, fill=tk.X)
        select_button = tk.Button(output_folder_frame, width=10, text="Choose...", command=self._choose_output_folder)
        select_button.pack(anchor=tk.W, side=tk.LEFT, padx=2, pady=2)

        output_folder_label = tk.Label(output_folder_frame, text="Output folder:")
        output_folder_label.pack(anchor=tk.E, side=tk.LEFT, padx=(8,0), pady=2)

        output_folder_entry = tk.Entry(output_folder_frame, exportselection=False, width=64)
        output_folder_entry.pack(anchor=tk.E, padx=(2,4), pady=2, side=tk.TOP, fill=tk.X)
        output_folder_entry['textvariable'] = self.output_folder

        self.next_button = tk.Button(lower_frame_0, text="Continue...", state=tk.DISABLED, command=self._goto_step1)
        self.next_button.pack(anchor=tk.SW, side=tk.LEFT)
        self.nc_file.trace_variable('w', self._update_next_button)
        self.output_folder.trace_variable('w', self._update_next_button)

        self.quit_button = tk.Button(lower_frame_0, text="Quit", command=self._quit_app)
        self.quit_button.pack(anchor=tk.SE, side=tk.RIGHT)

    def _quit_app(self):
        self._destroy_figures()
        self.quit()

    def _extract_code_from_filename(self, src_file):
        _, basename = os.path.split(src_file)
        _, file_type = os.path.splitext(src_file)
        name_parts = basename.split('.')
        netcode = name_parts[0]
        statcode = name_parts[1]
        full_code = '.'.join([netcode, statcode])
        return full_code

    def _update_next_button(self, *_unused):
        nc_file_value = self.nc_file.get()
        nc_file_valid = bool(nc_file_value) and os.path.isfile(nc_file_value)
        if nc_file_valid:
            self.station_code = self._extract_code_from_filename(nc_file_value)
        else:
            self.station_code = ''
        self.stn_code_label['text'] = "Station code: " + self.station_code
        output_folder_value = self.output_folder.get()
        output_folder_valid = bool(output_folder_value) and os.path.isdir(output_folder_value)
        self.next_button['state'] = tk.NORMAL if (nc_file_valid and output_folder_valid) else tk.DISABLED

    def _open_nc_file(self):
        initial_dir = self.nc_file.get()
        if not initial_dir:
            initial_dir = '/g/data/ha3/Passive/SHARED_DATA/GPS_Clock/xcorr'
        file_name = filedialog.askopenfilename(initialdir=initial_dir,
                                               title='Select .nc file to analyze',
                                               filetypes=(("nc files", "*.nc"), ("all files", "*")))
        if file_name:
            self.nc_file.set(file_name)

    def _choose_output_folder(self):
        initial_dir = self.output_folder.get()
        if not initial_dir:
            initial_dir = '/g/data/ha3/Passive/SHARED_DATA/GPS_Clock/corrections'
        folder = filedialog.askdirectory(initialdir=initial_dir, title='Select output folder', mustexist=True)
        if folder:
            self.output_folder.set(folder)

    def _goto_step1(self):
        if self.current_step == 0:
            # for child in self.TOP_PANE_0.winfo_children():
            #     self.child.destroy()
            self.ROOT_FRAME_0.destroy()
            self.ROOT_FRAME_0 = None
            self._create_step1_widgets()

    def _create_step1_widgets(self):
        self.busy()
        self.current_step = 1

        self.ROOT_FRAME_1 = tk.Frame(self)
        self.ROOT_FRAME_1.pack(fill=tk.BOTH, expand=1)

        left_frame_1 = tk.LabelFrame(self.ROOT_FRAME_1, width=800)
        left_frame_1.pack(anchor=tk.NW, side=tk.LEFT, fill=tk.X, padx=2, pady=2)

        self.stn_code_label = tk.Label(left_frame_1, text="Station code: " + self.station_code, font=8)
        self.stn_code_label.pack(anchor=tk.W, side=tk.TOP)

        time_window_frame = tk.LabelFrame(left_frame_1, borderwidth=2)
        time_window_frame.pack(anchor=tk.W, side=tk.TOP, fill=tk.X, padx=2, pady=2)

        time_window_label = tk.Label(time_window_frame, text="Time window (s)")
        time_window_label.pack(anchor=tk.NW, side=tk.TOP)
        time_window_entry = tk.Scale(time_window_frame, from_=30, to=1800, resolution=10,
                                     orient=tk.HORIZONTAL, length=500)
        time_window_entry.set(self.time_window.get())
        time_window_entry['variable'] = self.time_window
        time_window_entry['command'] = self._enable_refresh
        time_window_entry.pack(anchor=tk.SW, side=tk.BOTTOM, fill=tk.X)

        snr_frame = tk.LabelFrame(left_frame_1, borderwidth=2)
        snr_frame.pack(anchor=tk.W, side=tk.TOP, fill=tk.X, padx=2, pady=2)

        snr_label = tk.Label(snr_frame, text="SNR")
        snr_label.pack(anchor=tk.NW, side=tk.TOP)
        snr_entry = tk.Scale(snr_frame, from_=0, to=100, resolution=1, orient=tk.HORIZONTAL, length=500)
        snr_entry.set(self.snr_threshold.get())
        snr_entry['variable'] = self.snr_threshold
        snr_entry['command'] = self._enable_refresh
        snr_entry.pack(anchor=tk.SW, side=tk.BOTTOM, fill=tk.X)

        self.right_figure_canvas_1 = tk.Canvas(self.ROOT_FRAME_1)
        self.right_figure_canvas_1.pack(anchor=tk.NE, side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=4, pady=4)

        self.refresh_button = tk.Button(left_frame_1, text="Refresh", state=tk.DISABLED,
                                        command=self._update_step1_canvas)
        self.refresh_button.pack(anchor=tk.NW, side=tk.TOP, fill=tk.X, padx=2, pady=2)

        self.next_button = tk.Button(left_frame_1, text="Continue...", command=self._goto_step2)
        self.next_button.pack(anchor=tk.W, side=tk.LEFT, padx=2, pady=(16, 2))

        self.quit_button = tk.Button(left_frame_1, text="Quit", command=self._quit_app)
        self.quit_button.pack(anchor=tk.E, side=tk.RIGHT, padx=2, pady=(16, 2))

        self.xcorr_settings, self.xcorr_title_tag = read_correlator_config(self.nc_file.get())
        self.not_busy()

        self.next_button['state'] = tk.DISABLED
        self.next_button.update()
        self.quit_button['state'] = tk.DISABLED
        self.quit_button.update()
        self._update_step1_canvas()
        self.next_button['state'] = tk.NORMAL
        self.quit_button['state'] = tk.NORMAL

    def _enable_refresh(self, _new_val):
        self.refresh_button['state'] = tk.NORMAL
        self.next_button['state'] = tk.DISABLED

    def _destroy_figures(self):
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

    def _update_step1_canvas(self):
        self.busy()
        self.refresh_button['state'] = tk.DISABLED
        self.next_button['state'] = tk.DISABLED
        self._destroy_figures()
        self.right_figure_canvas_1.delete(tk.ALL)
        self.xcorr_ca, self.xcorr_fig = \
            plot_xcorr_file_clock_analysis(self.nc_file.get(), self.fds, self.time_window.get(),
                                           self.snr_threshold.get(), self.pearson_cutoff_factor, show=False,
                                           title_tag=self.xcorr_title_tag, settings=self.xcorr_settings)
        self.fig_canv = FigureCanvasTkAgg(self.xcorr_fig, master=self.right_figure_canvas_1)
        self.fig_canv.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.right_figure_canvas_1.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.next_button['state'] = tk.NORMAL
        self.not_busy()
        self.update()

    def _goto_step2(self):
        if self.current_step == 1:
            self.next_button['state'] = tk.DISABLED

            self._destroy_figures()

            self.right_figure_canvas_1.delete(tk.ALL)
            self.right_figure_canvas_1 = None
            self.ROOT_FRAME_1.destroy()
            self.ROOT_FRAME_1 = None
            self.update()

            info_label = tk.Label(self, text="Saving xcorr plot...")
            info_label.pack(pady=30)
            self.update()

            # Generate PNG file for the .nc file using the current settings.
            self.busy()
            batch_process_xcorr([self.nc_file.get()], self.fds, self.time_window.get(), self.snr_threshold.get(),
                                pearson_cutoff_factor=self.pearson_cutoff_factor, save_plots=True, force_save=True)
            self.not_busy()

            info_label.destroy()

            self._create_step2_widgets()

    def _create_step2_widgets(self):
        self.current_step = 2

        assert self.xcorr_ca is not None

        self.ROOT_FRAME_2 = tk.Frame(self)
        self.ROOT_FRAME_2.pack(fill=tk.BOTH, expand=1)

        upper_frame_2 = tk.LabelFrame(self.ROOT_FRAME_2, text="Station code: " + self.station_code)
        upper_frame_2.pack(anchor=tk.N, side=tk.TOP)

        coeff_widgets_frame_2 = tk.Frame(upper_frame_2)
        coeff_widgets_frame_2.pack(anchor=tk.NW, side=tk.TOP, padx=2, pady=2)

        coeff_widgets_label = tk.Label(coeff_widgets_frame_2, font=8,
                                       text="Tune coefficients to optimize clustering and filtering of outliers")
        coeff_widgets_label.pack(anchor=tk.W, side=tk.TOP, pady=2)

        coeff0_label = tk.LabelFrame(coeff_widgets_frame_2, text="Coefficient 0 (x-separation)")
        coeff0_label.pack(anchor=tk.W, side=tk.LEFT, padx=2)
        coeff0_spinbox = tk.Spinbox(coeff0_label, from_=0.0, to=100.0, increment=0.1,
                                    textvariable=self.cluster_coeff0)
        coeff0_spinbox.pack(padx=2, pady=2)

        coeff1_label = tk.LabelFrame(coeff_widgets_frame_2, text="Coefficient 1 (y-separation)")
        coeff1_label.pack(anchor=tk.W, side=tk.LEFT, padx=2)
        coeff1_spinbox = tk.Spinbox(coeff1_label, from_=0.0, to=100.0, increment=0.1,
                                         textvariable=self.cluster_coeff1)
        coeff1_spinbox.pack(padx=2, pady=2)

        coeff2_label = tk.LabelFrame(coeff_widgets_frame_2, text="Coefficient 2 (slope)")
        coeff2_label.pack(anchor=tk.W, side=tk.LEFT, padx=2)
        coeff2_spinbox = tk.Spinbox(coeff2_label, from_=0.0, to=100.0, increment=0.1,
                                    textvariable=self.cluster_coeff2)
        coeff2_spinbox.pack(padx=2, pady=2)

        control_buttons_frame_2 = tk.Frame(upper_frame_2)
        control_buttons_frame_2.pack(anchor=tk.NW, side=tk.TOP, padx=2, pady=2, fill=tk.X)

        self.next_button = tk.Button(control_buttons_frame_2)
        self.next_button['text'] = "Continue..."
        self.next_button['command'] = self._goto_step3
        self.next_button.pack(anchor=tk.NW, side=tk.LEFT)

        self.quit_button = tk.Button(control_buttons_frame_2)
        self.quit_button['text'] = "Quit"
        self.quit_button['command'] = self._quit_app
        self.quit_button.pack(anchor=tk.SE, side=tk.RIGHT)

        lower_frame_2 = tk.LabelFrame(self.ROOT_FRAME_2, text="Clustering Result")
        lower_frame_2.pack(anchor=tk.N, side=tk.TOP)

        self.cluster_figure_canvas_2 = tk.Canvas(lower_frame_2)
        self.cluster_figure_canvas_2.pack(anchor=tk.N, side=tk.TOP, fill=tk.BOTH, expand=True, padx=4, pady=4)

        self._refresh_clustering_canvas()

        self.update()

        coeff0_spinbox['command'] = self._refresh_clustering_canvas
        coeff1_spinbox['command'] = self._refresh_clustering_canvas
        coeff2_spinbox['command'] = self._refresh_clustering_canvas
        coeff0_spinbox.bind('<FocusOut>', self._refresh_clustering_canvas)
        coeff0_spinbox.bind('<Return>', self._refresh_clustering_canvas)
        coeff1_spinbox.bind('<FocusOut>', self._refresh_clustering_canvas)
        coeff1_spinbox.bind('<Return>', self._refresh_clustering_canvas)
        coeff2_spinbox.bind('<FocusOut>', self._refresh_clustering_canvas)
        coeff2_spinbox.bind('<Return>', self._refresh_clustering_canvas)

    def _refresh_clustering_canvas(self, *_unused):
        if self.cluster_fig_canv:
            self.cluster_fig_canv.get_tk_widget().destroy()
        self.cluster_figure_canvas_2.delete(tk.ALL)
        self._redraw_clustering_figure()
        self.cluster_fig_canv = FigureCanvasTkAgg(self.cluster_fig, master=self.cluster_figure_canvas_2)
        self.cluster_fig_canv.get_tk_widget().pack(anchor=tk.N, side=tk.TOP, fill=tk.BOTH, expand=True)
        self.cluster_figure_canvas_2.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.update()

    def _redraw_clustering_figure(self):
        cluster_coeffs = (self.cluster_coeff0.get(), self.cluster_coeff1.get(), self.cluster_coeff2.get())
        _, self.cluster_ids = self.xcorr_ca.do_clustering(cluster_coeffs)

        if self.cluster_fig is not None:
            self.cluster_fig.clear()
        else:
            self.cluster_fig = plt.figure(figsize=(16,9))

        self.xcorr_ca.plot_clusters(self.cluster_fig.gca(), self.cluster_ids, cluster_coeffs, self.station_code)
        self.cluster_fig.tight_layout()
        self.cluster_fig.autofmt_xdate()

    def _goto_step3(self):
        if self.current_step == 2:
            self.next_button['state'] = tk.DISABLED

            self.cluster_fig_file = os.path.join(self.output_folder.get(), self.station_code + "_1_clustering_profile.png")
            self.cluster_fig.savefig(self.cluster_fig_file, dpi=300)
            # Change permissions so that others in same group can read and write over this file.
            os.chmod(self.cluster_fig_file, (stat.S_IRGRP | stat.S_IWGRP))

            assert self.cluster_ids is not None
            num_clusters = len(set(self.cluster_ids[self.cluster_ids != -1]))
            self.degrees = dict(zip(range(num_clusters), [1]*num_clusters))

            self._destroy_figures()

            self.cluster_figure_canvas_2.delete(tk.ALL)
            self.cluster_figure_canvas_2 = None
            self.ROOT_FRAME_2.destroy()
            self.ROOT_FRAME_2 = None
            self.update()

            self._create_step3_widgets()

    def _create_step3_widgets(self):
        self.current_step = 3

        assert self.degrees is not None

        self.ROOT_FRAME_3 = tk.Frame(self)
        self.ROOT_FRAME_3.pack(fill=tk.BOTH, expand=1)

        left_frame_3 = tk.LabelFrame(self.ROOT_FRAME_3, text="Station code: " + self.station_code)
        left_frame_3.pack(anchor=tk.NW, side=tk.LEFT, fill=tk.X, padx=2, pady=2)

        self.degree_controls = []
        for i, d in self.degrees.items():
            dc = SplineDegreeWidget(i, d, left_frame_3, self._refresh_spline_canvas)
            dc.pack(anchor=tk.NW, side=tk.TOP, fill=tk.X, padx=4)
            self.degree_controls.append(dc)

        resampling_period_label = tk.Label(left_frame_3, text="Enter resampling period (days):")
        resampling_period_label.pack(anchor=tk.W, padx=2, pady=2, side=tk.TOP)

        resampling_period_entry = tk.Entry(left_frame_3, textvariable=self.resampling_period_days,
                                           exportselection=False)
        resampling_period_entry['width'] = 8
        resampling_period_entry.pack(anchor=tk.E, padx=5, pady=2, side=tk.TOP)
        self.resampling_period_days.trace_variable('w', self._resample_rate_changed)

        self.export_button = tk.Button(left_frame_3, text="Export...", command=self._export_corrections)
        self.export_button.pack(anchor=tk.SW, side=tk.LEFT, padx=2, pady=(16, 2))

        self.quit_button = tk.Button(left_frame_3, text="Quit", command=self._quit_app)
        self.quit_button.pack(anchor=tk.SE, side=tk.RIGHT, padx=2, pady=(16, 2))

        right_frame_3 = tk.Frame(self.ROOT_FRAME_3)
        right_frame_3.pack(anchor=tk.NE, side=tk.RIGHT, fill=tk.BOTH)

        right_upper_frame_3 = tk.LabelFrame(right_frame_3, text="Regression curves overlay")
        right_upper_frame_3.pack(anchor=tk.NE, side=tk.TOP, fill=tk.BOTH, padx=2, pady=2)
        self.curve_figure_canvas_3 = tk.Canvas(right_upper_frame_3)
        self.curve_figure_canvas_3.pack(anchor=tk.N, side=tk.TOP, fill=tk.BOTH, expand=True, padx=2, pady=2)

        right_lower_frame_3 = tk.LabelFrame(right_frame_3, text="Resampled regression curves")
        right_lower_frame_3.pack(anchor=tk.NE, side=tk.TOP, fill=tk.BOTH, padx=2, pady=2)
        self.resample_figure_canvas_3 = tk.Canvas(right_lower_frame_3)
        self.resample_figure_canvas_3.pack(anchor=tk.N, side=tk.TOP, fill=tk.BOTH, expand=True, padx=2, pady=2)

        self.regression_fig_lims = None
        self.resampling_fig_lims = None
        self._refresh_spline_canvas(None)

    def _resample_rate_changed(self, *_unused):
        # if float(self.resampling_period_days.get()) <= 0:
        #     self.resampling_period_days.set(1.0)
        self._refresh_resampling_canvas()

    def _update_export_button_state(self):
        export_enabled = False
        for c in self.degree_controls:
            export_enabled = export_enabled or c.is_enabled
        self.export_button['state'] = tk.NORMAL if export_enabled else tk.DISABLED

    def _refresh_spline_canvas(self, _new_val_unused):
        self._update_export_button_state()
        try:
            self._refresh_regression_canvas()
        except ValueError:
            pass
        try:
            self._refresh_resampling_canvas()
        except ValueError:
            pass

    def _refresh_regression_canvas(self):
        if self.regression_fig_canv:
            self.regression_fig_canv.get_tk_widget().destroy()
        self.curve_figure_canvas_3.delete(tk.ALL)
        self._redraw_regression_figure()
        if self.regression_fig_lims is not None:
            self.regression_fig.gca().set_xlim(self.regression_fig_lims[0])
            self.regression_fig.gca().set_ylim(self.regression_fig_lims[1])
        else:
            ax = self.regression_fig.gca()
            self.regression_fig_lims = (ax.get_xlim(), ax.get_ylim())
        self.regression_fig_canv = FigureCanvasTkAgg(self.regression_fig, master=self.curve_figure_canvas_3)
        self.regression_fig_canv.get_tk_widget().pack(anchor=tk.N, side=tk.TOP, fill=tk.BOTH, expand=True)
        self.curve_figure_canvas_3.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.update()

    def _redraw_regression_figure(self):
        # Harvest settings from widgets, then run regression
        self.degrees = {}
        self.selected_cluster_ids = copy.copy(self.cluster_ids)
        for i, c in enumerate(self.degree_controls):
            if c.is_enabled:
                self.degrees[i] = c.spline_degree
            else:
                self.selected_cluster_ids[self.selected_cluster_ids == i] = -1
        self.regression_curves = self.xcorr_ca.do_spline_regression(self.selected_cluster_ids, self.degrees)

        if self.regression_fig is not None:
            self.regression_fig.clear()
        else:
            self.regression_fig = plt.figure(figsize=(16, 9), dpi=self.display_dpi)

        if self.regression_curves:
            self.xcorr_ca.plot_regressors(self.regression_fig.gca(), self.selected_cluster_ids, self.regression_curves,
                                          self.station_code)
            self.regression_fig.tight_layout()
            self.regression_fig.autofmt_xdate()

    def _refresh_resampling_canvas(self):
        if self.resampling_fig_canv:
            self.resampling_fig_canv.get_tk_widget().destroy()
        self.resample_figure_canvas_3.delete(tk.ALL)
        self._redraw_resampling_figure()
        if self.resampling_fig_lims is not None:
            self.resampling_fig.gca().set_xlim(self.resampling_fig_lims[0])
            self.resampling_fig.gca().set_ylim(self.resampling_fig_lims[1])
        else:
            ax = self.resampling_fig.gca()
            self.resampling_fig_lims = (ax.get_xlim(), ax.get_ylim())
        self.resampling_fig_canv = FigureCanvasTkAgg(self.resampling_fig, master=self.resample_figure_canvas_3)
        self.resampling_fig_canv.get_tk_widget().pack(anchor=tk.N, side=tk.TOP, fill=tk.BOTH, expand=True)
        self.resample_figure_canvas_3.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.update()

    def _redraw_resampling_figure(self):
        # Harvest settings from widgets, then run resampling
        sec_per_day = 24*3600.0
        sample_period = self.resampling_period_days.get()*sec_per_day
        self.regular_corrections = self.xcorr_ca.do_spline_resampling(self.selected_cluster_ids, self.regression_curves,
                                                                      sample_period)

        if self.resampling_fig is not None:
            self.resampling_fig.clear()
        else:
            self.resampling_fig = plt.figure(figsize=(16, 9), dpi=self.display_dpi)

        if self.regular_corrections:
            self.xcorr_ca.plot_resampled_clusters(self.resampling_fig.gca(), self.selected_cluster_ids,
                                                  self.regular_corrections, self.station_code)
            self.resampling_fig.tight_layout()
            self.resampling_fig.autofmt_xdate()

    def _export_corrections(self):
        # Save plots to PNG files.
        try:
            self.export_button['state'] = tk.DISABLED
            self.update()
            self.busy()
            regression_fig_file = os.path.join(self.output_folder.get(), self.station_code + "_2_regression_profile.png")
            self.regression_fig.savefig(regression_fig_file, dpi=300)
            # Change permissions so that others in same group can read and write over this file.
            os.chmod(regression_fig_file, (stat.S_IRGRP | stat.S_IWGRP))
            resampling_fig_file = os.path.join(self.output_folder.get(),
                                               self.station_code + "_3_clock_correction_profile.png")
            self.resampling_fig.savefig(resampling_fig_file, dpi=300)
            os.chmod(resampling_fig_file, (stat.S_IRGRP | stat.S_IWGRP))

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
            # Change permissions so that others in same group can read and write over this file.
            os.chmod(output_file, (stat.S_IRGRP | stat.S_IWGRP))

            self.not_busy()
            messagebox.showinfo("Export success", "Saved following files:\n{}\n{}\n{}\n{}".format(
                self.cluster_fig_file, regression_fig_file, resampling_fig_file, output_file))
        except Exception as e:
            self.not_busy()
            messagebox.showerror("Export error", "Error occurred during exported - results might not be saved!\n"
                                                 "Error:\n{}".format(str(e)))
            return

        self.curve_figure_canvas_3.delete(tk.ALL)
        self.curve_figure_canvas_3 = None
        self.resample_figure_canvas_3.delete(tk.ALL)
        self.resample_figure_canvas_3 = None
        self.ROOT_FRAME_3.destroy()
        self.ROOT_FRAME_3 = None
        self._quit_app()

# end class


class SplineDegreeWidget(tk.LabelFrame):  # pragma: no cover
    def __init__(self, index, initial_value, master=None, command=None):
        assert initial_value >= 1 and initial_value <= 5
        tk.LabelFrame.__init__(self, master, text="Spline {}".format(index))
        self._enabled = tk.BooleanVar(self, True)
        self._degree = tk.IntVar(self, initial_value)
        self.ENABLE_TOGGLE = tk.Checkbutton(self, variable=self._enabled)
        self.ENABLE_TOGGLE['command'] = self._toggled
        self.ENABLE_TOGGLE.pack(anchor=tk.E, side=tk.LEFT)
        self.DEGREE_LABEL = tk.Label(self, text="Spline order:")
        self.DEGREE_LABEL.pack(anchor=tk.E, side=tk.LEFT, padx=2)
        self.DEGREE_CHOOSER = tk.OptionMenu(self, self._degree, 1, 2, 3, 4, 5, command=command)
        self.DEGREE_CHOOSER.pack(anchor=tk.E, side=tk.LEFT, padx=2, fill=tk.X, expand=1)
        self.command_handler = command

    def _toggled(self):
        if self._enabled.get():
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
        return self._enabled.get()

    @property
    def spline_degree(self):
        return self._degree.get()


if __name__ == '__main__': # pragma: no cover
    tk_root = tk.Tk()
    app = GpsClockCorrectionApp(master=tk_root)
    app.master.title("GPS Clock Correction Workflow")
    app.master.minsize(320, 120)
    app.master.columnconfigure(0, weight=1)
    app.master.rowconfigure(0, weight=1)
    app.mainloop()
    # tk_root.destroy()
# end if
