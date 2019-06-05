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


class GpsClockCorrectionApp(tk.Frame):

    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        self.nc_file = tk.StringVar()
        self.pack()
        self._last_dir = None
        self._createInitialWidgets()

    def _createInitialWidgets(self):
        self.NC_FILE_LABEL = tk.Label(self, text="Cross-correlation file:")
        self.NC_FILE_LABEL.pack(side=tk.LEFT)

        self.NC_FILE = tk.Entry(self)
        self.NC_FILE['width'] = 64
        self.NC_FILE.pack(padx=5, side=tk.LEFT)
        self.NC_FILE['textvariable'] = self.nc_file

        self.OPEN = tk.Button(self)
        self.OPEN['text'] = "Open..."
        self.OPEN['command'] = self.openNcFile
        self.OPEN.pack(side=tk.LEFT)

        self.NEXT = tk.Button(self)
        self.NEXT['text'] = "Next..."
        self.NEXT['state'] = tk.DISABLED
        self.NEXT.pack(side=tk.LEFT)
        self.nc_file.trace_variable('w', self.updateNextButton)

        self.QUIT = tk.Button(self)
        self.QUIT['text'] = "Quit"
        self.QUIT['command'] = self.quit
        self.QUIT.pack(side=tk.RIGHT)

    def updateNextButton(self, _1, _2, _3):
        nc_file_value = self.nc_file.get()
        nc_file_valid = bool(nc_file_value) and os.path.exists(nc_file_value)
        self.NEXT['state'] = tk.NORMAL if nc_file_valid else tk.DISABLED

    def openNcFile(self):
        initial_dir = '/g/data/ha3/Passive/SHARED_DATA/GPS_Clock/xcorr' if self._last_dir is None else self._last_dir
        file_name = tkFileDialog.askopenfilename(initialdir=initial_dir,
                                                 title='Select .nc file to analyze',
                                                 filetypes=(("nc files", "*.nc"), ("all files", "*")))
        self._last_dir = os.path.split(file_name)[0]
        self.nc_file.set(file_name)


tk_root = tk.Tk()
app = GpsClockCorrectionApp(master=tk_root)
app.master.title("GPS Clock Correction Workflow")
app.master.minsize(320, 160)
app.master.columnconfigure(0, weight=1)
app.master.rowconfigure(0, weight=1)
app.mainloop()
# tk_root.destroy()
