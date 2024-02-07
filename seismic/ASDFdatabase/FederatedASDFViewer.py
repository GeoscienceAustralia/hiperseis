"""
Description:
    Implements a web browser-based, nominally interactive GUI for viewing data holdings
    in a collection of ASDF files.

References:

CreationDate:   12/12/18
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     08/09/23   RH
"""

import remi.gui as gui
from remi import start, App
import os, sys
import numpy as np
from obspy import UTCDateTime
import click
import uuid
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
import io
import time
import threading
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from obspy.core import Trace, Stream
from obspy.signal import PPSD
import cartopy.crs as ccrs
from scipy.stats import circmean as cmean
from collections import defaultdict
from shapely import geometry

class CustomPPSD(PPSD):
    def __init__(self, stats, skip_on_gaps=False,
                 db_bins=(-150, 50, 1.), ppsd_length=3600.0, overlap=0.5,
                 special_handling=None, period_smoothing_width_octaves=1.0,
                 period_step_octaves=0.125, period_limits=None,
                 **kwargs):
        
        # flat response
        metadata = paz={'sensitivity': 1.0,
                        'gain': 1.0,
                        'poles': [0 + 1j],
                        'zeros': [0 + 1j]}
        
        super(CustomPPSD, self).__init__(stats, metadata, skip_on_gaps=skip_on_gaps,
                 db_bins=db_bins, ppsd_length=ppsd_length, overlap=overlap,
                 special_handling=special_handling, 
                 period_smoothing_width_octaves=period_smoothing_width_octaves,
                 period_step_octaves=period_step_octaves, 
                 period_limits=period_limits,
                 **kwargs)
    # end func
    
    def add(self, stream, verbose=False):
        if isinstance(stream, Trace):
            stream = Stream([stream])
        # end if
        
        # normalize streams
        stream = stream.copy()
        for tr in stream:
            if(tr.stats.npts > 0):
                tr.data = tr.data / np.max(np.fabs(tr.data))
            # end if
        # end for
        
        super(CustomPPSD, self).add(stream, verbose=verbose)
    # end func
    
    def _plot_histogram(self, fig, draw=False, filename=None):
        """
        Reuse a previously created figure returned by `plot(show=False)`
        and plot the current histogram stack (pre-computed using
        :meth:`calculate_histogram()`) into the figure. If a filename is
        provided, the figure will be saved to a local file.
        Note that many aspects of the plot are statically set during the first
        :meth:`plot()` call, so this routine can only be used to update with
        data from a new stack.
        """
        import matplotlib.pyplot as plt
        ax = fig.axes[0]
        xlim = ax.get_xlim()
        if "quadmesh" in fig.ppsd:
            fig.ppsd.pop("quadmesh").remove()

        if fig.ppsd.cumulative:
            data = self.current_histogram_cumulative * 100.0
        else:
            # avoid divison with zero in case of empty stack
            data = (
                self.current_histogram * 100.0 /
                (self.current_histogram_count or 1))

        xedges = self.period_xedges
        if fig.ppsd.xaxis_frequency:
            xedges = 1.0 / xedges

        if "meshgrid" not in fig.ppsd:
            fig.ppsd.meshgrid = np.meshgrid(xedges, self.db_bin_edges)
        ppsd = ax.pcolormesh(
            fig.ppsd.meshgrid[0], fig.ppsd.meshgrid[1], data.T,
            cmap=fig.ppsd.cmap, zorder=-1)
        fig.ppsd.quadmesh = ppsd

        if "colorbar" not in fig.ppsd:
            cb = plt.colorbar(ppsd, ax=ax)
            cb.mappable.set_clim(*fig.ppsd.color_limits)
            cb.set_label(fig.ppsd.label)
            fig.ppsd.colorbar = cb

        if fig.ppsd.max_percentage is not None:
            ppsd.set_clim(*fig.ppsd.color_limits)

        if fig.ppsd.grid:
            if fig.ppsd.cmap.name == "jet":
                color = {"color": "0.7"}
            else:
                color = {}
            ax.grid(True, which="major", **color)
            ax.grid(True, which="minor", **color)

        ax.set_xlim(*xlim)

        if filename is not None:
            plt.savefig(filename)
        elif draw:
            with np.errstate(under="ignore"):
                plt.draw()
        return fig    
# end class

def printException(e: Exception):
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, e, fname, exc_tb.tb_lineno)
# end func

class FigureImage(gui.Image):
    def __init__(self, **kwargs):
        super(FigureImage, self).__init__("/%s/get_image_data?update_index=0" % id(self), **kwargs)
        self._buf = None
        #self._buflock = threading.Lock()

        self._fig = kwargs.pop('fig')
        self.redraw()
    # end func

    def redraw(self):
        canv = FigureCanvasAgg(self._fig)
        buf = io.BytesIO()
        canv.print_figure(buf, format='png')
        #with self._buflock:
        if self._buf is not None:
            self._buf.close()
        self._buf = buf

        i = int(time.time() * 1e6)
        self.attributes['src'] = "/%s/get_image_data?update_index=%d" % (id(self), i)

        super(FigureImage, self).redraw()
    # end func

    def get_image_data(self, update_index):
        #with self._buflock:
        if self._buf is None:
            return None
        self._buf.seek(0)
        data = self._buf.read()

        return [data, {'Content-type': 'image/png'}]
    # end func
# end class

ROW_WIDGET_WIDTH = 1000
ROW_WIDGET_HEIGHT = 500
MAP_WIDGET_WIDTH = 1000
MAP_WIDGET_HEIGHT = 500
MAP_WIDGET_PADDING = 700
PADDING_FACTOR = 1.1
TRACE_FIG_WIDTH = 7
TRACE_FIG_HEIGHT = 3.5
MAP_FIG_WIDTH = 7
MAP_FIG_HEIGHT = 3.5
DEFAULT_TRC_LENGTH = '600'
font = {'family' : 'normal',
        'size'   : 8}

matplotlib.rc('font', **font)
matplotlib.rc('xtick', labelsize=8) 
matplotlib.rc('ytick', labelsize=8) 

class DataViewer(App):
    def __int__(self, *args, **kwargs):
        super().__int__(*args, **kwargs)
    # end func

    def getNetworks(self):
        result = set([item.split('.')[0] for item in self.fds.unique_coordinates.keys()])
        return sorted(list(result))
    # end func

    def getStations(self, net):
        result = set()
        for item in self.fds.unique_coordinates.keys():
            nc, sc = item.split('.')
            if(nc == net): result.add(sc)
        # end for
        return sorted(list(result))
    # end func

    def getLocations(self, net, sta):
        rows = self.fds.get_stations(UTCDateTime('1900-01-01'), UTCDateTime('2100-01-01'),
                                     network=net, station=sta)
        locs = set()
        for row in rows: locs.add(row[2])

        return sorted(list(locs))
    # end func

    def getChannels(self, net, sta, loc):
        rows = self.fds.get_stations(UTCDateTime('1900-01-01'), UTCDateTime('2100-01-01'),
                                     network=net, station=sta, location=loc)
        comps = set()
        for row in rows: comps.add(row[3])

        return sorted(list(comps))
    # end func

    def mapWidget(self):
        def setMapImage(nc, minLon=None, minLat=None, maxLon=None, maxLat=None):
            try:
                fig = Figure(figsize=(MAP_FIG_WIDTH, MAP_FIG_HEIGHT))
                stations = np.array(self.getStations(nc))

                if(len(stations)):
                    nstart, nend = self.fds.get_global_time_range(nc)

                    lons = []
                    lats = []
                    scodes = []
                    for sc in stations:
                        netsta = '{}.{}'.format(nc, sc)
                        lon, lat = self.fds.unique_coordinates[netsta]

                        lons.append(lon)
                        lats.append(lat)
                    # end for
                    lons = np.array(lons)
                    lats = np.array(lats)
                    
                    boundsProvided = False
                    if(minLon is None and  minLat is None and maxLon is None and  maxLat is None):
                        minLon = np.min(lons)
                        maxLon = np.max(lons)
                        minLat = np.min(lats)
                        maxLat = np.max(lats)
                        minLon -= 0.5
                        maxLon += 0.5
                        minLat -= 0.5
                        maxLat += 0.5
                    else:
                        boundsProvided = True
                        polygon = [(minLon,minLat), (maxLon, minLat), (maxLon, maxLat), (minLon, maxLat)]
                        polygon = geometry.Polygon(polygon)
                        
                        imask = np.zeros(stations.shape, dtype='?')
                        for i in np.arange(len(stations)):
                            p = geometry.Point(lons[i], lats[i])
                            if (polygon.contains(p)): imask[i] = 1
                        # end for
                        
                        stations = stations[imask]
                        lons = lons[imask]
                        lats = lats[imask]
                    # end if
                    #print('Bounds: [{}, {}] - [{}, {}]'.format(minLon, minLat, maxLon, maxLat))
                    
                    clon = np.mean(np.array([minLon, maxLon]))
                    if(len(lons)>0):
                        clon = cmean(lons, high=180, low=-180)
                    # end if
                    
                    crs = ccrs.PlateCarree(central_longitude=clon)
                    left = 0.05
                    bottom = 0.05
                    width = 0.9
                    height = 0.8
                    ax = fig.add_axes([left, bottom, width, height], projection=crs)
                    # draw coastlines.
                    ax.coastlines('50m')
                    gl = ax.gridlines(crs=crs, draw_labels=True,
                                      dms=False,
                                      linewidth=1, color='gray',
                                      alpha=0.5, linestyle='--')

                    ax.xaxis.set_major_formatter(gl.xformatter)

                    gl.xlabel_style = {'size':6, 'color': 'k'}
                    gl.ylabel_style = {'size':6, 'color': 'k'}
                    
                    # plot stations
                    for i, sc in enumerate(stations):
                        lon, lat = lons[i], lats[i]

                        px, py = lon - clon, lat
                        pxl, pyl = lon - clon + 0.02, lat - 0.1
                        ax.scatter(px, py, 20, transform=crs, marker='v', c='r', edgecolor='none', zorder=1)
                        ax.annotate(sc, xy=(pxl, pyl), fontsize=7, zorder=2)
                    # end for
                    
                    # set extents only when bounds have been provided
                    if(boundsProvided):
                        ax.set_extent([minLon, maxLon, minLat, maxLat], crs=ccrs.PlateCarree())
                    # end if
                    
                    title = 'Availability: {} - {}'.format(nstart.strftime('%Y-%m-%d'),
                                                           nend.strftime('%Y-%m-%d'))
                    fig.suptitle(title)
                    fig.tight_layout()
                # end if

                mi = FigureImage(fig=fig)

                # update widgets
                self.wrapperContainer.children[key].children['rightContainer'].children['plot'] = mi
                self.wrapperContainer.children[key].children['leftContainer']. \
                    children['lonBoundsBox'].children['min'].set_value(minLon)
                self.wrapperContainer.children[key].children['leftContainer']. \
                    children['lonBoundsBox'].children['max'].set_value(maxLon)
                self.wrapperContainer.children[key].children['leftContainer']. \
                    children['latBoundsBox'].children['min'].set_value(minLat)
                self.wrapperContainer.children[key].children['leftContainer']. \
                    children['latBoundsBox'].children['max'].set_value(maxLat)

            except Exception as e:
                printException(e)
            # end try
        # end func

        def mapNetChanged(emitter, value=None):
            nc = self.wrapperContainer.children[key].children['leftContainer'].children['nBox'].children['net'].get_value()
            
            # update plot
            self.wrapperContainer.children[key].children['rightContainer'].children['plot'] = gui.Label('Loading..')
            t = threading.Thread(target=setMapImage,
                                 args=(nc,))
            t.start()
        # end func

        def boundsChanged(emitter, value):
            nc = self.wrapperContainer.children[key].children['leftContainer'].children['nBox'].children['net'].get_value()
            minLon = float(self.wrapperContainer.children[key].children['leftContainer'].children['lonBoundsBox'].children['min'].get_value())
            maxLon = float(self.wrapperContainer.children[key].children['leftContainer'].children['lonBoundsBox'].children['max'].get_value())
            minLat = float(self.wrapperContainer.children[key].children['leftContainer'].children['latBoundsBox'].children['min'].get_value())
            maxLat = float(self.wrapperContainer.children[key].children['leftContainer'].children['latBoundsBox'].children['max'].get_value())

            # update plot
            self.wrapperContainer.children[key].children['rightContainer'].children['plot'] = gui.Label('Loading..')
            t = threading.Thread(target=setMapImage,
                                 args=(nc, minLon, minLat, maxLon, maxLat,))
            t.start()
        # end func

        def writeCoordinates(emitter, value=None):
            def pathConfirmed(pe, pv):
                coords = defaultdict(list)
                for k, v in self.fds.unique_coordinates.items():
                    net, sta = k.split('.')

                    if(net == nc): coords[sta] = [v[0], v[1]]
                # end for
                
                try:
                    with open(pv, 'w') as fh:
                        fh.write('# net, sta, lon, lat\n')
                        for sta in sorted(coords.keys()):
                            line = '{},{},{:3.4f},{:3.4f}\n'.format(net, sta, coords[sta][0], coords[sta][1])
                            fh.write(line)
                        # end for
                    # end with
                except Exception as e:
                    print('Failed to write coordinates to {}, with error {}'.format(pv, e))
                # end try
            # end func
            
            nc = self.wrapperContainer.children[key].children['leftContainer'].children['nBox'].children['net'].get_value()

            ofn = os.path.join(os.getcwd(), '{}.txt'.format(nc))
            pathDialog = gui.InputDialog('Select File', 'Output file name: ',
                                         initial_value=ofn, width=500)
            pathDialog.confirm_value.do(pathConfirmed)
            pathDialog.show(self)
        # end func

        key = str(uuid.uuid4())
        container = gui.HBox(width=MAP_WIDGET_WIDTH, height=MAP_WIDGET_HEIGHT, style={'margin': '0px auto'})
        leftContainer = gui.VBox(width=MAP_WIDGET_WIDTH*0.25, height=MAP_WIDGET_HEIGHT*0.45,
                                 style={'border': '1px solid blue', 'margin': '5px'})
        rightContainer = gui.VBox(width=MAP_WIDGET_WIDTH*0.75, height=MAP_WIDGET_HEIGHT,
                                  style={'border': '1px solid blue', 'margin': '0px'})

        #############################################################
        # populate leftContainer
        #############################################################
        nBox = gui.HBox(width=MAP_WIDGET_WIDTH*0.1, height=MAP_WIDGET_HEIGHT*0.1, style={'margin': '0px auto'})
        nLabelBox = gui.HBox(width=MAP_WIDGET_WIDTH*0.1, height=MAP_WIDGET_HEIGHT*0.1, style={'margin': '0px auto'})
        boundsLabelBox = gui.HBox(width=MAP_WIDGET_WIDTH*0.2, height=MAP_WIDGET_HEIGHT*0.1, style={'margin': '0px auto'})
        lonBoundsBox = gui.HBox(width=MAP_WIDGET_WIDTH*0.2, height=MAP_WIDGET_HEIGHT*0.1, style={'margin': '0px auto'})
        latBoundsBox = gui.HBox(width=MAP_WIDGET_WIDTH*0.2, height=MAP_WIDGET_HEIGHT*0.1, style={'margin': '0px auto'})

        # network selection
        net = gui.DropDown.new_from_list(self.getNetworks())
        net.onchange.do(mapNetChanged)
        nBox.append({'net':net})
        nLabelBox.append(gui.Label('Network'))

        # map-controls
        minLon = gui.Input(width=MAP_WIDGET_WIDTH*0.05)
        maxLon = gui.Input(width=MAP_WIDGET_WIDTH*0.05)
        minLat = gui.Input(width=MAP_WIDGET_WIDTH*0.05)
        maxLat = gui.Input(width=MAP_WIDGET_WIDTH*0.05)

        minLon.onchange.do(boundsChanged)
        maxLon.onchange.do(boundsChanged)
        minLat.onchange.do(boundsChanged)
        maxLat.onchange.do(boundsChanged)

        resetBounds = gui.Button('Reset Bounds', height=30, margin='1px auto')
        resetBounds.onclick.do(mapNetChanged)
        exportCoordinates = gui.Button('Export Coordinates', height=30, margin='1px auto')
        exportCoordinates.onclick.do(writeCoordinates)

        lonBoundsBox.append({'lonBoundsLabel': gui.Label('Lon: '),
                             'min':minLon, 'max': maxLon})
        latBoundsBox.append({'latBoundsLabel': gui.Label('Lat: '),
                             'min': minLat,
                             'max': maxLat})
        boundsLabelBox.append({'d1':gui.Label(''), 'minLabel': gui.Label('Min'),
                               'd2':gui.Label(''),'maxLabel': gui.Label('Max')})

        leftContainer.append({'nLabelBox': nLabelBox, 'nBox':nBox,
                              'boundsLabelBox':boundsLabelBox,
                              'lonBoundsBox': lonBoundsBox,
                              'latBoundsBox': latBoundsBox,
                              'resetBounds': resetBounds,
                              'exportCoordinates': exportCoordinates})

        #############################################################
        # populate rightContainer
        #############################################################
        rightContainer.append({'plot': gui.Label('Loading..')})

        t = threading.Thread(target=setMapImage,
                             args=(net.get_value(),))
        t.start()
        container.append({'leftContainer': leftContainer, 'rightContainer': rightContainer})

        return container, key
    # end func

    def rowWidget(self):
        # Define onclick callbacks
        def netChanged(emitter, value):
            nc = self.rowContainer.children[key].children['leftContainer'].children['nslcBox'].children['net'].get_value()
            sta = gui.DropDown.new_from_list(self.getStations(nc))
            sta.onchange.do(staChanged)
            self.rowContainer.children[key].children['leftContainer'].children['nslcBox'].children['sta'] = sta

            #print('netChanged')
            staChanged(None, None)
        # end func

        def staChanged(emitter, value):
            nc = self.rowContainer.children[key].children['leftContainer'].children['nslcBox'].children['net'].get_value()
            sc = self.rowContainer.children[key].children['leftContainer'].children['nslcBox'].children['sta'].get_value()

            loc = gui.DropDown.new_from_list(self.getLocations(nc, sc))
            loc.onchange.do(locChanged)
            self.rowContainer.children[key].children['leftContainer'].children['nslcBox'].children['loc'] = loc
            
            #print('staChanged')
            locChanged(None, None)
        # end func

        def locChanged(emitter, value):
            nc = self.rowContainer.children[key].children['leftContainer'].children['nslcBox'].children['net'].get_value()
            sc = self.rowContainer.children[key].children['leftContainer'].children['nslcBox'].children['sta'].get_value()
            lc = self.rowContainer.children[key].children['leftContainer'].children['nslcBox'].children['loc'].get_value()

            cha = gui.DropDown.new_from_list(self.getChannels(nc, sc, lc))
            cha.onchange.do(chaChanged)
            self.rowContainer.children[key].children['leftContainer'].children['nslcBox'].children['cha'] = cha

            #print('locChanged')
            chaChanged(None, None)
        # end func

        def getMeta(nc, sc, lc, cc):
            lon, lat = self.fds.unique_coordinates['{}.{}'.format(nc, sc)]
            st, et = self.fds.get_global_time_range(nc, sc, lc, cc)
            #print(nc, sc, lc, cc, st, et)
            locStr = "Lon: {:.2f}, Lat: {:.2f}".format(lon, lat)
            availStr = "Availability: {} - {}".format(st.strftime('%Y-%m-%d'), et.strftime('%Y-%m-%d'))

            return locStr, availStr
        # end func

        def setTraceImage(nc, sc, lc, cc, st=None, et=None):
            try:
                if(st is None and et is None):
                    st, et = self.fds.get_global_time_range(nc, sc, lc, cc)

                    # set start-time 
                    self.rowContainer.children[key].children['rightContainer'].children['startStepBox']. \
                        children['start'].set_value(st.strftime('%Y-%m-%dT%H:%M:%S'))
                # end if
                step = int(self.rowContainer.children[key].children['rightContainer'].children['startStepBox']. \
                    children['step'].get_value())
                isPPSD = int(self.rowContainer.children[key].children['rightContainer'].children['startStepBox']. \
                    children['ppsd'].get_value())

                fig = Figure(figsize=(TRACE_FIG_WIDTH, TRACE_FIG_HEIGHT))
                stream = self.fds.get_waveforms(nc, sc, lc, cc, st, st + step)
                
                if(len(stream)):                     
                    if(not isPPSD):
                        fig = stream.plot(fig=fig, handle=True, type='relative')
                    else:
                        ppsd = CustomPPSD(stream[0].stats)

                        ppsd.add(stream)
                        fig = ppsd.plot(show_percentiles=False, 
                                        show_coverage=False,
                                        show_noise_models=False, show=False)
                        fig.set_size_inches(TRACE_FIG_WIDTH, TRACE_FIG_HEIGHT)
                    # end if
                    fig.axes[0].text(0.01, 0.01,
                                     'SR: {} Hz'.format(stream[0].stats.sampling_rate),
                                     bbox=dict(facecolor='white', linewidth=0, alpha=0.7),
                                     color='k',
                                     fontsize=7, weight='bold',
                                     transform=fig.axes[0].transAxes)
                # end if

                ti = FigureImage(fig=fig)
                self.rowContainer.children[key].children['rightContainer'].children['plot'] = ti
            except Exception as e:
                printException(e)
            # end try
        # end func

        def chaChanged(emitter, value):
            nc = self.rowContainer.children[key].children['leftContainer'].children['nslcBox'].children['net'].get_value()
            sc = self.rowContainer.children[key].children['leftContainer'].children['nslcBox'].children['sta'].get_value()
            lc = self.rowContainer.children[key].children['leftContainer'].children['nslcBox'].children['loc'].get_value()
            cc = self.rowContainer.children[key].children['leftContainer'].children['nslcBox'].children['cha'].get_value()

            # update metadata
            locStr, availStr = getMeta(nc, sc, lc, cc)
            self.rowContainer.children[key].children['leftContainer'].\
                children['locLabel'].set_text(locStr)
            self.rowContainer.children[key].children['leftContainer'].\
                children['availLabel'].set_text(availStr)

            # update plot
            self.rowContainer.children[key].children['rightContainer'].children['plot'] = gui.Label('Loading..')
            t = threading.Thread(target=setTraceImage,
                                 args=(nc, sc, lc, cc))
            t.start()
            #print('chaChanged')
        # end func

        def startStepChanged(emitter, value=None):
            nc = self.rowContainer.children[key].children['leftContainer'].children['nslcBox'].children[
                'net'].get_value()
            sc = self.rowContainer.children[key].children['leftContainer'].children['nslcBox'].children[
                'sta'].get_value()
            lc = self.rowContainer.children[key].children['leftContainer'].children['nslcBox'].children[
                'loc'].get_value()
            cc = self.rowContainer.children[key].children['leftContainer'].children['nslcBox'].children[
                'cha'].get_value()

            try:
                st = self.rowContainer.children[key].children['rightContainer'].children['startStepBox']. \
                    children['start'].get_value()
                st = UTCDateTime(st)
                step = int(self.rowContainer.children[key].children['rightContainer'].children['startStepBox']. \
                           children['step'].get_value())
                et = st + step

                # enable ppsd button if step > 3600
                self.rowContainer.children[key].children['rightContainer'].children['startStepBox']. \
                    children['ppsd'].set_enabled(step > 3600)
                self.rowContainer.children[key].children['rightContainer'].children['startStepBox']. \
                    children['ppsdLabel'].set_enabled(step > 3600)

                if (step <= 3600):
                    self.rowContainer.children[key].children['rightContainer'].children['startStepBox']. \
                        children['ppsd'].set_value(False)
                # end if

                # update plot
                self.rowContainer.children[key].children['rightContainer'].children['plot'] = gui.Label('Loading..')
                t = threading.Thread(target=setTraceImage,
                                     args=(nc, sc, lc, cc, st, et))
                t.start()
            except Exception as e:
                printException(e)
            # end try
        # end func

        def removeWidget(emitter):
            if(self.rowWidgetCount > 1):
                try:
                    self.rowContainer.remove_child(self.rowContainer.children[key])
                    self.rowWidgetCount -= 1
                    self.rowContainer.set_size(ROW_WIDGET_WIDTH * PADDING_FACTOR,
                                                  self.rowWidgetCount * ROW_WIDGET_HEIGHT * PADDING_FACTOR)
                    self.wrapperContainer.set_size(ROW_WIDGET_WIDTH * PADDING_FACTOR,
                                                   self.rowWidgetCount * ROW_WIDGET_HEIGHT * PADDING_FACTOR +
                                                   MAP_WIDGET_PADDING)
                except Exception as e:
                    printException(e)
                # end try
            # end if
        # end func

        key = str(uuid.uuid4())
        container = gui.HBox(width=ROW_WIDGET_WIDTH, height=ROW_WIDGET_HEIGHT, style={'margin': '0px auto'})
        leftContainer = gui.VBox(width=ROW_WIDGET_WIDTH * 0.25, height=ROW_WIDGET_HEIGHT * 0.3,
                                 style={'border': '1px solid blue', 'margin': '5px'})
        rightContainer = gui.VBox(width=ROW_WIDGET_WIDTH * 0.75, height=ROW_WIDGET_HEIGHT,
                                  style={'border': '1px solid blue', 'margin': '0px'})

        #############################################################
        # populate leftContainer
        #############################################################
        nslcBox = gui.HBox(width=ROW_WIDGET_WIDTH * 0.25, height=ROW_WIDGET_HEIGHT * 0.8, style={'margin': '0px auto'})
        nslcLabelBox = gui.HBox(width=ROW_WIDGET_WIDTH * 0.25, height=ROW_WIDGET_HEIGHT * 0.2,
                                style={'margin': '0px auto'})

        net = gui.DropDown.new_from_list(self.getNetworks())
        sta = gui.DropDown.new_from_list(self.getStations(net.get_value()))
        loc = gui.DropDown.new_from_list(self.getLocations(net.get_value(), sta.get_value()))
        cha = gui.DropDown.new_from_list(self.getChannels(net.get_value(), sta.get_value(), loc.get_value()))

        net.onchange.do(netChanged)
        sta.onchange.do(staChanged)
        loc.onchange.do(locChanged)
        cha.onchange.do(chaChanged)

        nslcBox.append({'net':net, 'sta':sta, 'loc':loc, 'cha':cha})
        nslcLabelBox.append([gui.Label('Network'),
                             gui.Label('Station'),
                             gui.Label('Location'),
                             gui.Label('Channel')])

        #############################################################
        # Add remove button
        #############################################################
        rmButton = gui.Button('Remove Row', height=100)

        rmButton.onclick.do(removeWidget)
        
        locStr, availStr = getMeta(net.get_value(), sta.get_value(), loc.get_value(), cha.get_value())
        locLabel = gui.Label(locStr)
        availLabel = gui.Label(availStr)
        leftContainer.append({'nslcLabelBox':nslcLabelBox, 'nslcBox':nslcBox,
                              'locLabel':locLabel, 'availLabel':availLabel, 'rmButton':rmButton})

        #############################################################
        # populate rightContainer
        #############################################################
        startStepBox = gui.HBox(width=ROW_WIDGET_WIDTH*0.40, height=ROW_WIDGET_HEIGHT*0.1,
                                style={'margin': '0px'})
        startStepLabelBox = gui.HBox(width=ROW_WIDGET_WIDTH*0.30, height=ROW_WIDGET_HEIGHT*0.1,
                                style={'margin': '0px'})
        startStepLabelBox.append([gui.Label('Start'), 
                                  gui.Label('Length (s)')])

        inp = gui.Input()
        step = gui.DropDown.new_from_list([DEFAULT_TRC_LENGTH, '3600', '14400', '86400'])
        step.select_by_value(DEFAULT_TRC_LENGTH)
        ppsd = gui.CheckBox()
        ppsd.set_enabled(False)
        ppsdLabel = gui.Label('PPSD')
        ppsdLabel.set_enabled(False)
        startStepBox.append({'start':inp, 'step':step, 'ppsdLabel':ppsdLabel, 'ppsd':ppsd})
        inp.onchange.do(startStepChanged)
        step.onchange.do(startStepChanged)
        ppsd.onchange.do(startStepChanged)

        rightContainer.append({'startStepLabelBox': startStepLabelBox, 'startStepBox': startStepBox,
                               'plot': gui.Label('Loading..')})

        t = threading.Thread(target=setTraceImage,
                             args=(net.get_value(), sta.get_value(), loc.get_value(), cha.get_value()))
        t.start()

        container.append({'leftContainer': leftContainer, 'rightContainer': rightContainer})

        self.rowWidgetCount += 1
        return container, key
    # end func

    def main(self, fds:FederatedASDFDataSet):
        def addWidget(emitter):
            row, key = self.rowWidget()
            self.rowContainer.append(row, key)
            self.rowContainer.set_size(ROW_WIDGET_WIDTH* PADDING_FACTOR,
                                          self.rowWidgetCount * ROW_WIDGET_HEIGHT * PADDING_FACTOR)
            self.wrapperContainer.set_size(ROW_WIDGET_WIDTH * PADDING_FACTOR,
                                           self.rowWidgetCount * ROW_WIDGET_HEIGHT * PADDING_FACTOR +
                                           MAP_WIDGET_PADDING)
        # end func

        self.fds = fds
        # create master container
        self.rowContainer = gui.VBox(width=ROW_WIDGET_WIDTH * PADDING_FACTOR,
                                        height=ROW_WIDGET_HEIGHT * PADDING_FACTOR,
                                        style={'margin': '0px auto', 'overflow': 'scroll'})
        self.wrapperContainer = None
        self.rowWidgetCount = 0

        # Create first row widget
        row, key = self.rowWidget()
        self.rowContainer.append(row, key)

        # Button for adding more rows
        addButtonBox = gui.HBox(width=ROW_WIDGET_WIDTH*0.6, height=ROW_WIDGET_HEIGHT*0.1, style={'margin': '0px'})
        addButton = gui.Button('Add Row')
        addButton.onclick.do(addWidget)
        addButtonBox.append({'dummyLabel': gui.Label(''), 'addButton': addButton})

        # Create map widget
        map, mapKey = self.mapWidget()

        # container returned contains the master-container
        container = gui.VBox(width=ROW_WIDGET_WIDTH * PADDING_FACTOR,
                             height=ROW_WIDGET_HEIGHT * PADDING_FACTOR + MAP_WIDGET_PADDING,
                             style={'margin': '5px auto', 'overflow': 'scroll'})
        # title
        title = gui.Label('HiPerSeis FederatedASDF Viewer', style={'font-size': '30px'})

        # notes
        notesBox = gui.HBox(width=ROW_WIDGET_WIDTH*0.4, height=ROW_WIDGET_HEIGHT*0.1, style={'margin': '0px'})
        notes = gui.Label('Notes: PPSD plots, generated based on a flat instrument response '
                          'and normalized trace amplitudes, require traces longer than 1 hr.',
                          style={'font-size': '10px'})
        notesBox.append([gui.Label(''), notes])

        container.append({'title': title,
                          mapKey: map,
                          'notesBox': notesBox,
                          'addButtonBox': addButtonBox,
                          'rowContainer': self.rowContainer})

        # returning the root widget
        self.wrapperContainer = container

        return container
    # end func
# end class


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('asdf-source', required=True,
                type=click.Path(exists=True))
def process(asdf_source):
    """
    ASDF_SOURCE: Path to text file containing paths to ASDF files\n

    Example: python FederatedASDFViewer.py /path/to/asdf_files.txt

    An http server will then be started with the address shown in the
    terminal prompt, e.g:

    INFO:remi.server:Started httpserver http://0.0.0.0:1122/

    Copy the http address and paste it in a web browser window to start
    the viewer.
    """

    fds = FederatedASDFDataSet(asdf_source, single_threaded_access=False)

    # starts the webserver
    start(DataViewer, address='0.0.0.0', port=1122, start_browser=False,
          update_interval=0.1, userdata=(fds,))
# end func

if (__name__ == '__main__'):
    process()
# end if
