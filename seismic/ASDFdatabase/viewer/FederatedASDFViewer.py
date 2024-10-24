"""
Description:
    Implements a web browser-based, nominally interactive GUI for viewing data holdings
    in a collection of ASDF files.

References:

CreationDate:   07/02/23
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     08/09/23   RH
"""

import remi
import remi.gui as gui
from remi import start, App
import os, sys
import numpy as np
from obspy import UTCDateTime
import click
import uuid
import threading
from matplotlib.figure import Figure
import cartopy.crs as ccrs
from scipy.stats import circmean as cmean
from collections import defaultdict
from shapely import geometry
from seismic.misc import print_exception
import plotly.express as px
import plotly.offline as plyo
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
from seismic.ASDFdatabase.viewer.custom_ppsd import CustomPPSD
from seismic.ASDFdatabase.viewer.gui_components import FigureImage, PlotlyGraph

NULL_CHANNEL_CODE = 'XXX'
NULL_AVAILABILITY = 'Availability: 2100-01-01 - 1900-01-01'

ROW_WIDGET_WIDTH = 1000
ROW_WIDGET_HEIGHT = 500
MAP_WIDGET_WIDTH = 1000
MAP_WIDGET_HEIGHT = 500
MAP_WIDGET_PADDING = 700
PADDING_FACTOR = 1.1
TRACE_FIG_WIDTH = 7
TRACE_FIG_HEIGHT = 3.5
MAP_FIG_WIDTH = 750
MAP_FIG_HEIGHT = 500
DEFAULT_TRC_LENGTH = '600'

class DataViewer(App):
    def __init__(self, *args):
        self.widgets_that_requires_javascript_after_update=[]
        super(DataViewer, self).__init__(*args)
    # end func

    def getNetworks(self):
        result = self.nslc_dict.keys()
        return sorted(list(result))
    # end func

    def getStations(self, net):
        result = self.nslc_dict[net].keys()
        return sorted(list(result))
    # end func

    def getLocations(self, net, sta):
        result = self.nslc_dict[net][sta].keys()
        return sorted(list(result))
    # end func

    def getChannels(self, net, sta, loc):
        result = list(self.nslc_dict[net][sta][loc])
        # add a null entry to users get a chance to select the channel
        # which channel they want to see data for
        result.append(NULL_CHANNEL_CODE)
        return sorted(result)
    # end func

    def mapWidget(self):
        def setMapImage(nc):
            def get_plotly_html():
                stations = np.array(self.getStations(nc))

                lons = []
                lats = []
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
                # end if
                lons = np.array(lons)
                lats = np.array(lats)

                data = {'Lon': lons, 'Lat': lats, 'Station': stations}
                fig = px.scatter_geo(data, lat='Lat', lon='Lon',
                                     text='Station',
                                     projection='equirectangular',
                                     color_discrete_sequence=['red'], width=MAP_FIG_WIDTH,
                                     height=MAP_FIG_HEIGHT)
                fig.update_traces(textposition="top center",
                                  mode='markers+text',
                                  marker={'size': 5, 'symbol': 'triangle-down'})
                fig.update_layout(font=dict(family="Tahoma",
                                            size=10,  # Set the font size here
                                            color="Black"))

                output = plyo.plot(fig, output_type="div", include_plotlyjs=False)
                return output
            # end func

            pg = PlotlyGraph(self, get_plotly_html())
            self.wrapperContainer.children[key].children['rightContainer'].children['plot'] = pg
        # end func

        def mapNetChanged(emitter, value=None):
            nc = self.wrapperContainer.children[key].children['leftContainer'].children['nBox'].children['net'].get_value()
            
            # update plot
            self.wrapperContainer.children[key].children['rightContainer'].children['plot'] = gui.Label('Loading..')

            t = threading.Thread(target=setMapImage,
                                 args=(nc,))
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

                mapNetChanged(emitter, value)
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

        # network selection
        net = gui.DropDown.new_from_list(self.getNetworks())
        net.onchange.do(mapNetChanged)
        nBox.append({'net':net})
        nLabelBox.append(gui.Label('Network'))

        exportCoordinates = gui.Button('Export Coordinates', height=30, margin='1px auto')
        exportCoordinates.onclick.do(writeCoordinates)

        leftContainer.append({'nLabelBox': nLabelBox, 'nBox':nBox,
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
            #print(nc, sc, lc, cc, st, et)
            locStr = "Lon: {:.2f}, Lat: {:.2f}".format(lon, lat)

            availStr = None
            if(cc == NULL_CHANNEL_CODE):
                availStr = NULL_AVAILABILITY
            else:
                st, et = self.fds.get_global_time_range(nc, sc, lc, cc)
                availStr = "Availability: {} - {}".format(st.strftime('%Y-%m-%d'), et.strftime('%Y-%m-%d'))
            # end if

            return locStr, availStr
        # end func

        def setTraceImage(nc, sc, lc, cc, st=None, et=None):
            if(cc == NULL_CHANNEL_CODE): return # nothing to do for null channel-code

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
                print_exception(e)
            # end try
        # end func

        def chaChanged(emitter, value):
            nc = self.rowContainer.children[key].children['leftContainer'].children['nslcBox'].children['net'].get_value()
            sc = self.rowContainer.children[key].children['leftContainer'].children['nslcBox'].children['sta'].get_value()
            lc = self.rowContainer.children[key].children['leftContainer'].children['nslcBox'].children['loc'].get_value()
            cc = self.rowContainer.children[key].children['leftContainer'].children['nslcBox'].children['cha'].get_value()

            #print('chaChanged')

            # update metadata
            locStr, availStr = getMeta(nc, sc, lc, cc)
            self.rowContainer.children[key].children['leftContainer'].\
                children['locLabel'].set_text(locStr)
            self.rowContainer.children[key].children['leftContainer'].\
                children['availLabel'].set_text(availStr)

            if(cc == NULL_CHANNEL_CODE): return

            # update plot
            self.rowContainer.children[key].children['rightContainer'].set_enabled(True)
            self.rowContainer.children[key].children['rightContainer'].children['plot'] = gui.Label('Loading..')
            t = threading.Thread(target=setTraceImage,
                                 args=(nc, sc, lc, cc))
            t.start()
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
                self.rowContainer.children[key].children['rightContainer'].set_enabled(True)
                self.rowContainer.children[key].children['rightContainer'].children['plot'] = gui.Label('Loading..')
                t = threading.Thread(target=setTraceImage,
                                     args=(nc, sc, lc, cc, st, et))
                t.start()
            except Exception as e:
                print_exception(e)
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
                    print_exception(e)
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
                               'plot': gui.Label('')})
        rightContainer.set_enabled(False)

        t = threading.Thread(target=setTraceImage,
                             args=(net.get_value(), sta.get_value(), loc.get_value(), cha.get_value()))
        t.start()

        container.append({'leftContainer': leftContainer, 'rightContainer': rightContainer})

        self.rowWidgetCount += 1
        return container, key
    # end func

    def main(self, fds:FederatedASDFDataSet):
        # import plotly library:
        self.page.children['head'].add_child("plotly_import",
                                             '<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>\n')
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
        # populate net, sta, loc, cha dict
        self.nslc_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

        nslc_list = self.fds.get_nslc_list()
        for row in nslc_list:
            net, sta, loc, cha = row
            self.nslc_dict[net][sta][loc].append(cha)
        # end for

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
    def do_gui_update(self):
        """ This method gets called also by Timer, a new thread, and so needs to lock the update
        """
        with self.update_lock:
            changed_widget_dict = {}
            self.root.repr(changed_widget_dict)
            for widget in changed_widget_dict.keys():
                html = changed_widget_dict[widget]
                __id = str(widget.identifier)
                self._send_spontaneous_websocket_message(
                    remi.server._MSG_UPDATE + __id + ',' + remi.server.to_websocket(html))
            # end for
        # end with

        self._need_update_flag = False

        for wu in self.widgets_that_requires_javascript_after_update:
            for cw_html in changed_widget_dict.values():
                if wu.attributes['id'] in cw_html:
                    self.execute_javascript(wu.javascript)
                # end if
            # end for
        # end for

    # end func

    def onpageshow(self, *args):
        """ WebPage Event that occurs on webpage gets shown """
        super(DataViewer, self).onpageshow(*args)
        for wu in self.widgets_that_requires_javascript_after_update:
            wu.refresh()
        # end for
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
          update_interval=0, userdata=(fds,))
# end func

if (__name__ == '__main__'):
    process()
# end if
