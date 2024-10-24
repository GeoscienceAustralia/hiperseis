"""
Description:
    Implements custom widgets based on standard REMI widgets

References:

CreationDate:   24/10/24
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     24/10/24   RH
"""

import io
import time
import matplotlib
from matplotlib.backends.backend_agg import FigureCanvasAgg
import re
import remi.gui as gui

font = {'family' : 'normal',
        'size'   : 8}

matplotlib.rc('font', **font)
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)

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

class PlotlyGraph(gui.Widget):
    def __init__(self, appInstance, plotlyHTML=None, *args, **kwargs):
        super(PlotlyGraph, self).__init__(*args, **kwargs)
        self.appInstance = appInstance
        self.javascript = ""
        self._classes = set()

        if (plotlyHTML):
            self.setPlotlyHTML(plotlyHTML)
        # end if
    # end func

    def setPlotlyHTML(self, plotlyHTML):
        self.HTML = plotlyHTML

        # Replace the div id randomly choosen by plotly by the one given in argument (=the one used by remi):
        div_id = str(id(self))
        plotly_original_div_id = re.search(r"<div id=\"(.*)\" class", self.HTML).groups()[0]
        self.HTML = self.HTML.replace(plotly_original_div_id, div_id)

        # Get the plotly div class and add it to the remi widget:
        # div=re.search(r"(<div id=.*></div>\n)",self.HTML).groups()[0]
        div = re.search(r"(<div id=.*></div>)", self.HTML).groups()[0]
        div_class = re.search(r'class="(.*)" style', div).groups()[0]

        if (div_class not in self.attributes['class'].split(' ')):
            self.add_class(div_class)

        # Get the div style from plotly and the rules to the remi widget:
        div_style = re.search(r'style="(.*)">', div).groups()[0]
        # (Remi requires a dict for styling, create it here:)
        div_style = div_style.replace(" ", "").replace(";", ":").split(":")  # a list of [name1,value1,name2,value2...]
        div_style_dict = {}
        for i in range(0, len(div_style) - 1, 2):  # fill the dict:
            div_style_dict[div_style[i]] = div_style[i + 1]
        self.set_style(div_style_dict)

        # Get javascript code actually creating the graph:
        self.javascript = self.HTML.split('<script type="text/javascript">')[1].split('</script>')[0]

        # Register this Widget to call its javascript after App.do_gui_update() if updated
        self.appInstance.widgets_that_requires_javascript_after_update.append(self)
    # end func

    def refresh(self):
        self.appInstance.execute_javascript(self.javascript)
    # end func
# end class
