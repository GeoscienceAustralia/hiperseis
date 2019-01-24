from rf import RFStream
from pstplot import pst_plot_rf

class PST_RFStream(RFStream):

    def __init__(self, rfstream):
        self.rfstr = rfstream

    def plot_rf(self, *args, **kwargs):
        return pst_plot_rf(self.rfstr, *args, **kwargs)
