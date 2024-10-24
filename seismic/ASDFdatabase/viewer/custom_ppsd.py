"""
Description:
    Implements PPSD based on a flat response

References:

CreationDate:   08/09/23
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     24/10/24   RH
"""

import numpy as np
from obspy.core import Trace, Stream
from obspy.signal import PPSD

class CustomPPSD(PPSD):
    def __init__(self, stats, skip_on_gaps=False,
                 db_bins=(-150, 50, 1.), ppsd_length=3600.0, overlap=0.5,
                 special_handling=None, period_smoothing_width_octaves=1.0,
                 period_step_octaves=0.125, period_limits=None,
                 **kwargs):

        # flat response
        metadata = paz = {'sensitivity': 1.0,
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
            if (tr.stats.npts > 0):
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
