# -*- coding: utf8 -*-
import matplotlib.pyplot as plt
from traceselector import TraceSelector
import numpy as np
from matplotlib.ticker import (AutoMinorLocator, MaxNLocator)

def pst_plot_rf(stream, fname=None, fig_width=7., trace_height=0.5,
            stack_height=0.5, scale=1, fillcolors=(None, None), trim=None,
            info=(('back_azimuth', u'baz (°)', 'b'),
                  ('distance', u'dist (°)', 'r'))):

    widget_arr = []
    if len(stream) == 0:
        return
    if trim:
        stream = stream.slice2(*trim, reftime='onset')
    N = len(stream)
    # calculate lag times
    stats = stream[0].stats
    times = stream[0].times() - (stats.onset - stats.starttime)
    # calculate axes and figure dimensions
    # big letters: inches, small letters: figure fraction
    H = trace_height
    HS = stack_height
    FB = 0.5
    FT = 0.2
    DW = 0.1
    FH = H * (N + 2) + HS + FB + FT + DW
    h = H / FH
    hs = HS / FH
    fb = FB / FH
    ft = FT / FH
    FL = 0.5
    FR = 0.2
    FW = fig_width
    FW3 = 0.8
    FW2 = FW - FL - FR - (DW + FW3) * bool(info)
    fl = FL / FW
    fr = FR / FW
    fw2 = FW2 / FW
    fw3 = FW3 / FW
    # init figure and axes
    fig = plt.figure(figsize=(FW, FH))
    ax1 = fig.add_axes([fl, fb, fw2, h * (N + 2)])
    if info:
        ax3 = fig.add_axes(
            [1 - fr - fw3, fb, fw3, h * (N + 2)], sharey=ax1)
        info = list(info)
        info[0] = [ax3] + list(info[0])
        if len(info) > 1:
            ax4 = ax3.twiny()
            info[1] = [ax4] + list(info[1])
    # plot individual receiver functions

    def _plot(ax, t, d, i):
        c1, c2 = fillcolors
        if c1:
            ax.fill_between(t, d + i, i, where=d >= 0, lw=0., facecolor=c1)
        if c2:
            ax.fill_between(t, d + i, i, where=d < 0, lw=0., facecolor=c2)
        if i > 0:
            if i%2 > 0:
                t1 = ax.text(-5.4, (i-.4), str(i).zfill(3), fontsize=5.5).set_bbox(dict(boxstyle='round', facecolor='green', alpha=0.7))
                trsel = TraceSelector(ax, t1)
                trsel.set_active(True)
                ax.plot(t, d + i, 'g')
            else:
                t1 = ax.text(20.1, (i-.4), str(i).zfill(3), fontsize=5.5).set_bbox(dict(boxstyle='round', facecolor='red', alpha=0.7))
                trsel = TraceSelector(ax, t1)
                trsel.set_active(True)
                ax.plot(t, d + i, 'r')
    max_ = max(np.max(np.abs(tr.data)) for tr in stream)
    for i, tr in enumerate(stream):
        _plot(ax1, times, tr.data / max_ * scale, i + 1)
    #clabels = [str(i) for i in range(len(stream))]
#    clabels = ('',)*len(stream)
#    cvals = (False,)*len(stream)
#    select_arr = CheckButtons(ax1, clabels, cvals)
    # plot right axes with header information
    for ax, header, label, color in info:
        data = [tr.stats[header] for tr in stream]
        ax.plot(data, 1 + np.arange(len(stream)), '.' + color, mec=color)
        ax.set_xlabel(label, color=color, size='small')
        if header == 'back_azimuth':
            ax.set_xticks(np.arange(5) * 90)
            ax.set_xticklabels(['0', '', '180', '', '360'], size='small')
        else:
            ax.xaxis.set_major_locator(MaxNLocator(4))
            for l in ax.get_xticklabels():
                l.set_fontsize('small')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
    # set x and y limits
    ax1.set_xlim(times[0], times[-1])
    ax1.set_ylim(-0.5, N + 1.5)
    ax1.set_yticklabels('')
    ax1.set_xlabel('time (s)')
    ax1.xaxis.set_minor_locator(AutoMinorLocator())

    # plot stack
    stack = stream.stack()
    if len(stack) > 1:
        warnings.warn('Different stations or channels in one RF plot.')
    elif len(stack) == 1:
        ax2 = fig.add_axes([fl, 1 - ft - hs, fw2, hs], sharex=ax1)
        _plot(ax2, times, stack[0].data, 0)
        for l in ax2.get_xticklabels():
            l.set_visible(False)
        ax2.yaxis.set_major_locator(MaxNLocator(4))
        for l in ax2.get_yticklabels():
            l.set_fontsize('small')
        # annotate plot with seed id
        bbox = dict(boxstyle='round', facecolor='white', alpha=0.8, lw=0)
        text = '%s traces  %s' % (len(stream), stack[0].id)
        ax2.annotate(text, (1 - 0.5 * fr, 1 - 0.5 * ft),
                     xycoords='figure fraction', va='top', ha='right',
                     bbox=bbox, clip_on=False)
    # save plot
    if fname:
        fig.savefig(fname)
        plt.close(fig)
    else:
        return fig

