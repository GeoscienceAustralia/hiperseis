#!/usr/bin/env python
# coding: utf-8
"""
Plotting helper functions for Wavefield Decomposition module.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb

# pylint: disable=invalid-name


def plot_Esu_space(H, k, Esu, title=None, savefile_name=None, show=True, c_range=(None, None), decorator=None):
    """
    Plot SU energy as function of H-k.

    :param H: Grid of H coordinates corresponding to Esu point values.
    :param k: Grid of k coordinates corresponding to Esu point values.
    :param Esu: Energy values on H-k grid.
    :param title: Plot title [OPTIONAL]
    :param savefile_name: Output file in which to save plot [OPTIONAL]
    :param show: If True, display the image and block until it is closed.
    :param c_range: Custom range of Esu to contour (min, max values)
    :param decorator: Callback function to customize plot.
    :return: None
    """
    colmap = 'plasma'
    plt.figure(figsize=(16, 12))
    min_E, max_E = (np.nanmin(Esu) if c_range[0] is None else c_range[0],
                    np.nanmax(Esu) if c_range[1] is None else c_range[1])
    plt.contourf(k, H, Esu, levels=np.linspace(min_E, max_E, 51), cmap=colmap)
    plt.colorbar()
    plt.contour(k, H, Esu, levels=np.linspace(min_E, max_E, 11), colors='k', linewidths=1, antialiased=True)
    plt.xlabel('Crustal $\\kappa$', fontsize=14)
    plt.ylabel('Crustal $H$ (km)', fontsize=14)
    plt.tick_params(right=True, labelright=True, which='both')
    plt.tick_params(top=True, labeltop=True, which='both')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.minorticks_on()
    plt.xlim(np.min(k), np.max(k))
    plt.ylim(np.min(H), np.max(H))
    plt.grid(linestyle=':', color="#80808080")
    if title is not None:
        plt.title(title, fontsize=20, y=1.05)
    # end if
    if decorator is not None:
        decorator(plt.gca())
    # end if
    if savefile_name is not None:
        plt.savefig(savefile_name, dpi=300)
    # end if
    if show:
        plt.show()
    # end if
    plt.close()
# end func


def plot_Nd(soln, title='', scale=1.0, vars=None):
    """Plotting routine for N-dimensional solution in grid format. Diagonal contains histograms of solution
    samples for each variable, and off-diagonal contains pair-wise scatter plots of sample points plus
    solution clusters. The distinct solutions coordinates are colour coded and labelled on the histograms.

    :param soln: Solution container returned from optimize_minimize_mhmcmc_cluster()
    :type soln: scipy.optimize.OptimizeResult with additional custom fields
    :param title: Overall plot title, defaults to '' (no title)
    :type title: str, optional
    :param scale: Scale size of overall plot, defaults to 1.0. Adjust this depending on dimensionality or desired size.
    :type scale: float, optional
    :param vars: List of names of variables, should be same length as dimensionality of solution
    :type vars: list(str)
    :return: Tuple containing the PairGrid, list of axes for the secondary (histogram) axes in the diagonal of
        the grid, and list of text items representing solution labels.
    :rtype: tuple(seaborn.PairGrid, list(matplotlib.axes.Axes), list(matplotlib.text.Text))
    """

    soln_alpha = 0.3
    samples_alpha = 0.05
    hist_alpha = 0.5
    axis_font_size = 12
    text_font_size = 10
    ndims = len(soln.bounds.lb)

    if vars is None:
        vars = ['x' + str(i) for i in range(ndims)]
    # end if

    # Use PairGrid to set up grid and useful attributes of plot.
    df = pd.DataFrame(soln.samples, columns=vars)
    p = sb.PairGrid(df, height=3.2*scale)
    # Plot samples (not actual solution, just samples of MCMC process) as grey background on off-diagonals.
    p = p.map_offdiag(plt.scatter, color='#808080', alpha=samples_alpha, s=2*scale**2, rasterized=True)

    diag_hist_ax = []
    row_idx, col_idx = np.indices((ndims, ndims))
    adjustable_text = []  # Collect line text labels
    for row, col in zip(row_idx.flat, col_idx.flat):
        if row == col:
            # Diagonal plots - use full sample histogram.
            # axd is the original diagonal axes created by PairGrid
            axd = p.axes[row, row]
            # Set label sizes
            axd.tick_params(labelsize=axis_font_size*scale)
            axd.xaxis.label.set_size(axis_font_size*scale)
            axd.yaxis.label.set_size(axis_font_size*scale)
            # Duplicate axes with separate, hidden vertical scale for the histogram.
            ax = axd.twinx()
            ax.set_axis_off()
            # Plot full samples histogram.
            deltas = np.diff(soln.bins[row])
            ax.bar(soln.bins[row, :-1] + 0.5*deltas, soln.distribution[row], color='#808080',
                   alpha=hist_alpha, width=np.min(deltas))
            ax.set_title('{} sample distribution'.format(vars[row]), y=0.9, color='#404040', fontsize=11*scale)
            # Lock axes ranges to parameter ranges
            ax.set_xlim(soln.bounds.lb[row], soln.bounds.ub[row])
            # Add vertical lines to histogram to indication solution locations and label value.
            for i, _x in enumerate(soln.x):
                color = 'C' + str(i)
                ax.axvline(_x[row], color=color, linestyle='--', linewidth=1.2*scale)
                # Sneakily use the axd axes for labelling, as it has same scale on x- and y- axes,
                # which we can use to make sure the labels for multiple solutions are at different heights.
                # Work out exact position on local x-axis.
                x_lim = ax.get_xlim()
                x_range = x_lim[1] - x_lim[0]
                if (_x[row] - x_lim[0])/x_range >= 0.5:
                    hjust = 'right'
                    hoffset = -0.02*x_range
                else:
                    hjust = 'left'
                    hoffset = 0.02*x_range
                # end if
                # Work out exact position on local y-axis, using full N-dimensional solution to minimize
                # overlap by project N-dimensional position onto the diagonal of the bounded space.
                bounds_diag = soln.bounds.ub - soln.bounds.lb
                denom = np.dot(bounds_diag, bounds_diag)
                y_pos_norm = np.dot(_x - soln.bounds.lb, bounds_diag)/denom
                assert 0.0 <= y_pos_norm <= 1.0
                y_pos = x_lim[0] + 0.9*y_pos_norm*x_range
                if y_pos_norm >= 0.5:
                    vjust = 'top'
                else:
                    vjust = 'bottom'
                # end if
                t = axd.text(_x[row] + hoffset, y_pos, '{:.3f}'.format(_x[row]), ha=hjust, va=vjust, color=color,
                             fontsize=text_font_size*scale, fontstyle='italic', fontweight='semibold', zorder=100+i)
                adjustable_text.append(t)
            # end for
            diag_hist_ax.append(ax)
        else:
            # Off-diagonal plots.
            ax = p.axes[row, col]
            # Set label sizes
            ax.tick_params(labelsize=axis_font_size*scale)
            ax.xaxis.label.set_size(axis_font_size*scale)
            ax.yaxis.label.set_size(axis_font_size*scale)
            # Plot distinct solution clusters
            for i, cluster in enumerate(soln.clusters):
                color = 'C' + str(i)
                ax.scatter(cluster[:, col], cluster[:, row], c=color, s=2*scale**2, alpha=soln_alpha, rasterized=True)
            # end for
            # Lock axes ranges to parameter ranges
            ax.set_xlim(soln.bounds.lb[col], soln.bounds.ub[col])
            ax.set_ylim(soln.bounds.lb[row], soln.bounds.ub[row])
            # Add dotted grid
            p.axes[row, col].grid(color='#80808080', linestyle=':')
        # end if
    # end for
    # Overall plot title
    if title:
        plt.suptitle(title, y=0.96, fontsize=16*scale)
    # end if

    plt.subplots_adjust(left=0.125, top=0.9, bottom=0.10, right=0.95)

    # TODO: figure out how to adjust line labels in the y-direction so not overlapping,
    #   see https://support.sisense.com/hc/en-us/community/posts/360037908374-Getting-Around-Overlapping-Data-Labels-With-Python

    return p, diag_hist_ax, adjustable_text
# end func
