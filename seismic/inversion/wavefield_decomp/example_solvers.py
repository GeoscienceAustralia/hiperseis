#!/usr/bin/env python
# coding: utf-8
"""
Examples of using the optimization (minimization) solver functions.
"""

import logging

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                    level=logging.INFO,
                    datefmt='%Y-%m-%d %H:%M:%S')

import numpy as np
import pandas as pd
from scipy.optimize import Bounds
import matplotlib.pyplot as plt
from landscapes.single_objective import sphere, himmelblau, easom, rosenbrock, rastrigin
import seaborn as sb

from seismic.inversion.wavefield_decomp.solvers import optimize_minimize_mhmcmc_cluster


# pylint: disable=invalid-name, missing-function-docstring, logging-format-interpolation, too-many-statements

def plot_2d(solution, title=''):
    # Visualize results.
    plt.figure(figsize=(16, 8))
    ndims = solution.x.shape[-1]
    for i in range(ndims):
        plt.subplot(1, ndims, i + 1)
        plt.bar(solution.bins[i, :-1] + 0.5*np.diff(solution.bins[i, :]), solution.distribution[i, :])
        # plt.xticks(hist.bins[i, :])
        plt.xlabel('x{}'.format(i))
        plt.ylabel('Counts')
    # end for
    plt.show()

    plt.figure(figsize=(8, 8))
    if solution.samples is not None:
        plt.scatter(solution.samples[:, 0], solution.samples[:, 1], c='#202020', alpha=0.05, s=5)
    # end if
    x = solution.x.copy()
    if len(x) == 1:
        np.expand_dims(x, axis=0)
    for i, _x in enumerate(x):
        color = 'C' + str(i)
        cluster = solution.clusters[i]
        plt.scatter(cluster[:, 0], cluster[:, 1], c=color, s=5, alpha=0.3)
        plt.scatter(_x[0], _x[1], marker='x', s=200, c=color, alpha=0.9)
        plt.scatter(_x[0], _x[1], marker='o', s=320, facecolors='none', edgecolors=color, alpha=0.9, linewidth=2)
    # end for
    plt.xlim(solution.bounds.lb[0], solution.bounds.ub[0])
    plt.ylim(solution.bounds.lb[1], solution.bounds.ub[1])
    plt.axis('equal')
    plt.xlabel('x0')
    plt.ylabel('x1')
    if title:
        plt.title(title)
    # end if
    plt.grid(color='#80808080', linestyle=':')
    # plt.savefig(objective.__name__.replace('<', '').replace('>', '') + '.png', dpi=300)
    plt.show()
# end func


def plot_Nd(soln, title='', scale=1.0):
    """Plotting routine for N-dimensional solution in grid format.

    :param soln: [description]
    :type soln: [type]
    :param title: [description], defaults to ''
    :type title: str, optional
    :param scale: [description], defaults to 1.0
    :type scale: float, optional
    :return: [description]
    :rtype: [type]
    """

    soln_alpha = 0.3
    samples_alpha = 0.05
    hist_alpha = 0.5
    axis_font_size = 12
    text_font_size = 10
    ndims = len(soln.bounds.lb)

    # Use PairGrid to set up grid and useful attributes of plot.
    df = pd.DataFrame(soln.samples, columns=['x' + str(i) for i in range(ndims)])
    p = sb.PairGrid(df, height=3.2*scale)
    # Plot samples (not actual solution, just samples of MCMC process) as grey background on off-diagonals.
    p = p.map_offdiag(plt.scatter, color='#808080', alpha=samples_alpha, s=2*scale**2)

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
            ax.set_title('x{} sample distribution'.format(row), y=0.9, color='#404040', fontsize=11*scale)
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
                    hjust = 'left'
                    hoffset = 0.02*x_range
                else:
                    hjust = 'right'
                    hoffset = -0.02*x_range
                # end if
                # Work out exact position on local y-axis, using full N-dimensional solution to minimize
                # overlap by project N-dimensional position onto the diagonal of the bounded space.
                bounds_diag = soln.bounds.ub - soln.bounds.lb
                denom = np.dot(bounds_diag, bounds_diag)
                y_pos_norm = np.dot(_x - soln.bounds.lb, bounds_diag)/denom
                assert 0.0 <= y_pos_norm <= 1.0
                y_pos = x_lim[0] + y_pos_norm*x_range
                if y_pos_norm >= 0.5:
                    vjust = 'bottom'
                else:
                    vjust = 'top'
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
                ax.scatter(cluster[:, col], cluster[:, row], c=color, s=2*scale**2, alpha=soln_alpha)
            # end for
            # Add dotted grid
            p.axes[row, col].grid(color='#80808080', linestyle=':')
        # end if
    # end for
    # Overall plot title
    if title:
        plt.suptitle(title, y=1.05, fontsize=16*scale)
    # end if
    # TODO: figure out how to adjust line labels in the y-direction so not overlapping,
    #   see https://support.sisense.com/hc/en-us/community/posts/360037908374-Getting-Around-Overlapping-Data-Labels-With-Python

    return p, diag_hist_ax, adjustable_text
# end func


def plot_Nd_joint(soln, title='', scale=1.0):
    # TBD: JointPlot format for N-dimensional solution.
    pass
# end func


def example_2d():
    # 2D test functions as per https://en.wikipedia.org/wiki/Test_functions_for_optimization

    def bi_quadratic(x, mu, cov):
        # The number returned from the function must be non-negative.
        # The exponential of the negative of this value is the probablity.
        sqrt_arg = np.matmul(np.matmul((x - mu).T, cov), x - mu)
        assert sqrt_arg >= 0
        x2fac = np.sqrt(sqrt_arg)
        return x2fac
    # end func

    logging.info("Solving sphere function")
    bounds = Bounds(np.array([-3, -3]), np.array([3, 3]))
    soln = optimize_minimize_mhmcmc_cluster(
        sphere, bounds, burnin=10000, maxiter=50000, collect_samples=10000, logger=logger)
    if soln.success:
        logging.info("Solution:\n{}".format(soln.x))
        plot_2d(soln, title='Sphere function minima')
    # end if

    logging.info("Solving himmelblau function")
    bounds = Bounds(np.array([-5, -5]), np.array([5, 5]))
    soln = optimize_minimize_mhmcmc_cluster(
        himmelblau, bounds, burnin=10000, maxiter=50000, collect_samples=10000, target_ar=0.3, T=10, logger=logger)
    if soln.success:
        logging.info("Solution:\n{}".format(soln.x))
        plot_2d(soln, title='Himmelblau function minima')
    # end if

    logging.info("Solving easom function")
    bounds = Bounds(np.pi + np.array([-3, -2]), np.pi + np.array([3, 4]))
    soln = optimize_minimize_mhmcmc_cluster(easom, bounds, burnin=10000, maxiter=50000, collect_samples=10000, T=0.2)
    if soln.success:
        logging.info("Solution:\n{}".format(soln.x))
        plot_2d(soln, title='Easom function minima')
    # end if

    logging.info("Solving rosenbrock function")
    bounds = Bounds(np.array([-4, -4]), np.array([4, 4]))
    soln = optimize_minimize_mhmcmc_cluster(rosenbrock, bounds, burnin=10000, maxiter=50000, collect_samples=10000,
                                            T=1000)
    if soln.success:
        logging.info("Solution:\n{}".format(soln.x))
        plot_2d(soln, title='Rosenbrock function minima')
    # end if

    logging.info("Solving rastrigin function")
    bounds = Bounds(np.array([-4, -4]), np.array([4, 4]))
    soln = optimize_minimize_mhmcmc_cluster(lambda xy: rastrigin(xy) - 2, bounds, burnin=10000, maxiter=50000,
                                            collect_samples=10000, T=5, N=5, logger=logger)
    if soln.success:
        logging.info("Solution:\n{}".format(soln.x))
        plot_2d(soln, title='Rastrigin function minima')
    # end if

    # Custom test function
    logging.info("Solving bi-variate quadratic function")
    mu = np.array([0, 1])
    cov = np.array([[5, -6.0], [-6.0, 20.0]])
    fixed_args = (mu, cov)
    bounds = Bounds(np.array([-3, -2]), np.array([3, 4]))
    soln = optimize_minimize_mhmcmc_cluster(bi_quadratic, bounds, fixed_args, burnin=10000, maxiter=50000,
                                            rnd_seed=20200220, logger=logger)
    if soln.success:
        logging.info("Solution:\n{}".format(soln.x))
        plot_2d(soln, title='Custom bi-variate quadratic function minima')
    # end if
# end func


def example_3d():
    logging.info("Solving 3D sphere function")
    translate = np.array([-1, 1.5, 2.5])
    bounds = Bounds(np.array([-3, -3, -3]) + translate, np.array([3, 3, 3]) + translate)
    soln = optimize_minimize_mhmcmc_cluster(lambda xyz: sphere(xyz - translate),
                                            bounds, burnin=10000, maxiter=50000, collect_samples=10000, logger=logger)
    if soln.success:
        logging.info("Solution:\n{}".format(soln.x))
        p, _, _ = plot_Nd(soln, title='Sphere 3D function minima', scale=1.2)
        p.savefig('sphere_3d_viz_example.png', dpi=300)
        plt.show()
    # end if

    logging.info("Solving 3D Rastrigin function")
    bounds = Bounds(np.array([-5, -5, -5]), np.array([5, 5, 5]))
    soln = optimize_minimize_mhmcmc_cluster(
        rastrigin, bounds, burnin=10000, maxiter=50000, collect_samples=10000, T=5, N=5, logger=logger)
    if soln.success:
        logging.info("Solution:\n{}".format(soln.x))
        p, _, _ = plot_Nd(soln, title='Rastrigin function minima')
        p.savefig('rastrigin_3d_viz_example.png', dpi=300)
        plt.show()
    # end if

# end func


def example_4d():
    pass
# end func


def main():

    # example_2d()
    example_3d()
    # example_4d()

# end func


if __name__ == "__main__":
    logger = logging.getLogger(__name__)
    main()
# end if
