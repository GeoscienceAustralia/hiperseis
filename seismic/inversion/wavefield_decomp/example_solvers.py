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
from scipy.optimize import Bounds
import matplotlib.pyplot as plt
from landscapes.single_objective import sphere, himmelblau, easom, rosenbrock, rastrigin
import seaborn as sb

from seismic.inversion.wavefield_decomp.solvers import optimize_minimize_mhmcmc_cluster


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


def plot_3d(soln, title=''):
    ndims = soln.x.shape[-1]
    height = 5
    aspect_ratio = 1
    figsize = (ndims*height*aspect_ratio, ndims*height)
    fig, axes = plt.subplots(ndims, ndims, figsize=figsize, sharex='col', sharey='row', squeeze=False)
    # Plot samples if available in grey. Then overlay with scatter plot of clusters.
    # Plot full histogram on the diagonal for each variable.
    # Label scatter plots and histograms with actual solution values colour coded by cluster.
    # Set axes limits from the pre-defined bounds.

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
        plot_3d(soln, title='Sphere 3D function minima')
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
