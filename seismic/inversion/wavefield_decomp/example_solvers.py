#!/usr/bin/env python
# coding: utf-8
"""
Examples of using the optimization (minimization) solver functions.
"""

import logging

#pylint: disable=wrong-import-position

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                    level=logging.INFO,
                    datefmt='%Y-%m-%d %H:%M:%S')

import numpy as np
from scipy.optimize import Bounds
import matplotlib.pyplot as plt
from landscapes.single_objective import sphere, himmelblau, easom, rosenbrock, rastrigin, styblinski_tang

from seismic.inversion.wavefield_decomp.solvers import optimize_minimize_mhmcmc_cluster
from seismic.inversion.wavefield_decomp.wfd_plot import plot_Nd

# pylint: disable=invalid-name, missing-function-docstring, logging-format-interpolation, too-many-statements


def plot_Nd_joint(soln, title='', scale=1.0):
    # TBD: Develop JointPlot format for N-dimensional solution.
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
        plot_Nd(soln, title='Sphere function minima', scale=1.3)
    # end if

    logging.info("Solving himmelblau function")
    bounds = Bounds(np.array([-5, -5]), np.array([5, 5]))
    soln = optimize_minimize_mhmcmc_cluster(
        himmelblau, bounds, burnin=10000, maxiter=50000, collect_samples=10000, target_ar=0.3, T=10, logger=logger)
    if soln.success:
        logging.info("Solution:\n{}".format(soln.x))
        plot_Nd(soln, title='Himmelblau function minima', scale=1.3)
    # end if

    logging.info("Solving easom function")
    bounds = Bounds(np.pi + np.array([-3, -2]), np.pi + np.array([3, 4]))
    soln = optimize_minimize_mhmcmc_cluster(easom, bounds, burnin=10000, maxiter=50000, collect_samples=10000, T=0.2)
    if soln.success:
        logging.info("Solution:\n{}".format(soln.x))
        plot_Nd(soln, title='Easom function minima', scale=1.3)
    # end if

    logging.info("Solving rosenbrock function")
    bounds = Bounds(np.array([-4, -4]), np.array([4, 4]))
    soln = optimize_minimize_mhmcmc_cluster(rosenbrock, bounds, burnin=10000, maxiter=50000, collect_samples=10000,
                                            T=1000)
    if soln.success:
        logging.info("Solution:\n{}".format(soln.x))
        plot_Nd(soln, title='Rosenbrock function minima', scale=1.3)
    # end if

    logging.info("Solving rastrigin function")
    bounds = Bounds(np.array([-4, -4]), np.array([4, 4]))
    soln = optimize_minimize_mhmcmc_cluster(lambda xy: rastrigin(xy) - 2, bounds, burnin=10000, maxiter=50000,
                                            collect_samples=10000, T=5, N=5, logger=logger)
    if soln.success:
        logging.info("Solution:\n{}".format(soln.x))
        plot_Nd(soln, title='Rastrigin function minima', scale=1.3)
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
        plot_Nd(soln, title='Custom bi-variate quadratic function minima', scale=1.3)
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
        p, _, _ = plot_Nd(soln, title='Sphere 3D function minima', scale=1.0)
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
    logging.info("Solving 4D Styblinski-Tang function")
    bounds = Bounds(np.array([-5, -5, -5, -5]), np.array([5, 5, 5, 5]))
    soln = optimize_minimize_mhmcmc_cluster(
        styblinski_tang, bounds, burnin=100000, maxiter=500000, collect_samples=20000, T=20, N=5, logger=logger,
        target_ar=0.5, ar_tolerance=0.04)
    if soln.success:
        logging.info("Solution:\n{}".format(soln.x))
        p, _, _ = plot_Nd(soln, title='Styblinski-Tang function minima', scale=0.7)
        p.savefig('styblinski_tang_4d_viz_example.png', dpi=300)
        plt.show()
    # end if

# end func


def main():

    example_2d()
    example_3d()
    example_4d()

# end func


if __name__ == "__main__":
    logger = logging.getLogger(__name__)
    main()
# end if
