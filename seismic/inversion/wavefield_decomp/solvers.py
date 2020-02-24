#!/usr/bin/env python
# coding: utf-8
"""
Objective function minimization solvers.
"""

import numpy as np
import copy
from collections import deque

from scipy.optimize import Bounds
# from scipy.optimize import OptimizeResult
from sortedcontainers import SortedList
from sklearn.cluster import dbscan

# DEBUG ONLY
import matplotlib.pyplot as plt


# pylint: disable=invalid-name

class SolverGlobalMhMcmc:
    """
    Drop-in custom solver for scipy.optimize.minimize, based on Metrolpolis-Hastings Monte Carlo
    Markov Chain random walk with burn-in and adaptive acceptance rate, followed by N-dimensional
    clustering.

    Rather than returning one global solution, the solver returns the N best ranked local solutions.
    It also returns a probability distribution of each unknown based on Monte Carlo statistics.

    """
    pass
# end class


class HistogramIncremental():
    """Class to incrementally accumulate N-dimensional histogram stats at runtime.
    """
    def __init__(self, bounds, nbins=20):
        self.ndims = len(bounds.lb)
        assert len(bounds.ub) == self.ndims
        self.bins = np.linspace(bounds.lb, bounds.ub, nbins + 1).T
        self.hist = np.zeros((self.ndims, nbins), dtype=int)
    # end func

    def __iadd__(self, x):
        assert len(x) == self.ndims
        for i, _x in enumerate(x):
            idx = np.digitize(_x, self.bins[i, :])
            self.hist[i, idx - 1] += 1
        # end for
        return self
    # end func

# end class


class BoundedRandNStepper():
    """Step one dimensional at a time using normal distribution of steps.
    """
    def __init__(self, bounds, initial_step=None):
        self.bounds = bounds
        self.range = bounds.ub - bounds.lb
        if initial_step is not None:
            self._stepsize = initial_step
        else:
            self._stepsize = 0.15*(bounds.ub - bounds.lb)
        # end if
    # end func

    def __call__(self, x):
        ndims = len(x)
        while True:
            dim = np.random.randint(0, ndims)
            x_new = x[dim] + self._stepsize[dim]*np.random.randn()
            if self.bounds.lb[dim] <= x_new <= self.bounds.ub[dim]:
                break
            # end if
        # end while
        result = x.copy()
        result[dim] = x_new
        return result
    # end func

    @property
    def stepsize(self):
        return self._stepsize.copy()
    # end func

    @stepsize.setter
    def stepsize(self, stepsize_new):
        # Limit step size adaption so that attempt to reach desired acceptance rate does not
        # diverge or vanish the stepping quanta.
        stepsize_new = np.maximum(np.minimum(stepsize_new, self.range), 1e-3*self.range)
        self._stepsize = stepsize_new
    # end func

# end class


class AdaptiveStepsize():
    """
    Class to implement adaptive stepsize. Pulled from scipy.optimize.basinhopping and adapted.
    """
    def __init__(self, stepper, accept_rate=0.5, interval=50, factor=0.9, ar_tolerance=0.05, verbose=True):
        self.takestep = stepper
        self.target_accept_rate = accept_rate
        self.interval = interval
        self.factor = factor
        self.tolerance = ar_tolerance
        self.verbose = verbose

        self.nstep = 0
        self.nstep_tot = 0
        self.naccept = 0
    # end func

    def __call__(self, x):
        return self.take_step(x)
    # end func

    def _adjust_step_size(self):
        old_stepsize = copy.copy(self.takestep.stepsize)
        accept_rate = float(self.naccept) / self.nstep
        scaled_step = False
        if accept_rate > self.target_accept_rate + self.tolerance:
            # We're accepting too many steps. Take bigger steps.
            self.takestep.stepsize /= self.factor
            scaled_step = True
        elif accept_rate < self.target_accept_rate - self.tolerance:
            # We're not accepting enough steps.  Take smaller steps.
            self.takestep.stepsize *= self.factor
            scaled_step = True
        # end if
        if scaled_step and self.verbose:
            # TODO: Replace with logger
            print("adaptive stepsize: acceptance rate {} target {} new stepsize {} old stepsize {}"
                  .format(accept_rate, self.target_accept_rate, self.takestep.stepsize, old_stepsize))
        # end if
    # end func

    def take_step(self, x):
        self.nstep += 1
        self.nstep_tot += 1
        if self.nstep % self.interval == 0:
            self._adjust_step_size()
        return self.takestep(x)
    # end func

    def notify_accept(self):
        self.naccept += 1
    # end func

# end class


# Compose as global function first, then abtract out into classes.
def optimize_minimize_mhmcmc_cluster(objective, bounds, args=(), x0=None, T=1, N=3, burnin=1000000, maxiter=10000000,
                                     target_ar=0.4):
    """
    Minimize objective function and return up to N solutions.

    :param objective: Objective function to minimize
    :param bounds: Bounds of the parameter space.
    :param args: Any additional fixed parameters needed to completely specify the objective function.
    :param x0: Initial guess. If None, will be selected randomly and uniformly within the parameter bounds.
    :param T: The "temperature" parameter for the accept or reject criterion.
    :param N: Maximum number of minima to return
    :param burnin: Number of random steps to discard before starting to accumulate statistics.
    :param maxiter: Maximum number of steps to take (including burnin).
    :param target_ar: Target acceptance rate of point samples generated by stepping.
    :return: OptimizeResult containing solution(s) and solver data.
    """
    assert maxiter > burnin, "maxiter {} not greater than burnin steps {}".format(maxiter, burnin)
    main_iter = maxiter - burnin

    beta = 1.0/T

    # DEBUG ONLY:
    np.random.seed(20200220)

    if x0 is None:
        x0 = np.random.uniform(bounds.lb, bounds.ub)
    # end if
    assert np.all((x0 >= bounds.lb) & (x0 <= bounds.ub))
    x = x0.copy()
    funval = objective(x, *args)

    # Set up stepper with adaptive acceptance rate
    stepper = BoundedRandNStepper(bounds)
    stepper = AdaptiveStepsize(stepper, accept_rate=target_ar, interval=50)

    x_queue = deque(maxlen=10000)
    x_queue.append(x)
    rejected_randomly = 0
    accepted_burnin = 0
    for _ in range(burnin):
        x_new = stepper(x)
        funval_new = objective(x_new, *args)
        log_alpha = -(funval_new - funval)*beta
        if log_alpha > 0 or np.log(np.random.rand()) <= log_alpha:
            x = x_new
            funval = funval_new
            x_queue.append(x)
            stepper.notify_accept()
            accepted_burnin += 1
        elif log_alpha <= 0:
            rejected_randomly += 1
        # end if
    # end for
    ar = float(accepted_burnin)/burnin
    print("Burnin acceptance rate: {}".format(ar))

    pts_burnin = np.array(x_queue)
    del x_queue

    accepted = 0
    rejected_randomly = 0

    minima = SortedList(key=lambda rec: rec[1])
    hist = HistogramIncremental(bounds, nbins=100)
    # Cached a lot of potential minimum values, as these need to be clustered before return N results
    N_cached = 100*N
    for _ in range(main_iter):
        x_new = stepper(x)
        funval_new = objective(x_new, *args)
        log_alpha = -(funval_new - funval)*beta
        if log_alpha > 0 or np.log(np.random.rand()) <= log_alpha:
            x = x_new
            funval = funval_new
            minima.add((x, funval))
            if len(minima) > N_cached:
                minima.pop()
            # end if
            stepper.notify_accept()
            hist += x
            accepted += 1
        elif log_alpha <= 0:
            rejected_randomly += 1
        # end if
    # end for
    ar = float(accepted)/main_iter
    print("Acceptance rate: {}".format(ar))
    print("Best minima: {}".format(np.array([_mx[0] for _mx in minima[0:10]])))

    # plt.figure(figsize=(12, 20))
    # for i in range(hist.ndims):
    #     plt.subplot(2, 1, i + 1)
    #     plt.bar(hist.bins[i, :-1] + 0.5*np.diff(hist.bins[i, :]), hist.hist[i, :])
    #     # plt.xticks(hist.bins[i, :])
    #     plt.xlabel('x{}'.format(i))
    #     plt.ylabel('Counts')
    # # end for
    # plt.show()

    ##-------------------------
    # Cluster minima and associate each cluster with a local minimum.
    # Using a normalized coordinate space for cluster detection.

    x_range = bounds.ub - bounds.lb
    pts = np.array([x[0] for x in minima])
    pts_norm = (pts - bounds.lb)/x_range
    _, labels = dbscan(pts_norm, eps=0.05, min_samples=21, n_jobs=-1)

    # Visual presentation for adding to Jira ticket.
    plt.figure(figsize=(8, 8))
    plt.scatter(pts_burnin[:, 0], pts_burnin[:, 1], c='#202020', alpha=0.1, s=5)
    plt.scatter(pts[:, 0], pts[:, 1], s=5, c='k')
    for grp in range(max(labels) + 1):
        mask = (labels == grp)
        color = 'C' + str(grp)
        plt.scatter(pts[mask, 0], pts[mask, 1], s=5, c=color)
    # end for
    plt.xlim(bounds.lb[0], bounds.ub[0])
    plt.ylim(bounds.lb[1], bounds.ub[1])
    plt.axis('equal')
    plt.xlabel('x0')
    plt.ylabel('x1')
    plt.title(objective.__name__)
    plt.grid(color='#80808080', linestyle=':')
    plt.savefig(objective.__name__.replace('<', '').replace('>', '') + '.png', dpi=300)
    plt.show()

# end func


def bi_quadratic(x, mu, cov):
    # The number returned from the function must be non-negative.
    # The exponential of the negative of this value is the probablity.
    x2fac = np.sqrt(np.matmul(np.matmul((x - mu).T, cov), x - mu))
    return x2fac
# end func


def main():
    # DEV TESTING

    # Test functions as per https://en.wikipedia.org/wiki/Test_functions_for_optimization
    from landscapes.single_objective import sphere, himmelblau, easom, rosenbrock, rastrigin

    bounds = Bounds(np.array([-3, -3]), np.array([3, 3]))
    optimize_minimize_mhmcmc_cluster(sphere, bounds, burnin=10000, maxiter=50000)

    bounds = Bounds(np.array([-5, -5]), np.array([5, 5]))
    optimize_minimize_mhmcmc_cluster(himmelblau, bounds, burnin=10000, maxiter=50000, target_ar=0.3, T=10)

    bounds = Bounds(np.pi + np.array([-3, -2]), np.pi + np.array([3, 4]))
    optimize_minimize_mhmcmc_cluster(easom, bounds, burnin=10000, maxiter=50000, T=0.2)

    bounds = Bounds(np.array([-4, -4]), np.array([4, 4]))
    optimize_minimize_mhmcmc_cluster(rosenbrock, bounds, burnin=10000, maxiter=50000, T=1000)

    bounds = Bounds(np.array([-4, -4]), np.array([4, 4]))
    optimize_minimize_mhmcmc_cluster(lambda xy: rastrigin(xy) - 2, bounds, burnin=10000, maxiter=50000, T=5)


    # Custom test function
    mu = np.array([0, 1])
    cov = np.array([[5, -6.0], [-6.0, 20.0]])
    fixed_args = (mu, cov)
    bounds = Bounds(np.array([-3, -2]), np.array([3, 4]))
    optimize_minimize_mhmcmc_cluster(bi_quadratic, bounds, fixed_args, burnin=10000, maxiter=50000)
# end func


if __name__ == "__main__":
    main()
# end if
