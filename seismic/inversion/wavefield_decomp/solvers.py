#!/usr/bin/env python
# coding: utf-8
"""
Objective function minimization solvers.
"""

import copy
import time

import numpy as np
from tqdm.auto import tqdm
from scipy.optimize import OptimizeResult
from sortedcontainers import SortedList
from sklearn.cluster import dbscan

from seismic.inversion.wavefield_decomp.call_count_decorator import call_counter


# pylint: disable=invalid-name

DEFAULT_CLUSTER_EPS = 0.05


class SolverGlobalMhMcmc:
    """
    Drop-in custom solver for scipy.optimize.minimize, based on Metrolpolis-Hastings Monte Carlo
    Markov Chain random walk with burn-in and adaptive acceptance rate, followed by N-dimensional
    clustering.

    Rather than returning one global solution, the solver returns the N best ranked local solutions.
    It also returns a probability distribution of each unknown based on Monte Carlo statistics.

    """
    # TODO: Migrate free functions and classes into this solver class.
    #  But keep free function optimize_minimize_mhmcmc_cluster() for user convenience.
    pass
# end class


class HistogramIncremental():
    """Class to incrementally accumulate N-dimensional histogram stats at runtime.
    """
    def __init__(self, bounds, nbins=20):
        self._ndims = len(bounds.lb)
        assert len(bounds.ub) == self._ndims
        self._bins = np.linspace(bounds.lb, bounds.ub, nbins + 1).T
        self._hist = np.zeros((self._ndims, nbins), dtype=int)
    # end func

    def __iadd__(self, x):
        assert len(x) == self._ndims
        for i, _x in enumerate(x):
            idx = np.digitize(_x, self._bins[i, :])  # pylint: disable=unsubscriptable-object
            self._hist[i, idx - 1] += 1
        # end for
        return self
    # end func

    @property
    def dims(self):
        return self._ndims
    # end func

    @property
    def bins(self):
        if self.dims == 1:
            return self._bins.flatten()
        else:
            return self._bins.copy()
        # end if
    # end func

    @property
    def histograms(self):
        if self.dims == 1:
            return self._hist.flatten()
        else:
            return self._hist.copy()
        # end if
    # end func

# end class


class BoundedRandNStepper():
    """Step one dimensional at a time using normal distribution of steps within defined parameter bounds.
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
    def __init__(self, stepper, accept_rate=0.5, interval=50, factor=0.9, ar_tolerance=0.05):
        self.takestep = stepper
        self.target_accept_rate = accept_rate
        self.interval = interval
        self.factor = factor
        self.tolerance = ar_tolerance
        self.logger = None  # Assign directly to customize logger callable

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
        if scaled_step and self.logger:
            self.logger("adaptive stepsize: acceptance rate {} target {} new stepsize {} old stepsize {}"
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


def optimize_minimize_mhmcmc_cluster(objective, bounds, args=(), x0=None, T=1, N=3, burnin=100000, maxiter=1000000,
                                     target_ar=0.4, ar_tolerance=0.05,
                                     cluster_eps=DEFAULT_CLUSTER_EPS, rnd_seed=None,
                                     collect_samples=None, logger=None):
    """
    Minimize objective function and return up to N local minima solutions.

    :param objective: Objective function to minimize. Takes unpacked args as function call arguments and returns
        a float.
    :type objective: Callable(\*args) -> float
    :param bounds: Bounds of the parameter space.
    :type bounds: scipy.optimize.Bounds
    :param args: Any additional fixed parameters needed to completely specify the objective function.
    :type args: tuple or list
    :param x0: Initial guess. If None, will be selected randomly and uniformly within the parameter bounds.
    :type x0: numpy.array with same shape as elements of bounds
    :param T: The "temperature" parameter for the accept or reject criterion. To sample the domain well,
        should be in the order of the typical difference in local minima objective valuations.
    :type T: float
    :param N: Maximum number of minima to return
    :type N: int
    :param burnin: Number of random steps to discard before starting to accumulate statistics.
    :type burnin: int
    :param maxiter: Maximum number of steps to take (including burnin).
    :type maxiter: int
    :param target_ar: Target acceptance rate of point samples generated by stepping.
    :type target_ar: float between 0 and 1
    :param ar_tolerance: Tolerance on the acceptance rate before actively adapting the step size.
    :type ar_tolerance: float
    :param cluster_eps: Point proximity tolerance for DBSCAN clustering, in normalized bounds coordinates.
    :type cluster_eps: float
    :param rnd_seed: Random seed to force deterministic behaviour
    :type rnd_seed: int
    :param collect_samples: If not None and integral type, collect collect_samples at regular intervals
        and return as part of solution.
    :type collect_samples: int or NoneType
    :param logger: Logger instance for outputting log messages.
    :return: OptimizeResult containing solution(s) and solver data.
    :rtype: scipy.optimize.OptimizeResult with additional attributes
    """

    @call_counter
    def obj_counted(*args):
        return objective(*args)
    # end func

    assert maxiter >= 2*burnin, "maxiter {} should be at least twice burnin steps {}".format(maxiter, burnin)
    main_iter = maxiter - burnin

    if collect_samples is not None:
        assert isinstance(collect_samples, int), "collect_samples expected to be integral type"
        assert collect_samples > 0, "collect_samples expected to be positive"
    # end if

    beta = 1.0/T

    if rnd_seed is None:
        rnd_seed = int(time.time()*1000) % (1 << 31)
    # end if
    np.random.seed(rnd_seed)
    if logger:
        logger.info('Using random seed {}'.format(rnd_seed))
    # end

    if x0 is None:
        x0 = np.random.uniform(bounds.lb, bounds.ub)
    # end if
    assert np.all((x0 >= bounds.lb) & (x0 <= bounds.ub))
    x = x0.copy()
    funval = obj_counted(x, *args)

    # Set up stepper with adaptive acceptance rate
    stepper = BoundedRandNStepper(bounds)
    stepper = AdaptiveStepsize(stepper, accept_rate=target_ar, ar_tolerance=ar_tolerance, interval=50)

    # -------------------------------
    # DO BURN-IN
    rejected_randomly = 0
    accepted_burnin = 0
    tracked_range = tqdm(range(burnin), total=burnin, desc='BURN-IN')
    if logger:
        stepper.logger = lambda msg: tracked_range.write(logger.name + ':' + msg)
    else:
        stepper.logger = tracked_range.write
    # end if
    for _ in tracked_range:
        x_new = stepper(x)
        funval_new = obj_counted(x_new, *args)
        log_alpha = -(funval_new - funval)*beta
        if log_alpha > 0 or np.log(np.random.rand()) <= log_alpha:
            x = x_new
            funval = funval_new
            stepper.notify_accept()
            accepted_burnin += 1
        elif log_alpha <= 0:
            rejected_randomly += 1
        # end if
    # end for
    ar = float(accepted_burnin)/burnin
    if logger:
        logger.info("Burn-in acceptance rate: {}".format(ar))
    # end if

    # -------------------------------
    # DO MAIN LOOP
    if collect_samples is not None:
        nsamples = min(collect_samples, main_iter)
        sample_cadence = main_iter/nsamples
        samples = np.zeros((nsamples, len(x)))
        samples_fval = np.zeros(nsamples)
    # end if
    accepted = 0
    rejected_randomly = 0
    minima_sorted = SortedList(key=lambda rec: rec[1])  # Sort by objective function value
    hist = HistogramIncremental(bounds, nbins=100)
    # Cached a lot of potential minimum values, as these need to be clustered before return N results
    N_cached = int(np.ceil(N*main_iter/500))
    next_sample = 0.0
    sample_count = 0
    tracked_range = tqdm(range(main_iter), total=main_iter, desc='MAIN')
    if logger:
        stepper.logger = lambda msg: tracked_range.write(logger.name + ':' + msg)
    else:
        stepper.logger = tracked_range.write
    # end if
    for i in tracked_range:
        if collect_samples and i >= next_sample:
            assert sample_count < collect_samples
            samples[sample_count] = x
            samples_fval[sample_count] = funval
            sample_count += 1
            next_sample += sample_cadence
        # end if
        x_new = stepper(x)
        funval_new = obj_counted(x_new, *args)
        log_alpha = -(funval_new - funval)*beta
        if log_alpha > 0 or np.log(np.random.rand()) <= log_alpha:
            x = x_new
            funval = funval_new
            minima_sorted.add((x, funval))
            if len(minima_sorted) > N_cached:
                minima_sorted.pop()
            # end if
            stepper.notify_accept()
            hist += x
            accepted += 1
        elif log_alpha <= 0:
            rejected_randomly += 1
        # end if
    # end for
    stepper.logger = None
    ar = float(accepted)/main_iter
    if logger:
        logger.info("Acceptance rate: {}".format(ar))
        logger.info("Best minima (before clustering):\n{}".format(np.array([_mx[0] for _mx in minima_sorted[:10]])))
    # end if

    # -------------------------------
    # Cluster minima and associate each cluster with a local minimum.
    # Using a normalized coordinate space for cluster detection.
    x_range = bounds.ub - bounds.lb
    pts = np.array([x[0] for x in minima_sorted])
    fvals = np.array([x[1] for x in minima_sorted])
    pts_norm = (pts - bounds.lb)/x_range
    _, labels = dbscan(pts_norm, eps=cluster_eps, min_samples=21, n_jobs=-1)

    # Compute mean of each cluster and evaluate objective function at cluster mean locations.
    minima_candidates = []
    for grp in range(max(labels) + 1):
        mask = (labels == grp)
        mean_loc = np.mean(pts[mask, :], axis=0)
        # Evaluate objective function precisely at the mean location of each cluster
        fval = obj_counted(mean_loc, *args)
        minima_candidates.append((mean_loc, grp, fval))
    # end for

    # Rank minima locations by objective function.
    minima_candidates.sort(key=lambda c: c[2])

    # Pick up to N solutions
    solutions = minima_candidates[:N]

    # Put results into OptimizeResult container.
    # Add histograms to output result (in form of scipy.stats.rv_histogram)
    solution = OptimizeResult()
    solution.x = np.array([s[0] for s in solutions])
    solution.clusters = [pts[(labels == s[1])] for s in solutions]
    solution.cluster_funvals = [fvals[(labels == s[1])] for s in solutions]
    solution.bins = hist.bins
    solution.distribution = hist.histograms
    solution.acceptance_rate = ar
    solution.success = True
    solution.status = 0
    if len(solutions) > 0:
        solution.message = 'SUCCESS: Found {} local minima'.format(len(solutions))
    else:
        solution.message = 'WARNING: Found no clusters within tolerance {}'.format(cluster_eps)
    # end if
    solution.fun = np.array([s[2] for s in solutions])
    solution.jac = None
    solution.nfev = obj_counted.counter
    solution.njev = 0
    solution.nit = main_iter
    solution.maxcv = None
    solution.samples = samples if collect_samples else None
    solution.sample_funvals = samples_fval if collect_samples else None
    solution.bounds = bounds
    solution.version = 's0.3'  # Solution version for future traceability
    solution.rnd_seed = rnd_seed

    return solution

# end func
