#!/usr/bin/env python

import logging

import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import correlate
import seismic.receiver_fn.rf_util as rf_util
from seismic.units_utils import KM_PER_DEG
from rf.util import DEG2KM
from scipy.ndimage import gaussian_filter
from scipy import interpolate
from sklearn.cluster import dbscan
from scipy.integrate import simps as simpson
from scipy.optimize import minimize
from scipy.optimize import Bounds
import obspy
from obspy.taup.velocity_model import VelocityModel
import os

# pylint: disable=invalid-name,logging-format-interpolation

logging.basicConfig()

DEFAULT_H_RANGE = tuple(np.linspace(20.0, 70.0, 501))
DEFAULT_k_RANGE = tuple(np.linspace(1.5, 2.0, 301))
DEFAULT_WEIGHTS = np.array([0.5, 0.4, 0.1])

DEFAULT_SED_H_RANGE = tuple(np.linspace(0.01, 6, 21))
DEFAULT_SED_k_RANGE = tuple(np.linspace(1.0, 5.0, 21))

def compute_hk_stack(cha_data, h_range=None, k_range=None,
                     weights=DEFAULT_WEIGHTS, root_order=1):
    """Compute H-k stacking array on a dataset of receiver functions.

    :param cha_data: List or iterable of RF traces to use for H-k stacking.
    :type cha_data: Iterable(rf.RFTrace)
    :param h_range: Range of h values (Moho depth) values to cover, defaults to np.linspace(20.0, 70.0, 251)
    :type h_range: numpy.array [1D], optional
    :param k_range: Range of k values to cover, defaults to np.linspace(1.4, 2.0, 301)
    :type k_range: numpy.array [1D], optional
    :param weights: numpy array of length 3, containing weights for the three phases (Ps, PpPs and (PpSs + PsPs))
    :type weights: numpy.array [1D], optional
    :param root_order: Exponent for nth root stacking as per K.J.Muirhead (1968), defaults to 1
    :type root_order: int, optional
    :return: k-grid values [2D], h-grid values [2D], H-k stack values [2D]
    :rtype: numpy.array [2D], numpy.array [2D], numpy.array [2D]
    """
    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)

    if h_range is None:
        h_range = DEFAULT_H_RANGE
    # end if
    if k_range is None:
        k_range = DEFAULT_k_RANGE
    # end if

    # Pre-compute grid quantities
    k_grid, h_grid = np.meshgrid(k_range, h_range)

    tphase_amps = []
    for itrc, trc in enumerate(cha_data):
        lead_time = trc.stats.onset - trc.stats.starttime
        p = trc.stats.slowness / DEG2KM
        incl_deg = trc.stats.inclination
        incl_rad = np.deg2rad(incl_deg)
        Vp_inv = p / np.sin(incl_rad)
        Vs_inv = k_grid * Vp_inv

        term1 = np.sqrt(Vs_inv ** 2 - p ** 2)
        term2 = np.sqrt(Vp_inv ** 2 - p ** 2)

        t1 = h_grid * (term1 - term2)
        t2 = h_grid * (term1 + term2)
        t3 = h_grid * 2 * term1

        try:
            t1 += trc.stats.t1_offset
            t2 += trc.stats.t2_offset
            t3 += trc.stats.t3_offset
        except:
            pass
        # end try

        times = trc.times() - lead_time
        times_min = np.min(times)
        times_max = np.max(times)
        if(np.min(t1) < times_min or \
           np.max(t1) > times_max or \
           np.min(t2) < times_min or \
           np.max(t2) > times_max or \
           np.min(t3) < times_min or \
           np.max(t3) > times_max):

            nsl = '.'.join([trc.stats.network, trc.stats.station, trc.stats.location])
            log.warning('\nCorrected times for a trace in {} fall outside the available time-range'.format(nsl))
        # end if

        tio = interp1d(times, trc.data, fill_value=0, bounds_error=False)

        a, b, c = tio(t1), tio(t2), -tio(t3)
        tphase_amps.append([np.sign(a) * np.power(np.fabs(a), 1. / root_order),
                            np.sign(b) * np.power(np.fabs(b), 1. / root_order),
                            np.sign(c) * np.power(np.fabs(c), 1. / root_order)])
    # end for
    tphase_amps = np.array(tphase_amps)
    hk_stack = np.sum(np.dot(np.moveaxis(tphase_amps, 1, -1), weights), axis=0)
    hk_stack = np.sign(hk_stack) * np.power(np.fabs(hk_stack), root_order)

    return k_grid, h_grid, hk_stack
# end func

def compute_sediment_hk_stack(cha_data, H_c, k_c, h_range=None, k_range=None, root_order=9):
    """Compute H-k stacking array on a dataset of receiver functions.

    :param cha_data: List or iterable of RF traces to use for H-k stacking.
    :type cha_data: Iterable(rf.RFTrace)
    :param h_range: Range of h values (Moho depth) values to cover, defaults to np.linspace(20.0, 70.0, 251)
    :type h_range: numpy.array [1D], optional
    :param k_range: Range of k values to cover, defaults to np.linspace(1.4, 2.0, 301)
    :type k_range: numpy.array [1D], optional
    :param root_order: Exponent for nth root stacking as per K.J.Muirhead (1968), defaults to 1
    :type root_order: int, optional
    :return: k-grid values [2D], h-grid values [2D], H-k stack values [2D]
    :rtype: numpy.array [2D], numpy.array [2D], numpy.array [2D]
    """
    def obj_func(x0, amps):
        curr_stack = np.sum(np.dot(np.moveaxis(amps, 1, -1), x0), axis=0)
        idx = np.unravel_index(np.argmax(curr_stack), np.array(curr_stack).shape)
        return -curr_stack[idx]
    # end func

    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)

    if h_range is None:
        h_range = DEFAULT_SED_H_RANGE
    # end if
    if k_range is None:
        k_range = DEFAULT_SED_k_RANGE
    # end if

    iasp91_path = os.path.join(os.path.abspath(os.path.dirname(obspy.__file__)), 'taup/data/iasp91.tvel')
    vmodel = VelocityModel.read_tvel_file(iasp91_path)
    d = np.linspace(-80, -1, 120)
    v = vmodel.evaluate_above(depth=-d, prop='P')
    d[-1] = 0
    v[-1:-7:-1] = np.linspace(2, v[-6], 6)
    vfunc = interp1d(d, v)

    # Pre-compute grid quantities
    k_grid, h_grid = np.meshgrid(k_range, h_range)

    tphase_amps = []
    for itrc, trc in enumerate(cha_data):
        lead_time = trc.stats.onset - trc.stats.starttime
        p = trc.stats.slowness / DEG2KM
        incl_deg = trc.stats.inclination
        incl_rad = np.deg2rad(incl_deg)
        Vp_inv = p / np.sin(incl_rad)

        t4 = np.zeros(h_grid.shape)
        t2 = np.zeros(h_grid.shape)
        t3 = np.zeros(h_grid.shape)

        for i in np.arange(h_grid.shape[0]):
            interval1 = np.linspace(-h_grid[i, 0], 0, 50)
            interval2 = np.linspace(-(h_grid[i, 0] + H_c), -h_grid[i, 0], 50)
            v_interval1 = vfunc(interval1)
            v_interval2 = vfunc(interval2)

            for j in np.arange(h_grid.shape[1]):
                A = simpson(np.sqrt(np.power(v_interval1 / k_grid[i, j], -2.) - p ** 2), interval1)
                B = simpson(np.sqrt(np.power(v_interval1, -2.) - p ** 2), interval1)
                C = simpson(np.sqrt(np.power(v_interval2 / k_c, -2.) - p ** 2), interval2)
                D = simpson(np.sqrt(np.power(v_interval2, -2.) - p ** 2), interval2)

                t4[i, j] = A - B
                t2[i, j] = A + B + C + D
                t3[i, j] = 2 * A + 2 * C
            # end for
        # end for

        tio = interp1d(trc.times() - lead_time, trc.data)

        a, b, c = tio(t4), tio(t2), -tio(t3)
        tphase_amps.append([np.sign(a) * np.power(np.fabs(a), 1. / root_order),
                            np.sign(b) * np.power(np.fabs(b), 1. / root_order),
                            np.sign(c) * np.power(np.fabs(c), 1. / root_order)])
    # end for
    tphase_amps = np.array(tphase_amps)

    bounds = Bounds(np.zeros(3) + 0.01, np.array([0.3, 1, 1]))
    constraints = [{'type': 'ineq', 'fun': lambda x: x},
                   {'type': 'eq', 'fun': lambda x: np.sum(x) - 1}]

    starting_weights = np.array([0.2, 0.4, 0.4])
    opt_result = minimize(obj_func, starting_weights, tphase_amps, method='SLSQP',
                          bounds=bounds, constraints=constraints)
    x = opt_result['x']
    weights = x

    #print('Weights: ', weights)
    hk_stack = np.sum(np.dot(np.moveaxis(tphase_amps, 1, -1), weights), axis=0)
    hk_stack = np.sign(hk_stack) * np.power(np.fabs(hk_stack), root_order)

    return k_grid, h_grid, hk_stack, weights
# end func


def find_global_hk_maximum(k_grid, h_grid, hk_weighted_stack):
    """Given the weighted stack computed from function `compute_weighted_stack` and the corresponding
    k-grid and h-grid, find the location in H-k space of the global maximum.

    :param k_grid: Grid of k-values
    :type k_grid: Two-dimensional numpy.array
    :param h_grid: Grid of H-values
    :type h_grid: Two-dimensional numpy.array
    :param hk_weighted_stack: Grid of stacked RF sample values produced by function
        rf_stacking.computed_weighted_stack()
    :type hk_weighted_stack: Two-dimensional numpy.array
    :return: Location of global maximum on the H-k grid of the maximum stack value.
    :rtype: tuple(float, float)
    """
    max_loc = np.unravel_index(np.argmax(hk_weighted_stack), hk_weighted_stack.shape)
    h_max = h_grid[max_loc[0], 0]
    k_max = k_grid[0, max_loc[1]]
    return (h_max, k_max)
# end func


def _find_local_hk_maxima_helper(k_grid, h_grid, hk_stack, min_rel_value=0.5):
    """Given the weighted stack computed from function `compute_weighted_stack` and the corresponding
    k-grid and h-grid, find the locations in H-k space of all local maxima above a certain threshold.

    :param k_grid: Grid of k-values
    :type k_grid: Two-dimensional numpy.array
    :param h_grid: Grid of H-values
    :type h_grid: Two-dimensional numpy.array
    :param hk_stack: Grid of stacked RF sample values produced by function
        rf_stacking.computed_weighted_stack()
    :type hk_stack: Two-dimensional numpy.array
    :param min_rel_value: Minimum value required relative to the largest value in the stack, defaults to 0.5
    :type min_rel_value: float, optional
    :return: List of tuples containing parameters of local maxima solutions, with values in the following
        order: (H, k, stack_value, row_index, col_index)
    :rtype: list(tuple(float, float, float, int, int))
    """
    # Determine global maximum, as we only keep local maxima that are at least min_rel_value
    # proportion of this value.
    global_max = np.nanmax(hk_stack)
    # Compute 2D mask of all values that are greater than or equal to their 4 neighbours.
    m = ((hk_stack[1:-1, 1:-1] > hk_stack[1:-1, 0:-2]) & (hk_stack[1:-1, 1:-1] > hk_stack[1:-1, 2:]) &
         (hk_stack[1:-1, 1:-1] > hk_stack[0:-2, 1:-1]) & (hk_stack[1:-1, 1:-1] > hk_stack[2:, 1:-1]) &
         (hk_stack[1:-1, 1:-1] >= min_rel_value*global_max))
    # Find the row and column indices where m is True
    m_idx = np.nonzero(m)
    # Determine the stack values at the identified local maxima
    stack_vals = hk_stack[1:-1, 1:-1][m_idx]
    # Determine the k-values at the identified local maxima
    k_vals = k_grid[1:-1, 1:-1][m_idx]
    # Determine the h-values at the identified local maxima
    h_vals = h_grid[1:-1, 1:-1][m_idx]
    # Zip the candidate solutions into a tuple (H, k, stack, row_index, col_index).
    # Note that the row and column index here are in the original k- and h-grid, so must have +1 added
    # since the masking was done on the interior grid points only.
    solutions = tuple(zip(h_vals, k_vals, stack_vals, m_idx[0] + 1, m_idx[1] + 1))
    # Sort the solutions from highest stack value to the lowest.
    solutions = sorted(solutions, key=lambda v: v[2], reverse=True)
    return solutions
# end func

def find_local_hk_maxima(k_grid, h_grid, hk_stack_sum, max_number=3):
    # Method here is:
    # 1) find all local maxima
    # 2) cluster local maxima and compute centroid of each cluster
    # The centroids of the top max_number clusters are returned, ranked by stack amplitude.

    # Only consider positive stack regions.
    hk_stack = hk_stack_sum.copy()
    hk_stack[hk_stack < 0] = 0
    hk_stack_max = np.nanmax(hk_stack)

    # Smooth the stack, as we're not interested in high frequency local maxima
    hk_stack = gaussian_filter(hk_stack, sigma=3, mode='nearest')

    # This method only returns locations in the interior, not on the boundary of the domain
    local_maxima = _find_local_hk_maxima_helper(k_grid, h_grid, hk_stack, min_rel_value=0.7)
    if len(local_maxima) <= 1:
        return [(h, k) for h, k, _, _, _ in local_maxima]
    # end if

    # Perform clustering in normalized coordinates
    k_min, k_max = (np.nanmin(k_grid), np.nanmax(k_grid))
    k_range = k_max - k_min
    h_min, h_max = (np.nanmin(h_grid), np.nanmax(h_grid))
    h_range = h_max - h_min

    # Use DBSCAN to cluster nearby pointwise local maxima
    eps = 0.05
    pts_norm = np.array([[(k - k_min)/k_range, (h - h_min)/h_range, v/hk_stack_max] for h, k, v, _, _ in local_maxima])
    pts_hk = np.array([[h, k] for h, k, _, _, _ in local_maxima])
    _, labels = dbscan(pts_norm, eps, min_samples=2, metric='euclidean')

    # Collect group-based local maxima
    maxima_coords = []
    group_ids = set(labels[labels >= 0])
    for grp_id in group_ids:
        maxima_coords.append(np.mean(pts_hk[labels == grp_id], axis=0))
    # end for

    # Collect remaining non-grouped points and add them to list of local maxima
    loners = pts_hk[labels < 0]
    if np.any(loners):
        maxima_coords.extend(loners)
    # end if

    # Sort the maxima by amplitude of stack, and then discard values with amplitude too weak compared to strongest peak
    if len(maxima_coords) > 1:
        finterp = interpolate.interp2d(k_grid[0, :], h_grid[:, 0], hk_stack_sum)
        maxima_coords.sort(key=lambda p: finterp(p[1], p[0]), reverse=True)
        strongest = finterp(maxima_coords[0][1], maxima_coords[0][0])
        while finterp(maxima_coords[-1][1], maxima_coords[-1][0]) < 0.8*strongest:
            maxima_coords.pop()
        # end while
    # end if

    return maxima_coords[:max_number]
# end func
