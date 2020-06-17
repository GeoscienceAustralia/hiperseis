#!/usr/bin/env python

import logging

import numpy as np
from scipy.interpolate import interp1d

import seismic.receiver_fn.rf_util as rf_util
from seismic.units_utils import KM_PER_DEG

# pylint: disable=invalid-name,logging-format-interpolation

logging.basicConfig()


def compute_hk_stack(cha_data, V_p=None, h_range=np.linspace(20.0, 70.0, 251),
                     k_range=np.linspace(1.4, 2.0, 301), root_order=1, include_t3=True):
    """Compute H-k stacking array on a dataset of receiver functions.

    :param cha_data: List or iterable of RF traces to use for H-k stacking.
    :type cha_data: iterable(rf.RFTrace)
    :param V_p: P-wave velocity in crustal layer, defaults to None in which case it is inferred from trace metadata
    :type V_p: float, optional
    :param h_range: Range of h values (Moho depth) values to cover, defaults to np.linspace(20.0, 70.0, 251)
    :type h_range: numpy.array [1D], optional
    :param k_range: Range of k values to cover, defaults to np.linspace(1.4, 2.0, 301)
    :type k_range: numpy.array [1D], optional
    :param root_order: Exponent for nth root stacking as per K.J.Muirhead (1968), defaults to 1
    :type root_order: int, optional
    :param include_t3: If True, include the t3 (PpSs+PsPs) multiple in the stacking, defaults to True
    :type include_t3: bool, optional
    :return: k-grid values [2D], h-grid values [2D], H-k stack in series of 2D layers having one layer per multiple
    :rtype: numpy.array [2D], numpy.array [2D], numpy.array [3D]
    """
    # It is *critical* for the correctness of the output of this function that the V_p passed in is the same as the
    # V_p of the shallowest layer of the velocity model (from tau-py model) that was used to compute the arrival
    # inclination of each trace. Leave as None to have this function infer V_p from the trace metadata.

    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)

    infer_Vp = (V_p is None)
    V_p_inferred = infer_Vp_from_traces(cha_data, log)
    if infer_Vp:
        V_p = V_p_inferred
        log.debug("Inferred V_p = {}".format(V_p))
    else:
        log.debug("Using V_p = {} (inferred V_p from shallowest layer = {})".format(V_p, V_p_inferred))
    # end if

    # Pre-compute grid quantities
    k_grid, h_grid = np.meshgrid(k_range, h_range)
    # hk_stack = np.zeros_like(k_grid)

    # Whether to use RF slowness as the ray parameter (after unit conversion) and as source of V_p,
    # or else use specific V_p given by user.
    stream_stack = []
    # Loop over traces, compute times, and stack interpolated values at those times
    channel_id = None
    for tr in cha_data:
        if channel_id is None:
            channel_id = tr.stats.channel
        else:
            assert tr.stats.channel == channel_id, \
                "Stacking mismatching channel data: expected {}, found {}".format(channel_id, tr.stats.channel)
        # end if

        t1, t2, t3 = compute_theoretical_phase_times(tr, h_grid, k_grid, V_p, include_t3=include_t3)

        # Apply custom time offsets if given.
        for tval, offset_name in ((t1, 't1_offset'), (t2, 't2_offset'), (t3, 't3_offset')):
            try:
                if tval is not None:
                    tval += tr.stats.sediment[offset_name]
            except AttributeError:
                pass
            # end try
        # end for

        # Subtract lead time so that primary P-wave arrival is at time zero.
        lead_time = tr.stats.onset - tr.stats.starttime
        times = tr.times() - lead_time
        # Create interpolator from stream signal for accurate time sampling.
        interpolator = interp1d(times, tr.data, kind='linear', copy=False, bounds_error=False, assume_sorted=True)

        phase_sum = []
        phase_sum.append(rf_util.signed_nth_root(interpolator(t1), root_order))
        phase_sum.append(rf_util.signed_nth_root(interpolator(t2), root_order))
        if include_t3:
            # Negative sign on the third term is intentional, see Chen et al. (2010) and Zhu & Kanamori (2000).
            # It needs to be negative because the PpSs + PsPs peak has negative phase,
            # see http://eqseis.geosc.psu.edu/~cammon/HTML/RftnDocs/rftn01.html
            # Apply nth root technique to reduce uncorrelated noise (Chen et al. (2010))
            phase_sum.append(-rf_util.signed_nth_root(interpolator(t3), root_order))
        # end if

        stream_stack.append(phase_sum)
    # end for

    # Perform the stacking (sum) across streams. hk_stack retains separate t1, t2, and t3 components here.
    hk_stack = np.nanmean(np.array(stream_stack), axis=0)

    # This inversion of the nth root is different to Sippl and Chen, but consistent with Muirhead
    # who proposed the nth root technique. It improves the contrast of the resulting plot.
    if root_order != 1:
        hk_stack = rf_util.signed_nth_power(hk_stack, root_order)

    return k_grid, h_grid, hk_stack
# end func


def compute_weighted_stack(hk_components, weighting=(0.5, 0.5, 0.0)):
    """Given stack components from function `compute_hk_stack`, compute the overall weighted stack.

    :param hk_components: H-k stack layers returned from `compute_hk_stack`
    :type hk_components: np.array
    :param weighting: Weightings for (t1, t2, t3) layers respectively, defaults to (0.5, 0.5, 0.0)
    :type weighting: tuple, optional
    :return: Weighted stack in H-k space
    :rtype: numpy.array
    """
    assert hk_components.shape[0] == len(weighting), hk_components.shape
    hk_phase_stacked = np.dot(np.moveaxis(hk_components, 0, -1), np.array(weighting))
    return hk_phase_stacked
# end func


def infer_Vp_from_traces(cha_data, log=None):
    # Determine the internal V_p consistent with the trace ray parameters and inclinations.
    V_p_values = []
    for tr in cha_data:
        p = tr.stats.slowness / KM_PER_DEG
        incl_deg = tr.stats.inclination
        incl_rad = np.deg2rad(incl_deg)
        V_p_value = np.sin(incl_rad) / p
        V_p_values.append(V_p_value)
    # end for
    V_p_values = np.array(V_p_values)
    if  log is not None and not np.allclose(V_p_values, V_p_values, rtol=1e-3, atol=1e-4):
        log.error("Inconsistent V_p values inferred from traces, H-k stacking results may be unreliable!")
    # end if
    V_p = np.mean(V_p_values)

    return V_p
# end func


def compute_theoretical_phase_times(tr, H, k, V_p, include_t3=True):
    incl_deg = tr.stats.inclination
    incl_rad = np.deg2rad(incl_deg)
    sin_i = np.sin(incl_rad)
    p = sin_i / V_p

    H_on_V_p = H/V_p
    k2 = k*k

    p2_Vp2 = p * p * V_p * V_p
    term1 = H_on_V_p * np.sqrt(k2 - p2_Vp2)
    term2 = H_on_V_p * np.sqrt(1 - p2_Vp2)

    # Time for Ps
    t1 = term1 - term2

    # Time for PpPs
    t2 = term1 + term2

    if include_t3:
        # Time for PpSs + PsPs
        t3 = 2 * term1
    else:
        t3 = None
    # end if

    return t1, t2, t3
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


def find_local_hk_maxima(k_grid, h_grid, hk_stack, min_rel_value=0.5):
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
