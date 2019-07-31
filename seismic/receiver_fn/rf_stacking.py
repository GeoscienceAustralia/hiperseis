#!/usr/bin/env python

import logging

import numpy as np
from scipy.interpolate import interp1d

import seismic.receiver_fn.rf_util as rf_util

# pylint: disable=invalid-name

logging.basicConfig()


def compute_hk_stack(db_station, cha, h_range=np.linspace(10.0, 70.0, 301), k_range = np.linspace(1.3, 2.1, 201),
                     V_p=6.4, root_order=1, include_t3=True):

    # Pre-compute grid quantities
    k_grid, h_grid = np.meshgrid(k_range, h_range)
    hk_stack = np.zeros_like(k_grid)
    H_on_V_p = h_grid/V_p
    k2 = k_grid*k_grid

    stream_stack = []
    cha_data = db_station[cha]
    # Loop over traces, compute times, and stack interpolated values at those times
    for tr in cha_data:
        incl = tr.stats.inclination
        incl_rad = incl*np.pi/180.0
        sin_i = np.sin(incl_rad)
        slowness_sec_per_km = tr.stats.slowness/rf_util.KM_PER_DEG
        # p is the ray parameter
        p = sin_i*slowness_sec_per_km
        p2_Vp2 = p*p*V_p*V_p
        term1 = H_on_V_p*np.sqrt(k2 - p2_Vp2)
        term2 = H_on_V_p*np.sqrt(1 - p2_Vp2)
        # Time for Ps
        t1 = term1 - term2
        # Time for PpPs
        t2 = term1 + term2
        if include_t3:
            # Time for PpSs + PsPs
            t3 = 2*term1

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

        stream_stack.append(phase_sum)
    # end for

    # Perform the stacking (sum) across streams. hk_stack retains separate t1, t2, and t3 components here.
    hk_stack = np.nanmean(np.array(stream_stack), axis=0)

    # This inversion of the nth root is different to Sippl and Chen, but consistent with Muirhead
    # who proposed the nth root technique. It improves the contrast of the resulting plot.
    if root_order != 1:
        hk_stack = rf_util.signed_nth_power(hk_stack, root_order)

    return k_grid, h_grid, hk_stack


def compute_weighted_stack(hk_components, weighting=(0.5, 0.5, 0.0)):
    assert hk_components.shape[0] == len(weighting), hk_components.shape
    hk_phase_stacked = np.dot(np.moveaxis(hk_components, 0, -1), np.array(weighting))
    return hk_phase_stacked
