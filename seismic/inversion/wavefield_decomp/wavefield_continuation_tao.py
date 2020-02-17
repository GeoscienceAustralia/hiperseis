#!/usr/bin/env python
# coding: utf-8
"""Class encapsulating algorithm for 1D inversion using wavefield continuation.

    Based on reference:
    Kai Tao, Tianze Liu, Jieyuan Ning, Fenglin Niu, "Estimating sedimentary and crustal structure
    using wavefield continuation: theory, techniques and applications", *Geophysical Journal International*,
    Volume 197, Issue 1, April, 2014, Pages 443-457, https://doi.org/10.1093/gji/ggt515
"""

import numpy as np

from seismic.receiver_fn.rf_util import KM_PER_DEG
from seismic.receiver_fn.rf_util import sinc_resampling


class WfContinuationSuFluxComputer:
    """
    Implements computation of the upwards mean S-wave energy flux at the top of the mantle
    for an ensemble of events for one station.

    Class instance can be loaded with a dataset and then evaluated for arbitrary 1D earth models.

    Process:
    1. Load data from a station event dataset.  Copy of the data is buffered by this
        class in efficient format for energy flux calculation.
    2. Define 1D earth model and mantle half-space material properties (in external code).
    3. Call instance with models and receive energy flux results.
    """
    def __init__(self, station_event_dataset, f_s, time_window, cut_window):
        """

        :param station_event_dataset: Iterable container of obspy.Stream objects.
        :param f_s: Processing sample rate. Usuually much less than the sampling rate of the input raw seismic traces.
        :param time_window:
        :param cut_window:
        """

        if not station_event_dataset:
            return
        # end if

        # Check input streams are all in ZRT coordinates.
        assert np.all([''.join([tr.stats.channel[-1] for tr in st]).upper() == 'ZRT' for st in station_event_dataset])

        # Check input streams all have same sampling rate
        traces_fs = set(tr.stats.sampling_rate for st in station_event_dataset for tr in st)
        assert len(traces_fs) == 1, 'Inconsistent sampling rates {}!'.format(traces_fs)

        self._time_window = time_window
        self._cut_window = cut_window
        # Sampling rate and time step
        self._f_s = f_s
        self._dt = 1.0 / f_s
        # Number of sample points. Since the last point is considered the end of a closed interval,
        # we +1 to the numer of points to ensure f_s == 1/dt.
        self._npts = int(np.rint((time_window[1] - time_window[0])*f_s)) + 1
        self._time_axis = np.linspace(*self._time_window, self._npts)
        assert np.isclose(np.mean(np.diff(self._time_axis)), self._dt)
        self._station_eventdataset_to_v0(station_event_dataset)
        self._station_event_dataset_extract_p(station_event_dataset)
        self._precompute()
    # end func

    def _station_eventdataset_to_v0(self, data):
        """
        Convert dict of streams (indexed by event id) to numpy array in format required
        by function compute_su_energy(), including resampling to f_s and applying a cut
        window and sinc resampling.

        This conversion may be expensive and compute_su_energy() may need to be called
        many times, so preconverting the format to numpy array once only is important
        to overall performance.

        Any quality filtering needs to be performed prior to calling this function.

        Note: This function modifies in-place the traces in the values of data.

        :assigns self._v0: Numpy array of shape (N_events, 2, N_samples) containing the
            R- and Z-component traces for all events at sample rate f_s and covering
            duration of self._time_window.
        """

        # # Resample to f_s if any trace is not already at f_s
        for stream in data:
            if np.any(np.array([tr.stats.sampling_rate != self._f_s for tr in stream])):
                # Resampling lowpass only, as per Tao (anti-aliasing).
                # In order to not mutate input data, we filter a copy and replace the stream traces
                # with the filtered traces.
                stream_filt = stream.copy().filter('lowpass', freq=self._f_s/2.0, corners=2, zerophase=True)\
                                           .interpolate(self._f_s, method='lanczos', a=10)
                stream.clear()
                stream += stream_filt
            # end if
        # end for

        # Trim to time window
        for stream in data:
            stream.trim(stream[0].stats.onset + self._time_window[0],
                        stream[0].stats.onset + self._time_window[1])
        # end for

        # Cut central data segment and resample back to original length using sinc interpolation.
        times = self._time_axis
        for stream in data:
            for tr in stream:
                tr_cut = tr.copy().trim(tr.stats.onset + self._cut_window[0], tr.stats.onset + self._cut_window[1])
                tr_cut.detrend('linear')
                tr_cut.taper(0.10)
                cut_times = tr_cut.times() - (tr_cut.stats.onset - tr_cut.stats.starttime)
                cut_data = tr_cut.data
                resampled_data = sinc_resampling(cut_times, cut_data, times)
                # Replace trace data with cut resampled data
                tr.data = resampled_data
            # end for
        # end for

        # Check that all traces have same length
        trace_lengths = set(len(st[0]) for st in data)
        assert len(trace_lengths) == 1, 'Inconsistent trace lengths {} after sinc interpolation!'.format(trace_lengths)

        # Pull data arrays out into matrix format
        self._v0 = np.array([[st.select(component='R')[0].data.tolist(),
                              (-st.select(component='Z')[0].data).tolist()] for st in data])

    # end func

    def _station_event_dataset_extract_p(self, data):
        """
        Extract ray parameters from input data.

        :param data:
        :type data:
        :return: None
        """
        self._p = np.array([stream[0].stats.slowness/KM_PER_DEG for stream in data])
    # end func

    def _precompute(self):
        # Computations that need only be done once per dataset.
        self._nevts = self._v0.shape[0]  # Number of events

        # *NORMALIZE v0*
        # Reshape to facilitate max_vz normalization using numpy broadcast rules.
        v0 = np.moveaxis(self._v0, 0, -1)
        # Normalize each event signal by the maximum z-component amplitude.
        # We perform this succinctly using numpy multidimensional broadcasting rules.
        max_vz = np.abs(v0[1, :, :]).max(axis=0)
        v0 = v0 / max_vz
        # Reshape back to original shape.
        self._v0 = np.moveaxis(v0, -1, 0)

        # Transform v0 to the spectral domain using real FFT
        self._fv0 = np.fft.fft(self._v0, axis=-1)

        # Compute discrete frequencies
        self._w = 2 * np.pi * np.fft.fftfreq(self._npts, self._dt)

    # end if

    def __call__(self, mantle_props, layer_props, flux_window=(-10, 20)):
        """Compute upgoing S-wave energy for a given set of seismic time series v0.

        :param mantle_props: LayerProps representing mantle properties.
        :param layer_props: List of LayerProps.
        """
        # This is the callable operator that performs computations of energy flux

        # Compute mode matrices for mantle
        M_m, Minv_m, _ = WfContinuationSuFluxComputer._mode_matrices(mantle_props.Vp, mantle_props.Vs, mantle_props.rho,
                                                                     self._p)

        # Propagate from surface
        fvm = WfContinuationSuFluxComputer._propagate_layers(self._fv0, self._w, layer_props, self._p)
        fvm = np.matmul(Minv_m, fvm)

        num_pos_freq_terms = (fvm.shape[2] + 1) // 2
        # Velocities at top of mantle
        vm = np.fft.irfft(fvm[:, :, :num_pos_freq_terms], self._npts, axis=2)

        # Compute coefficients of energy integral for upgoing S-wave
        qb_m = np.sqrt(1 / mantle_props.Vs ** 2 - self._p * self._p)
        Nsu = self._dt * mantle_props.rho * (mantle_props.Vs ** 2) * qb_m

        # Compute mask for the energy integral time window
        integral_mask = (self._time_axis >= flux_window[0]) & (self._time_axis <= flux_window[1])
        vm_windowed = vm[:, :, integral_mask]

        # Take the su component.
        su_windowed = vm_windowed[:, 3, :]

        # Integrate in time
        Esu_per_event = Nsu * np.sum(np.abs(su_windowed) ** 2, axis=1)

        # Compute mean over events
        Esu = np.mean(Esu_per_event)

        return Esu, Esu_per_event, vm
    # end func

    def _mode_matrices(Vp, Vs, rho, p):
        """Compute M, M_inv and Q for a single layer for a scalar or array of ray parameters p.

        :param Vp: P-wave body wave velocity (scalar, labeled α in Tao's paper)
        :type Vp:
        :param Vs: S-wave body wave velocity (scalar, labeled β in Tao's paper)
        :type Vs:
        :param rho: Bulk material density, ρ (scalar)
        :type rho:
        :param p: Scalar or array of ray parameters (one per event)
        :type p:
        """
        qa = np.sqrt((1 / Vp ** 2 - p * p).astype(np.complex))
        assert not np.any(np.isnan(qa)), qa
        qb = np.sqrt((1 / Vs ** 2 - p * p).astype(np.complex))
        assert not np.any(np.isnan(qb)), qb
        eta = 1 / Vs ** 2 - 2 * p * p
        mu = rho * Vs * Vs
        trp = 2 * mu * p * qa
        trs = 2 * mu * p * qb
        mu_eta = mu * eta
        # First compute without velocity factors for reduced operation count.
        M = np.array([
            [p, p, qb, qb],
            [qa, -qa, -p, p],
            [-trp, trp, -mu_eta, mu_eta],
            [-mu_eta, -mu_eta, trs, trs]
        ])
        # Then times by velocity factors
        Vfactors = np.diag([Vp, Vp, Vs, Vs])
        M = np.matmul(np.moveaxis(M, -1, 0), Vfactors)

        Q = np.dstack([np.expand_dims(np.array([-_1, _1, -_2, _2]), 1) for (_1, _2) in zip(qa, qb)])
        Q = np.moveaxis(Q, -1, 0)

        # First compute without velocity factors for reduced operation count.
        mu_p = mu * p
        Minv = (1.0 / rho) * np.array([
            [mu_p, mu_eta / 2 / qa, -p / 2 / qa, -0.5 * np.ones(p.shape)],
            [mu_p, -mu_eta / 2 / qa, p / 2 / qa, -0.5 * np.ones(p.shape)],
            [mu_eta / 2 / qb, -mu_p, -0.5 * np.ones(p.shape), p / 2 / qb],
            [mu_eta / 2 / qb, mu_p, 0.5 * np.ones(p.shape), p / 2 / qb]
        ])
        # Then times by velocity factors
        Vfactors_inv = np.diag([1 / Vp, 1 / Vp, 1 / Vs, 1 / Vs])
        Minv = np.matmul(Vfactors_inv, np.moveaxis(Minv, -1, 0))

        #     # DEBUG CHECK - verify M*Minv is close to identity
        #     for i in range(M.shape[0]):
        #         _M = M[i,:,:]
        #         _Minv = Minv[i,:,:]
        #         assert _M.shape[0] == _M.shape[1]
        #         assert np.allclose(np.matmul(_M, _Minv).flatten(), np.eye(_M.shape[0]).flatten()), i

        return (M, Minv, Q)
    # end func

    def _propagate_layers(fv0, w, layer_props, p):
        """

        :param w:
        :param layer_props:
        :param p:
        :return:
        """
        # layer_props is a list of LayerProps
        fz = np.hstack((fv0, np.zeros_like(fv0)))
        for layer in layer_props:
            M, Minv, Q = WfContinuationSuFluxComputer._mode_matrices(layer.Vp, layer.Vs, layer.rho, p)
            fz = np.matmul(Minv, fz)
            # Expanding dims on w here means that at each level of the stack, phase_args is np.outer(Q, w)
            phase_args = np.matmul(Q, np.expand_dims(np.expand_dims(w, 0), 0))
            assert np.allclose(np.outer(Q[0,:,:], w).flatten(), phase_args[0,:,:].flatten()), (Q, w)
            # phase_args = np.matmul(Q - Q[1], np.expand_dims(np.expand_dims(w, 0), 0))  # formulation seen in Tao code
            # assert np.allclose(np.outer(Q[0,:,:] - Q[1], w).flatten(), phase_args[0,:,:].flatten()), (Q, w)
            phase_factors = np.exp(1j*layer.H*phase_args)
            fz = phase_factors*fz  # point-wise multiplication
            fz = np.matmul(M, fz)
        # end for
        return fz
    # end func


# end class
