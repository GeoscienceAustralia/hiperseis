#!/usr/bin/env python
# coding: utf-8
"""
Batch plotting a MCMC solution for batch of stations to a single pdf file.
"""

import os
import logging
import json
import copy

import click
from tqdm.auto import tqdm

import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.table import table
import rf

from seismic.inversion.wavefield_decomp.runners import load_mcmc_solution
from seismic.inversion.wavefield_decomp.wfd_plot import plot_Nd
from seismic.stream_io import read_h5_stream
from seismic.receiver_fn.rf_deconvolution import rf_iter_deconv


default_fig_width_inches = 6.4
RF_TRIM_WINDOW = (-10.0, 25.0)

def convert_Vs_to_k(soln, config):
    """
    Transform Vs variable into k variable in MCMC solution. Modifies soln in-place.

    :param soln: Solution container
    :type soln: Customized scipy.optimize.OptimizeResult
    :param config: Solution configuration
    :type config: dict
    :return: None
    """
    layers = config['layers']
    for i, layer in enumerate(layers):
        Vp = layer['Vp']
        if 'k_range' not in layer:
            assert 'Vs_range' in layer
            Vs_range = layer['Vs_range']
            layer.update({'k_range': [Vp/Vs_range[1], Vp/Vs_range[0]]})
        # end if
        # Scale Vs content in bounds, x, clusters, samples, bins and distribution to Vp/Vs
        lb = soln.bounds.lb[2*i + 1]
        soln.bounds.lb[2*i + 1] = Vp/soln.bounds.ub[2*i + 1]
        soln.bounds.ub[2 * i + 1] = Vp/lb
        soln.x[:, 2*i + 1] = Vp/soln.x[:, 2*i + 1]
        for cluster in soln.clusters:
            cluster[:, 2*i + 1] = Vp/cluster[:, 2*i + 1]
        # end for
        soln.samples[:, 2*i + 1] = Vp/soln.samples[:, 2*i + 1]
        soln.bins[2*i + 1] = np.flip(Vp/soln.bins[2*i + 1])
        soln.distribution[2*i + 1] = np.flip(soln.distribution[2*i + 1])
    # end for
# end func


def _compute_rf(data, config, log):
    st = rf.RFStream()
    event_ids = config.get("event_ids")
    src_file = config.get("waveform_file")
    if event_ids is None:
        log.error("Unable to generate RF without event IDs")
        return st
    # end if
    if src_file is None:
        log.error("Unable to generate RF without path to source file")
        return st
    # end if
    if not os.path.isfile(src_file):
        log.error("Source file {} for trace metadata not found, cannot generate RF".format(src_file))
        return
    # end if
    net, sta, loc = config["station_id"].split('.')
    src_waveforms = read_h5_stream(src_file, net, sta, loc)
    assert data.shape[0] == len(event_ids)
    for i, event_data in enumerate(data):
        evid = event_ids[i]
        src_stream = rf.RFStream([tr for tr in src_waveforms if tr.stats.event_id == evid])
        # Z component
        z_header = src_stream.select(component='Z')[0].stats
        su_opts = config["su_energy_opts"]
        z_header.starttime = z_header.onset + su_opts["time_window"][0]
        z_header.sampling_rate = su_opts["sampling_rate"]
        z_header.delta = 1.0/z_header.sampling_rate
        z_header.npts = event_data.shape[1]
        assert np.isclose(float(z_header.endtime - z_header.starttime), su_opts["time_window"][1] - su_opts["time_window"][0])
        tr = rf.rfstream.RFTrace(event_data[1, :], header=z_header)
        st += tr
        # R component
        r_header = z_header.copy()
        r_header.channel = z_header.channel[:-1] + 'R'
        tr = rf.rfstream.RFTrace(event_data[0, :], header=r_header)
        st += tr
    # end for

    st.filter('bandpass', freqmin=0.05, freqmax=1.0, corners=2, zerophase=True)
    normalize = 0  # Use Z-component for normalization
    st.rf(rotate=None, method='P', deconvolve='func', func=rf_iter_deconv,
          normalize=normalize, min_fit_threshold=75.0)

    return st
# end func


def plot_aux_data(soln, config, log, scale):
    f = plt.figure(constrained_layout=False, figsize=(6.4*scale, 6.4*scale))
    f.suptitle(config["station_id"], y=0.96, fontsize=16)
    gs = f.add_gridspec(2, 1, left=0.1, right=0.9, bottom=0.1, top=0.87, hspace=0.3,
                        wspace=0.3, height_ratios=[1, 2])
    gs_top = gs[0].subgridspec(1, 2)
    ax0 = f.add_subplot(gs_top[0, 0])
    ax1 = f.add_subplot(gs_top[0, 1])

    hist_alpha = 0.5
    soln_alpha = 0.3
    axis_font_size = 6*scale
    title_font_size = 6*scale
    nbins = 100

    # Plot energy distribution of samples and solution clusters
    energy_hist, bins = np.histogram(soln.sample_funvals, bins=nbins)
    energy_hist = energy_hist.astype(float)/np.max(energy_hist)
    ax0.bar(bins[:-1], energy_hist, width=np.diff(bins), align='edge', color='#808080', alpha=hist_alpha)

    for i, cluster_energies in enumerate(soln.cluster_funvals):
        color = 'C' + str(i)
        cluster_hist, _ = np.histogram(cluster_energies, bins)
        cluster_hist = cluster_hist.astype(float)/np.max(cluster_hist)
        ax0.bar(bins[:-1], cluster_hist, width=np.diff(bins), align='edge', color=color, alpha=soln_alpha)
    # end for
    ax0.set_title('Energy distribution of random samples and solution clusters',
                  fontsize=title_font_size)
    ax0.set_xlabel('$E_{SU}$ energy (arb. units)')
    ax0.set_ylabel('Normalized counts')
    ax0.tick_params(labelsize=axis_font_size)
    ax0.xaxis.label.set_size(axis_font_size)
    ax0.yaxis.label.set_size(axis_font_size)

    # Plot sorted per-event upwards S-wave energy at top of mantle per solution.
    # Collect event IDs of worst fit traces and present as table of waveform IDs.
    event_ids = config["event_ids"]
    events_best3 = []
    events_worst3 = []
    for i, esu in enumerate(soln.esu):
        assert len(esu) == len(event_ids)
        color = 'C' + str(i)
        esu_sorted = sorted(zip(esu, event_ids))
        events_best3.extend(esu_sorted[:3])
        events_worst3.extend(esu_sorted[-3:])
        esu_sorted = [e[0] for e in esu_sorted]
        ax1.plot(esu_sorted, color=color, alpha=soln_alpha)
    # end for
    events_best3 = sorted(events_best3)
    events_worst3 = sorted(events_worst3, reverse=True)
    best_events_set = set()
    worst_events_set = set()
    for _, evid in events_best3:
        best_events_set.add(evid)
        if len(best_events_set) >= 3:
            break
        # end if
    # end for
    for _, evid in events_worst3:
        worst_events_set.add(evid)
        if len(worst_events_set) >= 3:
            break
        # end if
    # end for
    _tab1 = table(ax1, cellText=[[e] for e in best_events_set], colLabels=['BEST'],
                  cellLoc='left', colWidths=[0.35], loc='upper left',
                  edges='horizontal', fontsize=8, alpha=0.6)  # alpha broken in matplotlib.table!
    _tab2 = table(ax1, cellText=[[e] for e in worst_events_set], colLabels=['WORST'],
                  cellLoc='left', colWidths=[0.35], loc='upper right',
                  edges='horizontal', fontsize=8, alpha=0.6)
    ax1.set_title('Ranked per-event energy for each solution point',
                  fontsize=title_font_size)
    ax1.set_xlabel('Rank (out of # source events)')
    ax1.set_ylabel('Event $E_{SU}$ energy (arb. units)')
    ax1.tick_params(labelsize=axis_font_size)
    ax1.xaxis.label.set_size(axis_font_size)
    ax1.yaxis.label.set_size(axis_font_size)

    # Plot receiver function at base of selected layers
    axis_font_size = 6*scale
    max_solutions = config["solver"].get("max_solutions", 3)
    for layer in config["layers"]:
        lname = layer["name"]
        if soln.subsurface and lname in soln.subsurface:
            base_seismogms = soln.subsurface[lname]
            # Generate RF and plot.
            gs_bot = gs[1].subgridspec(max_solutions, 1, hspace=0.4)
            for i, seismogm in enumerate(base_seismogms):
                soln_rf = _compute_rf(seismogm, config, log)
                assert isinstance(soln_rf, rf.RFStream)
                # Remove any traces for which deconvolution failed.
                # First, find their unique ID. Then remove all traces with that ID.
                exclude_ids = set([tr.stats.event_id for tr in soln_rf if len(tr) == 0])
                soln_rf = rf.RFStream([tr for tr in soln_rf if tr.stats.event_id not in exclude_ids])
                axn = f.add_subplot(gs_bot[i])
                if soln_rf:
                    color = 'C' + str(i)
                    rf_R = soln_rf.select(component='R').trim2(RF_TRIM_WINDOW[0], RF_TRIM_WINDOW[1], reftime='onset')
                    num_RFs = len(rf_R)
                    times = rf_R[0].times() + RF_TRIM_WINDOW[0]
                    data = rf_R.stack()[0].data
                    axn.plot(times, data, color=color, alpha=soln_alpha, linewidth=2)
                    axn.text(0.95, 0.95, 'N = {}'.format(num_RFs), fontsize=10,
                             ha='right', va='top', transform=axn.transAxes)
                    axn.set_xlabel('Time (sec)')
                    axn.grid(color='#80808080', linestyle=':')
                else:
                    axn.annotate('Empty RF plot', (0.5, 0.5), xycoords='axes fraction', ha='center')
                # end if
                axn.set_title(' '.join([config["station_id"], lname, 'base RF', '(soln {})'.format(i)]),
                              fontsize=title_font_size, y=0.92, va='top')
                axn.tick_params(labelsize=axis_font_size)
                axn.xaxis.label.set_size(axis_font_size)
                axn.yaxis.label.set_size(axis_font_size)
            # end for
            break  # TODO: Figure out how to add more layers if needed
        # end if
    # end for

    return f
# end func


@click.command()
@click.argument('solution_file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.option('--output-file', type=click.Path(dir_okay=False), required=True,
              help='Name of the output PDF file in which to save plots')
def main(solution_file, output_file):
    """Plot all the solutions found in a batch run of N-dimensional solver.

    Example usage:
        python seismic/inversion/wavefield_decomp/plot_nd_batch.py --output-file OA_wfd_out.pdf OA_wfd_out.h5

    :param solution_file: Input solution filename
    :param output_file: Output filename
    """
    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)

    soln_config, job_id = load_mcmc_solution(solution_file, logger=log)

    out_basename, ext = os.path.splitext(output_file)
    out_basename += '_' + job_id
    output_file = out_basename + ext

    with PdfPages(output_file) as pdf:
        for soln, config in tqdm(soln_config):
            empty_soln = (soln.x.shape[0] == 0)
            if empty_soln:
                continue
            assert soln.x.shape[-1] == len(config['layers'])*2
            vars = []
            for layer in config['layers']:
                layer_name = layer['name']
                vars += ['$H_{{{}}}$'.format(layer_name), '$k_{{{}}}$'.format(layer_name)]
            # end for
            vars = tuple(vars)

            # Dump settings page (per station)
            config_no_evids = copy.deepcopy(config)
            config_no_evids.pop('event_ids', None)
            config_text = json.dumps(config_no_evids, indent=4)
            f = plt.figure(figsize=(default_fig_width_inches, default_fig_width_inches*1.414))
            plt.gca().xaxis.set_visible(False)
            plt.gca().yaxis.set_visible(False)
            plt.title(config['station_id'] + ' Processing Settings')
            f.text(0.02, 0.98, 'Settings:\n' + config_text, fontsize=6, va='top',
                   fontname='monospace', transform=plt.gca().transAxes)
            pdf.savefig(dpi=300, papertype='a3', orientation='portrait')
            plt.close()

            #  Convert Vs parameter to k = Vp/Vs
            convert_Vs_to_k(soln, config)
            p, _, _ = plot_Nd(soln, title=config['station_id'], vars=vars)
            scale = p.fig.get_size_inches()[0]/default_fig_width_inches
            # Annotate top left axes with number of events use in the solver.
            ndims = len(soln.bounds.lb)
            ax = p.axes[0, ndims - 1]
            ax.text(0.95, 0.95, 'N = {}'.format(soln.num_input_seismograms), fontsize=10,
                    ha='right', va='top', transform=ax.transAxes)
            pdf.savefig(dpi=300, papertype='a3', orientation='portrait')
            plt.close()

            _p = plot_aux_data(soln, config, log, scale)
            pdf.savefig(dpi=300, papertype='a3', orientation='portrait')
            plt.close()
        # end for
    # end with
    log.info('Produced file {}'.format(output_file))

# end func


if __name__ == "__main__":
    main()
# end if
