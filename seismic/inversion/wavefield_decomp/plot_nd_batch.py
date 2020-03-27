#!/usr/bin/env python
# coding: utf-8
"""
Batch plotting a MCMC solution for batch of stations to a single pdf file.
"""

import os
import logging

import click
from tqdm.auto import tqdm

import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import rf

from seismic.inversion.wavefield_decomp.runners import load_mcmc_solution
from seismic.inversion.wavefield_decomp.wfd_plot import plot_Nd
from seismic.stream_io import read_h5_stream


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
    st.rf()  # TODO: Filter RF and set to iterative deconv
    return st
# end func


def plot_aux_data(soln, config, log):
    f = plt.figure(constrained_layout=False, figsize=(12.8, 12.8))
    f.suptitle(config["station_id"], y=0.96, fontsize=16)
    gs = f.add_gridspec(2, 1, left=0.1, right=0.9, bottom=0.1, top=0.9, hspace=0.3,
                        height_ratios=[1, 2])
    gs_top = gs[0].subgridspec(1, 2)
    ax0 = f.add_subplot(gs_top[0, 0])
    ax1 = f.add_subplot(gs_top[0, 1])
    # ax2 = f.add_subplot(gs[1:, :])

    hist_alpha = 0.5
    soln_alpha = 0.3
    axis_font_size = 12
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
    ax0.set_title('Energy distribution of random samples and solution clusters')
    ax0.set_xlabel('$E_{SU}$ energy (arb. units)')
    ax0.set_ylabel('Normalized counts')
    ax0.tick_params(labelsize=axis_font_size)
    ax0.xaxis.label.set_size(axis_font_size)
    ax0.yaxis.label.set_size(axis_font_size)

    # Plot sorted per-event upwards S-wave energy at top of mantle per solution
    for i, esu in enumerate(soln.esu):
        color = 'C' + str(i)
        esu_sorted = sorted(esu)
        ax1.plot(esu_sorted, color=color, alpha=soln_alpha)
    # end for
    ax1.set_title('Ranked per-event energy for each solution point')
    ax1.set_xlabel('Rank (out of # source events)')
    ax1.set_ylabel('Event $E_{SU}$ energy (arb. units)')
    ax1.tick_params(labelsize=axis_font_size)
    ax1.xaxis.label.set_size(axis_font_size)
    ax1.yaxis.label.set_size(axis_font_size)

    # Plot receiver function at base of selected layers
    for layer in config["layers"]:
        lname = layer["name"]
        if soln.subsurface and lname in soln.subsurface:
            base_seismogms = soln.subsurface[lname]
            # Generate RF and plot.
            n_solutions = len(base_seismogms)
            gs_bot = gs[1].subgridspec(n_solutions, 1)
            for i, seismogm in enumerate(base_seismogms):
                soln_rf = _compute_rf(seismogm, config, log)
                assert isinstance(soln_rf, rf.RFStream)
                axn = f.add_subplot(gs_bot[i])
                if soln_rf:
                    color = 'C' + str(i)
                    rf_z = soln_rf.select(component='Z')
                    times = rf_z[0].times() + config["su_energy_opts"]["time_window"][0]
                    data = rf_z.stack()[0].data
                    plt.plot(times, data, color=color, alpha=soln_alpha, linewidth=2)
                else:
                    axn.annotate('Empty RF plot', (0.5, 0.5), xycoords='axes fraction', ha='center')
                # end if
                axn.set_title(' '.join([config["station_id"], lname, 'base RF']))
                axn.tick_params(labelsize=axis_font_size)
                axn.xaxis.label.set_size(axis_font_size)
                axn.yaxis.label.set_size(axis_font_size)
            # end for
            break  # TODO: Figure out how to add more layers if needed
        # end if
    # end for

    # f.tight_layout()

    return f
# end func


@click.command()
@click.argument('solution_file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.option('--output-file', type=click.Path(dir_okay=False), required=True,
              help='Name of the output PDF file in which to save plots')
def main(solution_file, output_file):
    """
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

    vars2 = ('$H_{crust}$', '$k_{s,crust}$')
    vars4 = ('$H_{sed}$', '$k_{s,sed}$', '$H_{crust}$', '$k_{s,crust}$')
    with PdfPages(output_file) as pdf:
        for soln, config in tqdm(soln_config):
            if soln.x.shape[-1] == 2:
                vars = vars2
            elif soln.x.shape[-1] == 4:
                vars = vars4
            else:
                assert False
            # end if
            #  Convert Vs parameter to k = Vp/Vs
            convert_Vs_to_k(soln, config)
            p, _, _ = plot_Nd(soln, title=config['station_id'], vars=vars)
            # Annotate top left axes with number of events use in the solver.
            ndims = len(soln.bounds.lb)
            ax = p.axes[0, ndims - 1]
            ax.text(0.95, 0.95, 'N = {}'.format(soln.num_input_seismograms), fontsize=10,
                    ha='right', va='top', transform=ax.transAxes)
            pdf.savefig(dpi=300, papertype='a3', orientation='portrait')
            plt.close()

            p = plot_aux_data(soln, config, log)
            pdf.savefig(dpi=300, papertype='a3', orientation='portrait')
            plt.close()
        # end for
    # end with
    log.info('Produced file {}'.format(output_file))

# end func


if __name__ == "__main__":
    main()
# end if
