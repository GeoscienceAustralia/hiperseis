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

from seismic.inversion.wavefield_decomp.runners import load_mcmc_solution
from seismic.inversion.wavefield_decomp.wfd_plot import plot_Nd


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
        # end for
    # end with

# end func


if __name__ == "__main__":
    main()
# end if
