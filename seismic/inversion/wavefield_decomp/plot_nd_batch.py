#!/usr/bin/env python
# coding: utf-8
"""
Batch plotting a MCMC solution for batch of stations to a single pdf file.
"""

import os
import logging

import click
from tqdm.auto import tqdm

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

from seismic.inversion.wavefield_decomp.runners import load_mcmc_solution
from seismic.inversion.wavefield_decomp.wfd_plot import plot_Nd


@click.command()
@click.argument('solution_file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.option('--output-file', type=click.Path(dir_okay=False), required=True,
              help='Name of the output PDF file in which to save plots')
def main(solution_file, output_file):
    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)

    soln_config, job_id = load_mcmc_solution(solution_file, logger=log)

    out_basename, ext = os.path.splitext(output_file)
    out_basename += '_' + job_id
    output_file = out_basename + ext

    vars2 = ('$H_{crust}$', '$V_{s,crust}$')
    vars4 = ('$H_{sed}$', '$V_{s,sed}$', '$H_{crust}$', '$V_{s,crust}$')
    with PdfPages(output_file) as pdf:
        for soln, config in tqdm(soln_config):
            if soln.x.shape[-1] == 2:
                vars = vars2
            elif soln.x.shape[-1] == 4:
                vars = vars4
            else:
                assert False
            # end if
            p, _, _ = plot_Nd(soln, title=config['station_id'], vars=vars)
            pdf.savefig(dpi=300, papertype='a3', orientation='portrait')
            plt.close()
        # end for
    # end with

# end func


if __name__ == "__main__":
    main()
# end if