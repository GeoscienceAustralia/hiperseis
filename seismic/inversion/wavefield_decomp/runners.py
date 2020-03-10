#!/usr/bin/env python
# coding: utf-8
"""
Batch execution interfaces for wavefield continuation methods and solvers.
"""

import json
import logging

import click
import numpy as np

from seismic.inversion.wavefield_decomp.model_properties import LayerProps


logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                    level=logging.INFO,
                    datefmt='%Y-%m-%d %H:%M:%S')


def run_mcmc(config):
    mantle_config = config['mantle_props']
    mantle_props = LayerProps(mantle_config['Vp'], mantle_config['Vs'], mantle_config['rho'], np.Infinity)
# end func


def mcmc_solver_wrapper(model, obj_fn, mantle, Vp, rho, flux_window):
    num_layers = len(model)//2
    earth_model = []
    for i in range(num_layers):
        earth_model.append(LayerProps(Vp[i], model[2*i + 1], rho[i], model[2*i]))
    # end for
    earth_model = np.array(earth_model)
    energy, _, _ = obj_fn(mantle, earth_model, flux_window=flux_window)
    return energy
# end func


@click.command()
@click.argument('config_file', type=click.File('r'), required=True)
@click.option('--waveform-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.option('--output-file', type=click.Path(dir_okay=False), required=True)
def main(config_file, waveform_file, output_file):
    """TBD

    The output file is in HDF5 format. The configuration details are added to the output file for traceability.

    :param config_file:
    :param waveform_file:
    :param output_file:
    :return: None
    """
    config = json.load(config_file)
    logging.info("Config:\n{}".format(json.dumps(config, indent=4)))
    logging.info("Waveform source: {}".format(waveform_file))
    logging.info("Destination file: {}".format(output_file))
# end func


if __name__ == '__main__':
    main()
# end if
