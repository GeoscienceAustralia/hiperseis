#!/usr/bin/env python
"""
Compute moho estimate and generate maps from converted data and 
a JSON config file.
"""
import os
import json
import logging

import click
from seismic.receiver_fn import (pointsets2grid, plot_spatial_map, write_gis_data, write_gmt_data)
from seismic.receiver_fn.moho_config import ConfigConstants as cc, validate, CORR_FUNC_MAP


def run_workflow(config_file):
    with open(config_file, 'r') as f:
        config = json.load(f)

    validate(config)

    # Data prep/correction
    data_prep = config.get(cc.DATA_PREP)
    if data_prep is not None:
        for params in data_prep:
            CORR_FUNC_MAP[params[cc.CORR_FUNC]](
                    params[cc.DATA_TO_PREP], params[cc.CORR_DATA], params[cc.CORR_OUT])

    # Moho interpolation 
    pointsets2grid.make_grid(config_file)

    # Plotting
    plotting = config.get(cc.PLOTTING, {})
    if plotting.get(cc.PLOT_FLAG, False):
        plot_spatial_map.from_config(config_file)
    if plotting.get(cc.GMT_FLAG, False):
        write_gmt_data.from_config(config_file)
    if plotting.get(cc.GIS_FLAG, False):
        write_gis_data.write_depth_grid(config_file)
        write_gis_data.write_gradient_grid(config_file)
        write_gis_data.write_sample_locations(config_file)

    
@click.command()
@click.argument('config-file', type=click.Path(exists=True, dir_okay=False), required=True)
def main(config_file):
    run_workflow(config_file)


if __name__ == '__main__':
    main()
