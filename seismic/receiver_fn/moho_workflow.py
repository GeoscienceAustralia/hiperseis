#!/usr/bin/env python
"""
Compute moho estimate and generate maps from converted data and 
a JSON config file.
"""
import os
import json
import logging
import tempfile

import click
from seismic.receiver_fn import (pointsets2grid, plot_spatial_map, write_gis_data, write_gmt_data)
from seismic.receiver_fn.moho_config import (
    ConfigConstants as cc, validate, CORR_FUNC_MAP, MethodDataset)


def run_workflow(config_file):
    with open(config_file, 'r') as f:
        config = json.load(f)

    validate(config)

    # Filter out disabled methods
    methods = config[cc.METHODS]
    disabled_methods = [param[cc.NAME] for param in methods if param.get(cc.DISABLE, False)]
    if disabled_methods:
        print(f"Disabled methods: {disabled_methods}")
        methods = [param for param in methods if not param.get(cc.DISABLE, False)]
        # Raise an error if all methods are disabled
        if not methods:
            raise ValueError("All methods have been disabled - there is no data to process!")
        config[cc.METHODS] = methods
        # Save the modified config in a temp file so changes are applied through the workflow
        _, config_file = tempfile.mkstemp()
        with open(config_file, 'w') as f:
            json.dump(config, f)

    # Data prep/correction
    data_prep = config.get(cc.DATA_PREP)
    if data_prep is not None:
        # Correct d1 by d2
        for params in data_prep:
            d1 = MethodDataset({cc.DATA: params[cc.DATA_TO_PREP], cc.VAL_NAME: params[cc.DATA_VAL]})
            d2 = MethodDataset({cc.DATA: params[cc.CORR_DATA], cc.VAL_NAME: params[cc.CORR_VAL]})
            CORR_FUNC_MAP[params[cc.CORR_FUNC]](d1, d2, params[cc.CORR_OUT])

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

    # If methods were disabled, remove the temporary config file
    if disabled_methods:
        os.remove(config_file)

    
@click.command()
@click.argument('config-file', type=click.Path(exists=True, dir_okay=False), required=True)
def main(config_file):
    run_workflow(config_file)


if __name__ == '__main__':
    main()
