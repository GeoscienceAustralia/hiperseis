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
from seismic.receiver_fn.moho_config import WorkflowParameters
    

def run_workflow(config_file):
    params = WorkflowParameters(config_file)
    
    # Moho interpolation 
    pointsets2grid.make_grid(params)

    # Plotting
    if params.plot_spatial_map:
        plot_spatial_map.from_params(params)
    if params.gmt_write:
        write_gmt_data.from_params(params)
    if params.gis_write:
        write_gis_data.from_params(params)


@click.command()
@click.argument('config-file', type=click.Path(exists=True, dir_okay=False), required=True)
def main(config_file):
    run_workflow(config_file)


if __name__ == '__main__':
    main()
