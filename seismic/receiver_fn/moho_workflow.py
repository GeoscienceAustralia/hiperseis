#!/usr/bin/env python
"""
Compute moho estimate and generate maps from converted data and 
a JSON config file.
"""
import os
import json

import click
from seismic.receiver_fn import (pointsets2grid, plot_spatial_map, write_gis_data, write_gmt_data)


@click.command()
@click.argument('config-file', type=click.Path(exists=True, dir_okay=False), required=True)
def main(config_file):
    """
    """
    with open(config_file, 'r') as f:
        config = json.load(f)

    pointsets2grid.make_grid(config_file)
    plotting = config.get('plotting')
    if plotting is not None:
        if plotting.get('output_cartopy_plot', False):
            plot_spatial_map.from_config(config_file)
        if plotting.get('output_gmt_data', False):
            write_gmt_data.from_config(config_file)
        if plotting.get('output_gis_data', False):
            write_gis_data.write_depth_grid(config_file)
            write_gis_data.write_gradient_grid(config_file)




if __name__ == '__main__':
    main()
