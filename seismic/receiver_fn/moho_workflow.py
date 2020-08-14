#!/usr/bin/env python
"""
Compute moho estimate and generate maps from converted data and 
a JSON config file.
"""
import os
import json

import click
from seismic.receiver_fn import (pointsets2grid, plot_spatial_map)


@click.command()
@click.argument('config-file', type=click.Path(exists=True, dir_okay=False), required=True)
def main(config_file):
    """
    """
    pointsets2grid.make_grid(config_file)
    plot_spatial_map.from_config(config_file)


if __name__ == '__main__':
    main()
