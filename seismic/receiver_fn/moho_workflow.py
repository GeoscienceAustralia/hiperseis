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

        
@click.command()
@click.argument('config-file', type=click.Path(exists=True, dir_okay=False), required=True)
def main(config_file):
    with open(config_file, 'r') as f:
        config = json.load(f)

    validate(config)

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
            write_gis_data.write_sample_locations(config_file)


# Config validation #

TOP_LEVEL_SUPPORTED_KEYS = ['methods', 'plotting', 'bounds', 'output_spacing_degrees', 'output_dir']
METHOD_SUPPORTED_KEYS = ['name', 'csv_file', 'weighting', 'scale_length_degrees', 
                            'enable_sample_weighting']
PLOTTING_SUPPORTED_KEYS = ['output_cartopy_plot', 'cartopy_parameters', 'output_gmt_data', 
                           'output_gis_data']
CP_PARAM_SUPPORTED_KEYS = ['scale', 'format', 'show']

def _try_lookup(d, f, msg):
    try:
        return d[f]
    except KeyError as e:
        raise Exception(msg) from e


def _check_type(x, types, msg):
    if type(x) not in types:
        raise TypeError(msg)


def _check_keys(keys, sup_keys): 
    for k in keys:
        if k not in sup_keys:
            raise ValueError(f"Config field '{k}' is not recognised, supported fields: {sup_keys}")


def validate(config):
    """
    Check that required parameters exist and are the correct type.
    Check that only supported keys are included (helps to catch typos
    causing silent errors).
    """
    _check_keys(config.keys(), TOP_LEVEL_SUPPORTED_KEYS)

    methods = _try_lookup(config, 'methods', "'methods' dictionary does not exist")
    if not methods:
        raise Exception("At least one entry in 'methods' is required")

    for params in methods:
        _check_type(params, [dict], "'method' parameters must be stored in dict")

        _check_keys(params.keys(), METHOD_SUPPORTED_KEYS)

        name = _try_lookup(params, 'name', "'method' entry requires a name as 'name' field")
        _check_type(name, [str], "method 'name' must be of type str")

        csv_file = _try_lookup(params, 'csv_file', 
                f"method '{name}' requires data file as 'csv_file' field")
        _check_type(csv_file, [str], "'csv_file' must be path to CSV data as str")
        if not os.path.exists(csv_file):
            raise FileNotFoundError(f"CSV data '{csv_file}' for method '{name}' does not exist")

        weighting = _try_lookup(params, 'weighting', 
                f"method '{name}' requires a weighting as 'weighting' field. For neutral weighting, "
                "use a value of 1.0")
        _check_type(weighting, [float, int], 
                f"weighting for method '{name}' must be of type float or int")

        scale_len = _try_lookup(params, 'scale_length_degrees',
                f"method '{name}' requires spatial spread scale length in decimal degrees as field "
                "'scale_length_degrees")
        _check_type(scale_len, [float, int],
                f"scale length for method '{name}' must be of type float or int")

    plt = config.get('plotting')
    if plt is not None:
        _check_type(plt, [dict], "'plotting' parameters must be stored in dict")
        _check_keys(plt.keys(), PLOTTING_SUPPORTED_KEYS)
        cartopy_flag = plt.get('output_cartopy_plot')
        if cartopy_flag is not None:
            _check_type(cartopy_flag, [bool], "'output_cartopy_plot' must be of type bool")

            cp_params = plt.get('cartopy_parameters')
            if cp_params is not None:
                _check_type(cp_params, [dict], "'cartopy_parameters' must be stored in dict")
                _check_keys(cp_params.keys(), CP_PARAM_SUPPORTED_KEYS)
                scale = cp_params.get('scale')
                if scale is not None:
                    _check_type(scale, [list, tuple], "'scale' must be of type list or tuple")
                    if len(scale) != 2:
                        raise ValueError("'scale' must have exactly two elements")
                    for s in scale:
                        _check_type(s, [int, float], "'scale' values must be of type int or float")
                fmt = cp_params.get('format')
                if fmt is not None:
                    _check_type(fmt, [str], 
                            "'format' must be desired output filetype extension as str")
                show = cp_params.get('show')
                if show is not None:
                    _check_type(show, [bool], "'show' must be of type bool")

            gmt_flag = plt.get('output_gmt_data')
            if gmt_flag is not None:
                _check_type(gmt_flag, [bool], "'output_gmt_data' must be of type bool")

            gis_flag = plt.get('output_gis_data')
            if gis_flag is not None:
                _check_type(gis_flag, [bool], "'output_gis_data' must be of type bool")

    bounds = config.get('bounds')
    if bounds is not None:
        _check_type(bounds, [list, tuple], "'bounds' must be of type of list or tuple")
        if len(bounds) != 4:
            raise ValueError("'bounds' must have exactly 4 elements (west, south, east, north)")
        for b in bounds:
            _check_type(b, [float, int], "bounds values must be of type int or float")
    else:
        print("No bounds provided, bounds will be derived from extent of data")

    spacing = _try_lookup(config, 'output_spacing_degrees', 
        "grid cell size in decimal degrees is required as field 'output_spacing_degrees'")
    _check_type(spacing, [float, int], "'output_spacing_degrees' must be of value float or int")

    output_dir = config.get('output_dir')
    if output_dir is not None:
        _check_type(output_dir, [str], "'output_dir' must be a path as type str")
    else:
        print("No output directory provided, current working directory will be used")


if __name__ == '__main__':
    main()
