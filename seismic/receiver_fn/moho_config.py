"""
This module contains functions and constants to assist with parsing
the Moho workflow config.
"""
import os


class ConfigConstants:
    # Field names
    METHODS = 'methods'
    PLOTTING = 'plotting'
    BOUNDS = 'bounds'
    GRID_INTERVAL = 'grid_interval'
    OUTPUT_DIR = 'output_directory'
    NAME = 'name'
    DATA = 'data'
    WEIGHT = 'weight'
    SCALE_LENGTH = 'scale_length'
    PLOT_FLAG = 'output_plot'
    PLOT_PARAMS = 'plot_parameters'
    PP_SCALE = 'scale'
    PP_FMT = 'format'
    PP_SHOW = 'show'
    GMT_FLAG = 'output_gmt'
    GIS_FLAG = 'output_gis'

    # Output filenames
    MOHO_GRID = 'moho_grid.csv'
    MOHO_GRAD = 'moho_gradient.csv'
    MOHO_PLOT = 'moho_plot'
    GMT_DIR = 'gmt_data'
    MOHO_GRID_GMT = 'moho_grid.txt'
    MOHO_GRAD_GMT = 'moho_gradient.txt'
    LOCATIONS_GMT = '_locations.txt'
    GIS_DIR = 'gis_data'
    MOHO_GRID_GIS = 'moho_grid.tif'
    MOHO_GRAD_GIS = 'moho_gradient.tif'
    LOCATIONS_GIS = '_locations'

_cc = ConfigConstants

# Config validation #

TOP_LEVEL_SUPPORTED_KEYS = [_cc.METHODS, _cc.PLOTTING, _cc.BOUNDS, 
                            _cc.GRID_INTERVAL, _cc.OUTPUT_DIR, 
                            _cc.NAME, _cc.DATA, _cc.WEIGHT,
                            _cc.SCALE_LENGTH]

METHOD_SUPPORTED_KEYS = [_cc.NAME, _cc.DATA, _cc.WEIGHT,
                         _cc.SCALE_LENGTH]

PLOTTING_SUPPORTED_KEYS = [_cc.PLOT_FLAG, _cc.PLOT_PARAMS, 
                           _cc.GMT_FLAG, _cc.GIS_FLAG]

CP_PARAM_SUPPORTED_KEYS = [_cc.PP_SCALE, _cc.PP_FMT, _cc.PP_SHOW]

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

    methods = _try_lookup(config, _cc.METHODS, 
            f"{_cc.METHODS} dictionary does not exist")
    if not methods:
        raise Exception(f"At least one entry in {_cc.METHODS}  is required")

    for params in methods:
        _check_type(params, [dict], f"'{_cc.METHODS}' parameters must be stored in dict")

        _check_keys(params.keys(), METHOD_SUPPORTED_KEYS)

        name = _try_lookup(params, _cc.NAME, f"'{_cc.METHODS} entry requires a name as '{_cc.NAME}' field")
        _check_type(name, [str], f"method name must be of type str")

        csv_file = _try_lookup(params, _cc.DATA, 
                f"method '{name}' requires data file as '{_cc.DATA}' field")
        _check_type(csv_file, [str], "'{_cc.DATA} must be path to CSV data as str")
        if not os.path.exists(csv_file):
            raise FileNotFoundError(f"CSV data '{_cc.DATA}' for method '{name}' does not exist")

        weighting = _try_lookup(params, _cc.WEIGHT, 
                f"method '{name}' requires a weighting as '{_cc.WEIGHT}' field. For neutral weighting, "
                "use a value of 1.0")
        _check_type(weighting, [float, int], 
                f"weighting for method '{name}' must be of type float or int")

        scale_len = _try_lookup(params, _cc.SCALE_LENGTH,
                f"method '{name}' requires spatial spread scale length in decimal degrees as field "
                f"'{_cc.SCALE_LENGTH}' field")
        _check_type(scale_len, [float, int],
                f"scale length for method '{name}' must be of type float or int")

    plt = config.get(_cc.PLOTTING)
    if plt is not None:
        _check_type(plt, [dict], "plotting parameters must be stored in dict")
        _check_keys(plt.keys(), PLOTTING_SUPPORTED_KEYS)
        cartopy_flag = plt.get(_cc.PLOT_FLAG)
        if cartopy_flag is not None:
            _check_type(cartopy_flag, [bool], f"'{_cc.PLOT_FLAG}' must be of type bool")

            cp_params = plt.get(_cc.PLOT_PARAMS)
            if cp_params is not None:
                _check_type(cp_params, [dict], f"'{_cc.PLOT_PARAMS}' must be stored in dict")
                _check_keys(cp_params.keys(), CP_PARAM_SUPPORTED_KEYS)
                scale = cp_params.get(_cc.PP_SCALE)
                if scale is not None:
                    _check_type(scale, [list, tuple], f"'{_cc.PP_SCALE}' must be of type list or tuple")
                    if len(scale) != 2:
                        raise ValueError(f"'{_cc.PP_SCALE}' must have exactly two elements")
                    for s in scale:
                        _check_type(s, [int, float], f"{_cc.PP_SCALE}' values must be of type int or float")
                fmt = cp_params.get(_cc.PP_FMT)
                if fmt is not None:
                    _check_type(fmt, [str], 
                            f"'{_cc.PP_FMT}' must be desired output filetype extension as str")
                show = cp_params.get(_cc.PP_SHOW)
                if show is not None:
                    _check_type(show, [bool], f"'{_cc.PP_SHOW}' must be of type bool")

            gmt_flag = plt.get(_cc.GMT_FLAG)
            if gmt_flag is not None:
                _check_type(gmt_flag, [bool], f"'{_cc.GMT_FLAG}' must be of type bool")

            gis_flag = plt.get(_cc.GIS_FLAG)
            if gis_flag is not None:
                _check_type(gis_flag, [bool], f"'{_cc.GIS_FLAG}' must be of type bool")

    bounds = config.get(_cc.BOUNDS)
    if bounds is not None:
        _check_type(bounds, [list, tuple], f"'{_cc.BOUNDS}' must be of type of list or tuple")
        if len(bounds) != 4:
            raise ValueError(f"'{_cc.BOUNDS}' must have exactly 4 elements (west, south, east, north)")
        for b in bounds:
            _check_type(b, [float, int], "bounds values must be of type int or float")
    else:
        print("No bounds provided, bounds will be derived from extent of data")

    spacing = _try_lookup(config, _cc.GRID_INTERVAL, 
        f"grid cell size in decimal degrees is required as field '{_cc.GRID_INTERVAL}'")
    _check_type(spacing, [float, int], f"'{_cc.GRID_INTERVAL}' must be of value float or int")

    output_dir = config.get(_cc.OUTPUT_DIR)
    if output_dir is not None:
        _check_type(output_dir, [str], f"'{_cc.OUTPUT_DIR}' must be a path as type str")
    else:
        print("No output directory provided, current working directory will be used")

