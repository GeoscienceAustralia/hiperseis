"""
This module contains functions and constants to assist with parsing
the Moho workflow config.
"""
import os
import json
import math

import numpy as np
from obspy import read_inventory

import seismic.receiver_fn.ccp_correction


try:
    import kennett_dist
    K_DIST = True
except ImportError:
    print("Kennett distance function not found. Using haversine.")
    K_DIST = False

def _kennett_dist(s, r, **kwargs):
    """
    Spherical distance. Requires 'kennett_dist' f2py module.
    Copyright B.L.N Kennett 1978, R.S.E.S A.N.U
    """
    clons, clats = s
    clonr, clatr = r
    max_dist = kwargs.get('max_dist', np.inf)
    # Distance in degrees
    delta, _, _ = kennett_dist.ydiz(clats, clons, clatr, clonr)
    if delta > max_dist:
        return np.nan
    else:
        return delta


def _haversine(s, r, **kwargs):
    """
    Haversine distance.
    """
    lon1, lat1 = s
    lon2, lat2 = r
    max_dist = kwargs.get('max_dist', np.inf)
    R = 6372.8

    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    lat1 = math.radians(lat1)
    lat2 = math.radians(lat2)

    a = math.sin(dlat / 2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2)**2
    c = 2 * math.asin(math.sqrt(a))

    delta_degrees = c * 180.0/math.pi

    if delta_degrees > max_dist:
        return np.nan
    else:
        return delta_degrees


DIST_METRIC = _kennett_dist if K_DIST else _haversine


class ConfigConstants:
    # Field names
    METHODS = 'methods'
    PLOTTING = 'plotting'
    BOUNDS = 'bounds'
    GRID_INTERVAL = 'grid_interval'
    INTERPOLATION= 'interpolation'
    WEIGHT_CUTOFF = 'weight_cutoff'
    OUTPUT_DIR = 'output_directory'
    NAME = 'name'
    DATA = 'data'
    PRIORITY = 'priority'
    WEIGHT = 'weight'
    VAL_NAME = 'val_name'
    SW_NAME = 'sw_name'
    DISABLE = 'disable'
    INV_FILE = 'inventory_file'
    SCALE_LENGTH = 'scale_length'
    PLOT_FLAG = 'output_plot'
    PLOT_PARAMS = 'plot_parameters'
    PP_SCALE = 'scale'
    PP_FMT = 'format'
    PP_SHOW = 'show'
    PP_TITLE = 'title'
    PP_CB_LABEL = 'cb_label'
    GMT_FLAG = 'output_gmt'
    GIS_FLAG = 'output_gis'
    DATA_PREP = 'data_preperation'
    DATA_TO_PREP = 'data'
    DATA_VAL = 'data_val'
    CORR_DATA = 'correction_data'
    CORR_VAL  = 'correction_val'
    CORR_FUNC = 'correction_func'
    CORR_OUT = 'output_file'

    # Correction functions
    CCP_CORR = 'ccp_correction'

    # Output filenames
    MOHO_GRID = 'grid.csv'
    MOHO_GRAD = 'gradient.csv'
    MOHO_PLOT = 'plot'
    GMT_DIR = 'gmt_data'
    MOHO_GRID_GMT = 'grid.txt'
    MOHO_GRAD_GMT = 'gradient.txt'
    LOCATIONS_GMT = '_locations.txt'
    GIS_DIR = 'gis_data'
    MOHO_GRID_GIS = 'grid.tif'
    MOHO_GRAD_GIS = 'gradient.tif'
    LOCATIONS_GIS = '_locations'
    

_cc = ConfigConstants


# Correction function mapping
CORR_FUNC_MAP = {_cc.CCP_CORR: seismic.receiver_fn.ccp_correction.correct}


# Config validation #

TOP_LEVEL_SUPPORTED_KEYS = [_cc.METHODS, _cc.PLOTTING, _cc.BOUNDS, 
                            _cc.GRID_INTERVAL, _cc.OUTPUT_DIR, _cc.INTERPOLATION, _cc.WEIGHT_CUTOFF,
                            _cc.NAME, _cc.DATA, _cc.WEIGHT,
                            _cc.SCALE_LENGTH, _cc.DATA_PREP]

DATA_PREP_SUPPORTED_KEYS = [_cc.DATA_TO_PREP, _cc.CORR_DATA, _cc.CORR_FUNC, _cc.CORR_OUT,
                            _cc.DATA_VAL, _cc.CORR_VAL]

METHOD_SUPPORTED_KEYS = [_cc.NAME, _cc.DATA, _cc.WEIGHT, _cc.DISABLE, _cc.PRIORITY,
                         _cc.SCALE_LENGTH, _cc.VAL_NAME, _cc.SW_NAME, _cc.INV_FILE]

PLOTTING_SUPPORTED_KEYS = [_cc.PLOT_FLAG, _cc.PLOT_PARAMS, 
                           _cc.GMT_FLAG, _cc.GIS_FLAG]

CP_PARAM_SUPPORTED_KEYS = [_cc.PP_SCALE, _cc.PP_FMT, _cc.PP_SHOW, _cc.PP_TITLE, _cc.PP_CB_LABEL]


def _try_lookup(d, f, msg):
    try:
        return d[f]
    except KeyError as e:
        raise Exception(msg) from e


def _check_type(x, types, msg):
    if type(x) not in types:
        raise TypeError(msg)


def _check_keys(keys, sup_keys, part): 
    for k in keys:
        if k not in sup_keys:
            raise ValueError(f"Config field {k} is not recognised as part of {part} parameters, "
                             f"supported fields: {sup_keys}")


def validate(config):
    """
    Check that required parameters exist and are the correct type.
    Check that only supported keys are included (helps to catch typos
    causing silent errors).
    """
    _check_keys(config.keys(), TOP_LEVEL_SUPPORTED_KEYS, 'top-level')

    data_prep = config.get(_cc.DATA_PREP, [])
    for params in data_prep:
        _check_keys(params.keys(), DATA_PREP_SUPPORTED_KEYS, _cc.DATA_PREP)
        df = _try_lookup(params, _cc.DATA_TO_PREP,
                    f"{_cc.DATA_PREP} requires a data file to process as {_cc.DATA_TO_PREP} field")
        _check_type(df, [str], f"{_cc.DATA_TO_PREP} should be path to data as str")
        cd = _try_lookup(params, _cc.CORR_DATA,
                f"{_cc.DATA_PREP} requires correction data to process as {_cc.CORR_DATA} field")
        _check_type(cd, [str], f"{_cc.CORR_DATA} should be path to data as str")
        cf = _try_lookup(params, _cc.CORR_FUNC,
                f"{_cc.DATA_PREP} requires name of correction function as {_cc.CORR_FUNC} field")
        _check_type(cf, [str], f"{_cc.CORR_FUNC} should be name of correction function as str")
        if cf not in CORR_FUNC_MAP.keys():
            raise KeyError(f"Correction function {_cc.CORR_FUNC} not recognised. Valid correction "
                           f"functions are: {CORR_FUNC_MAP.keys()}")

    methods = _try_lookup(config, _cc.METHODS, 
            f"{_cc.METHODS} list does not exist")
    if not methods:
        raise Exception(f"At least one entry in {_cc.METHODS} is required")

    for params in methods:
        _check_type(params, [dict], f"{_cc.METHODS} parameters must be stored in dict")

        _check_keys(params.keys(), METHOD_SUPPORTED_KEYS, _cc.METHODS)

        name = _try_lookup(params, _cc.NAME, f"{_cc.METHODS} entry requires a name as {_cc.NAME} field")
        _check_type(name, [str], f"method name must be of type str")

        csv_file = _try_lookup(params, _cc.DATA, 
                f"method {name} requires data file as {_cc.DATA} field")
        _check_type(csv_file, [str], "{_cc.DATA} must be path to data as str")

        weighting = _try_lookup(params, _cc.WEIGHT, 
                f"method {name} requires a weighting as {_cc.WEIGHT} field. For neutral weighting, "
                "use a value of 1.0")
        _check_type(weighting, [float, int], 
                f"weighting for method {name} must be of type float or int")

        scale_len = _try_lookup(params, _cc.SCALE_LENGTH,
                f"method {name} requires spatial spread scale length in decimal degrees as field "
                f"{_cc.SCALE_LENGTH} field")
        _check_type(scale_len, [float, int],
                f"scale length for method {name} must be of type float or int")

        disable_flag = params.get(_cc.DISABLE)
        if disable_flag is not None:
            _check_type(disable_flag, [bool], f"{_cc.DISABLE} flag must be of type bool")

    plt = config.get(_cc.PLOTTING)
    if plt is not None:
        _check_type(plt, [dict], "plotting parameters must be stored in dict")
        _check_keys(plt.keys(), PLOTTING_SUPPORTED_KEYS, _cc.PLOTTING)
        cartopy_flag = plt.get(_cc.PLOT_FLAG)
        if cartopy_flag is not None:
            _check_type(cartopy_flag, [bool], f"{_cc.PLOT_FLAG} must be of type bool")

            cp_params = plt.get(_cc.PLOT_PARAMS)
            if cp_params is not None:
                _check_type(cp_params, [dict], f"{_cc.PLOT_PARAMS} must be stored in dict")
                _check_keys(cp_params.keys(), CP_PARAM_SUPPORTED_KEYS, _cc.PLOT_PARAMS)
                scale = cp_params.get(_cc.PP_SCALE)
                if scale is not None:
                    _check_type(scale, [list, tuple], f"{_cc.PP_SCALE} must be of type list or tuple")
                    if len(scale) != 2:
                        raise ValueError(f"{_cc.PP_SCALE} must have exactly two elements")
                    for s in scale:
                        _check_type(s, [int, float], f"{_cc.PP_SCALE} values must be of type int or float")
                fmt = cp_params.get(_cc.PP_FMT)
                if fmt is not None:
                    _check_type(fmt, [str], 
                            f"{_cc.PP_FMT} must be desired output filetype extension as str")
                show = cp_params.get(_cc.PP_SHOW)
                if show is not None:
                    _check_type(show, [bool], f"{_cc.PP_SHOW} must be of type bool")

            gmt_flag = plt.get(_cc.GMT_FLAG)
            if gmt_flag is not None:
                _check_type(gmt_flag, [bool], f"{_cc.GMT_FLAG} must be of type bool")

            gis_flag = plt.get(_cc.GIS_FLAG)
            if gis_flag is not None:
                _check_type(gis_flag, [bool], f"{_cc.GIS_FLAG} must be of type bool")

    bounds = config.get(_cc.BOUNDS)
    if bounds is not None:
        _check_type(bounds, [list, tuple], f"{_cc.BOUNDS} must be of type of list or tuple")
        if len(bounds) != 4:
            raise ValueError(f"{_cc.BOUNDS} must have exactly 4 elements (west, south, east, north)")
        for b in bounds:
            _check_type(b, [float, int], "bounds values must be of type int or float")
        if bounds[0] >= bounds[2]:
            raise ValueError(f"Left bound ({bounds[0]}) must be less than right bound ({bounds[2]})")
        if bounds[1] >= bounds[3]:
            raise ValueError(f"Bottom bound ({bounds[1]}) must be less than top bound ({bounds[3]})")
    else:
        print("No bounds provided, bounds will be derived from extent of data")

    spacing = _try_lookup(config, _cc.GRID_INTERVAL, 
        f"grid cell size in decimal degrees is required as field {_cc.GRID_INTERVAL}")
    _check_type(spacing, [float, int], f"{_cc.GRID_INTERVAL} must be of value float or int")

    output_dir = config.get(_cc.OUTPUT_DIR)
    if output_dir is not None:
        _check_type(output_dir, [str], f"{_cc.OUTPUT_DIR} must be a path as type str")
    else:
        print("No output directory provided, current working directory will be used")

    interpolation = config.get(_cc.INTERPOLATION)
    if interpolation is not None:
        _check_type(interpolation, [str], f"{_cc.INTERPOLATION} must be of type str")
        if(interpolation not in ['bk_gaussian', 'bk_exponential', 'clough_tocher']):
            print("{_cc.INTERPOLATION} must be either 'bk_exponential', 'bk_gaussian' or clough_tocher")
    else:
        print(f"\n*** Parameter {_cc.INTERPOLATION} not found; the default value 'clough_tocher' will be used ***\n")

    weight_cutoff = config.get(_cc.WEIGHT_CUTOFF)
    if weight_cutoff is not None:
        _check_type(weight_cutoff, [float], f"{_cc.INTERPOLATION} must be of type float")
    else:
        if(interpolation != 'clough_tocher' and interpolation is not None):
            print(f"\n*** Parameter {_cc.WEIGHT_CUTOFF} not found; the default value 0 will be used ***\n")
# end func

class WorkflowParameters:
    def __init__(self, config_file):
        self.config_file = config_file
        with open(self.config_file) as f:
            config = json.load(f)

        validate(config)

        self.output_dir = config.get(_cc.OUTPUT_DIR, os.getcwd())

        self.bounds = config.get(_cc.BOUNDS)
        self.grid_interval = config[_cc.GRID_INTERVAL]
        self.grid_data = os.path.join(self.output_dir, _cc.MOHO_GRID)
        self.grad_data = os.path.join(self.output_dir, _cc.MOHO_GRAD)
        self.interpolation = config.get(_cc.INTERPOLATION, 'clough_tocher')
        self.weight_cutoff = config.get(_cc.WEIGHT_CUTOFF, 0)

        plotting = config.get(_cc.PLOTTING)
        if plotting is not None:
            self.plot_spatial_map = plotting.get(_cc.PLOT_FLAG, False)
            plot_params = plotting.get(_cc.PLOT_PARAMS)
            if plot_params is not None:
                self.plot_scale = plot_params.get(_cc.PP_SCALE)
                self.plot_fmt = plot_params.get(_cc.PP_FMT, 'png')
                self.plot_show = plot_params.get(_cc.PP_SHOW, False)
                self.plot_title = plot_params.get(_cc.PP_TITLE, 'Moho depth from blended data')
                self.plot_label = plot_params.get(_cc.PP_CB_LABEL, 'Moho depth (km)')
            self.plot_file = os.path.join(self.output_dir, _cc.MOHO_PLOT + f'.{self.plot_fmt}')

            self.gmt_write = plotting.get(_cc.GMT_FLAG, False)
            self.gmt_dir = os.path.join(self.output_dir, _cc.GMT_DIR)
            self.gmt_grid_file = os.path.join(self.gmt_dir, _cc.MOHO_GRID_GMT)
            self.gmt_grad_file = os.path.join(self.gmt_dir, _cc.MOHO_GRAD_GMT)
            self.gmt_loc_file = os.path.join(self.gmt_dir, '{}' + _cc.LOCATIONS_GMT)

            self.gis_write = plotting.get(_cc.GIS_FLAG, False)
            self.gis_dir = os.path.join(self.output_dir, _cc.GIS_DIR)
            self.gis_grid_file = os.path.join(self.gis_dir, _cc.MOHO_GRID_GIS)
            self.gis_grad_file = os.path.join(self.gis_dir, _cc.MOHO_GRAD_GIS)
            self.gis_loc_file = os.path.join(self.gis_dir, '{}' + _cc.LOCATIONS_GIS)

        # Correction/data prep
        data_prep = config.get(_cc.DATA_PREP)
        if data_prep is not None:
            # Correct d1 by d2
            for params in data_prep:
                d1 = MethodDataset({cc.DATA: params[cc.DATA_TO_PREP], cc.VAL_NAME: params[cc.DATA_VAL]})
                d2 = MethodDataset({cc.DATA: params[cc.CORR_DATA], cc.VAL_NAME: params[cc.CORR_VAL]})
                CORR_FUNC_MAP[params[cc.CORR_FUNC]](d1, d2, params[cc.CORR_OUT])
    
        # Filter out disabled methods
        methods = WorkflowParameters.filter_methods(config[_cc.METHODS])

        self.method_datasets = [MethodDataset(params) for params in methods]
        self.prioritise_samples()

    def prioritise_samples(self):
        """
        Sorts datasets by priority parameter, compares them from highest 
        to lowest. E.g. A.priority = 2, B.priority = 1, C.priority = 0 
        then comparison will be A vs [B, C], B vs C.

        The idea is that points in A close to those in B and C will
        supercede them (points in B and C will be removed) and again
        for points in B close to those in C.
        """
        pri_sorted = sorted(self.method_datasets, key=lambda x: x.priority, reverse=True)
        for i, p1 in enumerate(pri_sorted):
            for j, p2 in enumerate(pri_sorted):
                if -1 in [p1.priority, p2.priority]:
                    # Non-prioritised data, skip
                    continue
                if p1.priority > p2.priority:
                    if all(x is not None for x in [p1.net, p1.sta, p2.net, p2.sta]): 
                        # If network and station data is available, compare points by 
                        # station origin and distance
                        WorkflowParameters.prioritise_by_station(p1, p2)
                    else:
                        # Otherwise compare only by distance
                        WorkflowParameters.prioritise_by_location(p1, p2)
    
    @staticmethod
    def prioritise_by_station(m1, m2):
        """
        Prioritises samples by looking at points that share stations 
        between the two datasets. If two points in the datasets 
        share a station and are within 2km of each other, then 
        the higher priority dataset (m1) is prioritised (m2 point is
        removed).
        """
        to_remove = []
        for net, sta in set(zip(m1.net, m1.sta)):
            m1_station = np.where(np.logical_and(m1.net == net, m1.sta == sta))
            m2_station = np.where(np.logical_and(m2.net == net, m2.sta == sta))
            for lat_1, lon_1 in zip(m1.lat[m1_station], m1.lon[m1_station]):
                for i, (lat_2, lon_2) in enumerate(zip(m2.lat[m2_station], m2.lon[m2_station])):
                    dist = DIST_METRIC((lon_1, lat_1), (lon_2, lat_2))
                    if dist <= 0.018:
                        to_remove.append(m2_station[0][i])
        print(f'{len(to_remove)} samples from {m1.name} are being '
              f'prioritised over samples from {m2.name}')
        m2.remove(to_remove)

    @staticmethod
    def prioritise_by_location(m1, m2):
        """
        Prioritises samples. If two points in the datasets are
        within 2km of each other, then the higher priority dataset
        (m1) is prioritised (m2 point is removed).
        """
        to_remove = []
        for lat_1, lon_1 in zip(m1.lat, m1.lon):
            for i, (lat_2, lon_2) in enumerate(zip(m2.lat, m2.lon)):
                dist = DIST_METRIC((lon_1, lat_1), (lon_2, lat_2))
                if dist <= 0.018:
                    to_remove.append(i)
        print(f'{len(to_remove)} samples from {m1.name} are being '
              f'prioritised over samples from {m2.name}')
        m2.remove(to_remove)

    @staticmethod
    def filter_methods(methods):
        disabled_methods = [param[_cc.NAME] for param in methods if param.get(_cc.DISABLE, False)]
        if disabled_methods:
            print(f"Disabled methods: {disabled_methods}")
            methods = [param for param in methods if not param.get(_cc.DISABLE, False)]
            # Raise an error if all methods are disabled
            if not methods:
                raise ValueError("All methods have been disabled - there is no data to process!")
        return methods


class MethodDataset:
    start_names = ['Start', 'START', 'start']
    time_names = ['Time', 'TIME', 'time']
    net_names = ['net', 'Net', 'NET', 'network', 'Network', 'NETWORK']
    sta_names = ['sta', 'Sta', 'STA', 'station', 'Station', 'STATION']
    lat_names = ['lat', 'Lat', 'LAT', 'latitude', 'Latitude', 'LATITUDE']
    lon_names = ['lon', 'Lon', 'LON', 'longitude', 'Longitude', 'LONGITUDE']

    def __init__(self, method_params):
        self.name = method_params.get(_cc.NAME, 'undefined')
        self.scale_length = method_params.get(_cc.SCALE_LENGTH, None)
        self.priority = method_params.get(_cc.PRIORITY, -1)

        # Parse data file
        # Scan for header
        with open(method_params[_cc.DATA], 'r') as f:
            start_line = 0
            start_found = False
            while not start_found:
                start_line += 1
                line = f.readline()
                if line.startswith('#'):
                    if line.strip('#, \n') in MethodDataset.start_names:
                        start_found = True
                if line == '' and not start_found:
                    raise Exception(f"No 'START' flag found in data file {method_params[_cc.DATA]}")

            # Check if next line after start flag contains epoch time
            line = f.readline()
            line_parts = line.split(' ')
            # It's the header without any spaces
            if len(line_parts) == 1:
                self.datetime = None
                header_line = line
                data_start = start_line + 1
            elif line_parts[1] in MethodDataset.time_names:
                time = line_parts[2].strip('#, \n')
                # Need to pass epoch as int to datetime64, find better way to differentiate ISO/epoch
                # If time contains 'T', treat as ISO
                if 'T' in time:
                    self.datetime = np.datetime64(time, 's')
                else:
                    self.datetime = np.datetime64(int(time), 's')
                header_line = f.readline()
                data_start = start_line + 2
            else:
                self.datetime = None
                header_line = line
                data_start = start_line + 1
            header_parts = header_line.split(',')
            net_col, sta_col, lat_col, lon_col, val_col, sw_col = -1, -1, -1, -1, -1, -1
            col_names = []
            usecols = []
            for i, x in enumerate(header_parts):
                x = x.strip('#, \n')
                if x in MethodDataset.net_names:
                    net_col = i
                    col_names.append('net')
                    usecols.append(net_col)
                elif x in MethodDataset.sta_names:
                    sta_col = i
                    col_names.append('sta')
                    usecols.append(sta_col)
                elif x in MethodDataset.lat_names:
                    lat_col = i
                    col_names.append('lat')
                    usecols.append(lat_col)
                elif x in MethodDataset.lon_names:
                    lon_col = i
                    col_names.append('lon')
                    usecols.append(lon_col)
                elif x == method_params.get(_cc.VAL_NAME):
                    val_col = i
                    col_names.append('val')
                    usecols.append(val_col)
                elif x == method_params.get(_cc.SW_NAME):
                    sw_col = i
                    col_names.append('sw')
                    usecols.append(sw_col)

            if val_col == -1:
                raise Exception(f"No value column named '{_cc.VAL_NAME}' found in data file")

            if (lat_col == -1 or lon_col == -1) and (sta_col == -1 or net_col == -1):
                raise Exception(f"Neither lat/lon columns or net/sta columns were provided. "
                        "Either lat/lon needs to specified for each sample or net/sta and an "
                        "inventory file")
            
            data = np.genfromtxt(method_params[_cc.DATA], delimiter=',', dtype=None,
                    names=col_names, usecols=usecols, encoding=None, skip_header=data_start)

            # Filter out NaN values
            data = data[~np.isnan(data['val'])]

            self.val = data['val']

            if 'sw' in col_names:
                self.sw = data['sw']
            else:
                self.sw = np.ones_like(self.val)

            self.dataset_weight = method_params.get(_cc.WEIGHT, 1.0)
            self.total_weight = self.sw * self.dataset_weight

            if 'net' in col_names:
                self.net = data['net']
            else:
                self.net = None

            if 'sta' in col_names:
                self.sta = data['sta']
            else:
                self.sta = None

            if 'lat' in col_names and 'lon' in col_names:
                self.lat = data['lat']
                self.lon = data['lon']
            # Derive locations from inventory file
            else:   
                # Only have access to the Net.Sta code so we pick the first channel for this
                # station and take that location.
                lats, lons = [], []
                cache = {}
                try:
                    inv = read_inventory(method_params[_cc.INV_FILE])
                except IOError as e:
                    raise Exception("Inventory file doesn't exist") from e
                for net, sta in zip(self.net, self.sta):
                    coords = cache.get('.'.join((net, sta)))
                    if coords is None:
                        for c in inv.get_contents()['channels']:
                            if c.startswith('.'.join((net,sta))):
                                coords = inv.get_coordinates(c)
                                lats.append(float(coords['latitude']))
                                lons.append(float(coords['longitude']))
                                cache['.'.join((net,sta))] = (float(coords['latitude']), 
                                    float(coords['longitude']))
                                break
                    else:
                        lats.append(coords[0])
                        lons.append(coords[1])

                self.lat = np.array(lats)
                self.lon = np.array(lons)

                if self.lat.size == 0 or self.lon.size == 0:
                    raise ValueError(
                        f"Could't find any location metadata in inventory file for dataset "
                        f"{self.name}. Make sure you are using the correct inventory file and that "
                        f"the network and station codes in the CSV file are as they appear in the "
                        f"inventory. Exiting as we can't continue without location data.")
            
    def remove(self, indices):
        if not indices:
            return
        else:
            self.val = np.delete(self.val, indices)
            self.sw = np.delete(self.sw, indices)
            self.total_weight = self.sw * self.dataset_weight
            self.net = np.delete(self.net, indices) if self.net is not None else None
            self.sta = np.delete(self.sta, indices) if self.sta is not None else None
            self.lat = np.delete(self.lat, indices)
            self.lon = np.delete(self.lon, indices)

