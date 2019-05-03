#!/bin/env python
"""
Bulk analysis script for analysing relative traveltime residuals from a pick ensemble
for the purpose of identifying time periods of GPS clock error in specific stations.

Example usage, which plots 7X.MA11 and 7X.MA12 residuals relative to all common events on
AU network:

``relative_tt_residuals_plotter.py --network1=AU --networks2=7X --stations2="MA11,MA12" /c/data_cache/Picks/20190320/ensemble.p.txt``

"""

import os
import datetime
import numpy as np
import pandas as pd
from pandas.plotting import register_matplotlib_converters

import pytz
import matplotlib
import matplotlib.dates
import matplotlib.pyplot as plt

import click
# Progress bar helper to indicate that slow tasks have not stalled
import tqdm

import obspy

# pylint: disable=invalid-name, fixme, too-many-locals, too-many-statements
# pylint: disable=attribute-defined-outside-init

register_matplotlib_converters()

# Priority order of trusted channels
CHANNEL_PREF_NO_SHZ = ['HHZ', 'HHZ_10', 'H?Z', 'BHZ_00', 'BHZ', 'BHZ_10', 'B?Z']
CHANNEL_PREF_BALANCED = CHANNEL_PREF_NO_SHZ + ['S?Z', 'SHZ']
CHANNEL_PREF_GREEDY = CHANNEL_PREF_BALANCED + ['???', '?']
CHANNEL_PREF = CHANNEL_PREF_BALANCED

# Default filter options
DEFAULT_MIN_DISTANCE = 30.0
DEFAULT_MAX_DISTANCE = 90.0
DEFAULT_MIN_EVENT_SNR = 10
DEFAULT_CWT_CUTOFF = 15
DEFAULT_SLOPE_CUTOFF = 3
DEFAULT_NSIGMA_CUTOFF = 4
DEFAULT_MIN_EVENT_MAG = 5.5
DEFAULT_STRICT_FILTERING = True

# IRIS format station catalogs providing additional station metadata to merge with the target network.
# Useful for pulling in additional results from comprehensive networks such as GE, IR and IU.
IRIS_ALTERNATE_STATIONS_FILE = {
    "AU": "AU_irisws-fedcatalog_20190305T012747Z.txt"
    }

DATE_FILTER_TEMP_DEPLOYMENTS = True

TEMP_DEPLOYMENTS = {}  # TODO: Remove from global scope


class FilterOptions:
    """Simple container type for filtering options.
    """
    def __init(self):
        self.strict_filtering = DEFAULT_STRICT_FILTERING
        self.min_event_snr = DEFAULT_MIN_EVENT_SNR
        self.cwt_cutoff = DEFAULT_CWT_CUTOFF
        self.slope_cutoff = DEFAULT_SLOPE_CUTOFF
        self.nsigma_cutoff = DEFAULT_NSIGMA_CUTOFF
        self.min_event_mag = DEFAULT_MIN_EVENT_MAG


class DisplayOptions:
    """Simple container type for display options.
    """
    def __init(self):
        self.show_deployments = False
        # Historical events to add to the plot
        self.events = None


class BatchOptions:
    """
    Simple container type for run time options.
    """
    def __init__(self):
        self.save_file = True
        # Setting this to False will generate a lot more events per station chart, sometimes making it
        # easier to spot drift. But it may also add many events with significant non-zero residual.
        self.batch_label = ''
        # X-axis time range
        self.x_range = None
        # Path to folder in which to save the clock errors to csv file
        self.export_path = None


def get_network_stations(df, netcode):
    """
    Get the unique station codes belonging to a given network from a Pandas DataFrame

    :param df: Pandas dataframe loaded from a pick event ensemble
    :type df: pandas.DataFrame
    :param netcode: Network code for which station codes will be returned
    :type netcode: str
    :return: Sorted list of station code strings extracted from df
    :rtype: list(str)
    """
    return sorted(df[df['net'] == netcode]['sta'].unique().tolist())


def get_network_mean(df, netcode):
    """
    Get the mean station latitude and longitude coordinates for all stations in a given network.

    :param df: Pandas dataframe loaded from a pick event ensemble
    :type df: pandas.DataFrame
    :param netcode: Network code for which mean coordinates will be returned
    :type netcode: str
    :return: Mean (latitude, longitude) coordinates of stations in the network
    :rtype: tuple(float, float)
    """
    mean_lat = df[df['net'] == netcode]['stationLat'].mean()
    mean_lon = df[df['net'] == netcode]['stationLon'].mean()
    return (mean_lat, mean_lon)


def get_network_date_range(df, netcode):
    """
    Get the date range of pick events in df for a given network code

    :param df: Pandas dataframe loaded from a pick event ensemble
    :type df: pandas.DataFrame
    :param netcode: Network code whose pick event dates min and max will be returned
    :type netcode: str
    :return: Min and max dates of picks for given network
    :rtype: tuple(obspy.UTCDateTime, obspy.UTCDateTime)
    """
    mask = (df['net'] == netcode)
    df_net = df.loc[mask]
    min_date = df_net['originTimestamp'].min()
    max_date = df_net['originTimestamp'].max()
    return (obspy.UTCDateTime(min_date), obspy.UTCDateTime(max_date))


def get_station_date_range(df, netcode, statcode):
    """
    Get the date range of pick events in df for a given network and station code

    :param df: Pandas dataframe loaded from a pick event ensemble
    :type df: pandas.DataFrame
    :param netcode: Network code
    :type netcode: str
    :param statcode: Station code
    :type statcode: str
    :return: Min and max dates of picks for given network and station
    :rtype: tuple(obspy.UTCDateTime, obspy.UTCDateTime)
    """
    mask = (df['net'] == netcode)
    df_net = df.loc[mask]
    mask = (df_net['sta'] == statcode)
    df_sta = df_net.loc[mask]
    min_date = df_sta['originTimestamp'].min()
    max_date = df_sta['originTimestamp'].max()
    return (obspy.UTCDateTime(min_date), obspy.UTCDateTime(max_date))


def get_overlapping_date_range(df, ref_network, target_network):
    """
    Get the range of dates for which pick events in df from ref_network overlap with
    any picks from target_network

    :param df: Pandas dataframe loaded from a pick event ensemble
    :type df: pandas.DataFrame
    :param ref_network: Network and station codes for reference network
    :type ref_network: dict of corresponding network and station codes under keys 'net' and 'sta'
    :param target_network: Network and station codes for target network
    :type target_network: dict of corresponding network and station codes under keys 'net' and 'sta'
    :return: Start and end dates of datetime range during which pick events overlap for
        ref_network and target_network
    :rtype: tuple(obspy.UTCDateTime, obspy.UTCDateTime)
    """
    mask_ref = df[list(ref_network)].isin(ref_network).all(axis=1)
    mask_targ = df[list(target_network)].isin(target_network).all(axis=1)
    mask = mask_ref | mask_targ
    if not np.any(mask):
        return (None, None)
    df_nets = df.loc[mask]
    keep_events = [e for e, d in df_nets.groupby('#eventID')
                   if np.any(d[list(ref_network)].isin(ref_network).all(axis=1)) and
                   np.any(d[list(target_network)].isin(target_network).all(axis=1))]
    event_mask = df_nets['#eventID'].isin(keep_events)
    df_nets = df_nets[event_mask]
    return (obspy.UTCDateTime(df_nets['originTimestamp'].min()), obspy.UTCDateTime(df_nets['originTimestamp'].max()))


def get_iris_station_codes(src_file, original_network):
    """
    Extract the station codes for a given network code from a IRIS query result file.

    :param src_file: IRIS catalog query result file containing network and station information
    :type src_file: str or path
    :param original_network: Network code whose station codes will be extracted from IRIS file
    :type original_network: str
    :return: Pandas dataframe with station codes as the index and station latitude, longitude
        as columns
    :rtype: pandas.DataFrame
    """
    # Get station codes listed in IRIS whose network is original_network.
    # We need this to ensure we get complete coverage of chosen network, as it is
    # possible some such codes appear only under other network codes such as IR, GE, etc.. in the event catalog.
    # Returns a Pandas dataframe consisting of each station code and its mean (latitude, longitude) position.
    df = pd.read_csv(src_file, header=0, sep='|')
    df.columns = [c.strip() for c in df.columns.tolist()]
    orig_net_df = df.loc[(df['Network'] == original_network)]
    orig_net_df.columns = [c.strip() for c in orig_net_df.columns.tolist()]
    orig_net_perm_stations = sorted(orig_net_df['Station'].unique())
    mean_lat = []
    mean_lon = []
    for sta in orig_net_perm_stations:
        mean_lat.append(orig_net_df.loc[(orig_net_df['Station'] == sta), 'Latitude'].mean())
        std_dev = orig_net_df.loc[(orig_net_df['Station'] == sta), 'Latitude'].std(ddof=0)
        assert std_dev < 1.0, "{}: {}".format(sta, std_dev)
        mean_lon.append(orig_net_df.loc[(orig_net_df['Station'] == sta), 'Longitude'].mean())
        std_dev = orig_net_df.loc[(orig_net_df['Station'] == sta), 'Longitude'].std(ddof=0)
        assert std_dev < 1.0, "{}: {}".format(sta, std_dev)
    df_dict = {'sta': orig_net_perm_stations, 'lat': mean_lat, 'lon': mean_lon}
    result_df = pd.DataFrame(df_dict)
    return result_df.set_index(['sta'])


def determine_alternate_matching_codes(df, iris_file, original_network):
    """
    Find stations from other networks in df with the same station codes, but different
    network codes, whose positions match the stations of the same code in the original
    network.

    :param df: Pandas dataframe loaded from a pick event ensemble
    :type df: pandas.DataFrame
    :param iris_file: IRIS catalog query result file containing network and station information
    :type iris_file: str or path
    :param original_network: Network code whose station codes will be extracted from IRIS file
    :type original_network: str
    :return: Matching sequences of network and station codes
    :rtype: tuple(str), tuple(str)
    """
    # Get all the IRIS station codes for the given network
    matching_network_stn_iris_df = get_iris_station_codes(iris_file, original_network)

    # Generate mask for all station codes from IRIS that are present in df
    mask_iris_stns = df['sta'].isin(matching_network_stn_iris_df.index)
    # Generate mask for all records not belonging to given network
    mask_not_orig = (df['net'] != original_network)
    # Generate new dataframe of station codes from IRIS for given network that
    # are listed under a different network code.
    df_orig_stns_codes = df.loc[mask_iris_stns & mask_not_orig]

    # For each non-original network record, compute its distance from the known corresponding original station location
    from obspy.geodetics import locations2degrees

    def dist_to_orig_stn(row, orig_df):
        row_sta = row['sta']
        orig_df_sta = orig_df.loc[row_sta]
        return locations2degrees(row['stationLat'], row['stationLon'], orig_df_sta['lat'], orig_df_sta['lon'])

    # For each station code, compute its distance in degrees from the same station code location in IRIS.
    # If distance is too far away, station will be considered different from the AU network same station code.
    print("Computing distances to original network station locations...")
    distances_from_orig = df_orig_stns_codes.apply(lambda r: dist_to_orig_stn(r, matching_network_stn_iris_df), axis=1)

    # Only keep stations which are close to the same station code from original network.
    angle_tolerance_deg = 1.0
    df_orig_stns_codes_matching = df_orig_stns_codes.loc[(distances_from_orig < angle_tolerance_deg)]

    # Split out the network and matching station codes from the dataframe
    new_codes = [(n, s) for (n, s), _ in df_orig_stns_codes_matching.groupby(['net', 'sta'])]
    new_nets, new_stas = zip(*new_codes)
    return new_nets, new_stas


def pandas_timestamp_to_plottable_datetime(data):
    """
    Convert float UTC timestamp to equivalent type that is plottable by matplotlib

    :param data: Pandas series of float timestamps
    :type data: pandas.Series
    :return: Array of Python datetimes
    :rtype: numpy.array(datetime)
    """
    return data.transform(datetime.datetime.utcfromtimestamp).astype('datetime64[ms]').dt.to_pydatetime()


def _plot_target_network_rel_residuals(df, target, ref, batch_options, filter_options, tt_scale=50,
                                       snr_scale=(0, 60), annotator=None):

    file_label = batch_options.batch_label
    save_file = batch_options.save_file

    def _plot_dataset(ds, net_code, ref_code):
        # Sort ds rows by SNR, so that the weakest SNR points are drawn first and the high SNR point last,
        # to make sure high SNR point are in the top rendering layer.
        ds = ds.sort_values('snr')
        times = pandas_timestamp_to_plottable_datetime(ds['originTimestamp'])
        vals = ds[yaxis].values
        qual = ds['snr'].values
        min_mag = 4.0
        mag = ds['mag'].values - min_mag
        ylabel = 'Relative TT residual (sec)'
        title = r"Station {} TT residuals relative to network {} (filtering: ref SNR$\geq${}, CWT$\geq${}, "\
                r"slope$\geq${}, $n\sigma\geq{}$)".format(ref_code, net_code, str(filter_options.min_event_snr),
                                                          str(filter_options.cwt_cutoff),
                                                          str(filter_options.slope_cutoff),
                                                          str(filter_options.nsigma_cutoff))
        if vals.any():
            plt.figure(figsize=(32, 9))
            sc = plt.scatter(times, vals, c=qual, alpha=0.5, cmap='gnuplot_r', s=np.maximum(50 * mag, 10),
                             edgecolors=None, linewidths=0)
            time_formatter = matplotlib.dates.DateFormatter("%Y-%m-%d")
            plt.axes().xaxis.set_major_formatter(time_formatter)
            cb = plt.colorbar(sc, drawedges=False)
            cb.set_label('Signal to noise ratio', fontsize=12)
            plt.grid(color='#808080', linestyle=':', alpha=0.7)
            plt.xlabel(xlabel, fontsize=14)
            plt.ylabel(ylabel, fontsize=14)
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            plt.xlim(time_range)
            plt.ylim((-tt_scale, tt_scale))
            plt.clim(snr_scale)
            plt.title(title, fontsize=18)
            plt.legend(['Point size = Mag - {}, Color = SNR'.format(min_mag)], fontsize=12, loc=1)
            plt.text(0.01, 0.96, "Channel selection: {}".format(CHANNEL_PREF),
                     transform=plt.gca().transAxes, fontsize=12)
            plt.text(0.01, 0.92, "Start date: {}".format(str(time_range[0])),
                     transform=plt.gca().transAxes, fontsize=12)
            plt.text(0.01, 0.88, "  End date: {}".format(str(time_range[1])),
                     transform=plt.gca().transAxes, fontsize=12)
            plt.tight_layout(pad=1.05)
            if annotator is not None:
                annotator()
            if save_file:
                subfolder = os.path.join(ref_code.split('.')[0] + file_label, net_code)
                # subfolder = net_code
                os.makedirs(subfolder, exist_ok=True)
                plt_file = os.path.join(subfolder, '_'.join([ref_code, net_code]) + '_' +
                                        ylabel.replace(" ", "").replace(".*", "") + ".png")
                plt.savefig(plt_file, dpi=150)
                plt.close()
            else:
                plt.show()
    # end plot_dataset

    df_times = pandas_timestamp_to_plottable_datetime(df['originTimestamp'])
    if batch_options.x_range is not None:
        time_range = batch_options.x_range
    else:
        time_range = (df_times.min(), df_times.max())
    ref_code = ".".join([ref['net'][0], ref['sta'][0]])
    print("Plotting time range " + " to ".join([t.strftime("%Y-%m-%d %H:%M:%S") for t in time_range]) +
          " for " + ref_code)
    yaxis = 'relTtResidual'
    xlabel = 'Event Origin Timestamp'

    # Remove reference station from target set before producing composite image.
    # The reference station may not be there, but remove it if it is.
    mask_ref = df[list(ref)].isin(ref).all(axis=1)
    mask_targ = df[list(target)].isin(target).all(axis=1)
    df_agg = df[(mask_targ) & (~mask_ref)]
    _plot_dataset(df_agg, ','.join(np.unique(target['net'])), ref_code)

    # Export median error for each origin event if export path is provided
    if batch_options.export_path is not None:
        os.makedirs(batch_options.export_path, exist_ok=True)
        df_agg = df_agg.sort_values('originTimestamp')
        df_export = df_agg[['originTimestamp', 'relTtResidual']]
        median_errors = {'originTimestamp': [], 'rawClockError': []}
        for origin_ts, df_event in df_export.groupby('originTimestamp'):
            median_errors['originTimestamp'].append(origin_ts)
            median_errors['rawClockError'].append(df_event['relTtResidual'].median())
        df_export = pd.DataFrame(median_errors)
        fname = ref_code + "_raw_clock_error.csv"
        fname = os.path.join(batch_options.export_path, fname)
        df_export.to_csv(fname, index=False)


def plot_network_relative_to_ref_station(df_plot, ref, target_stns, batch_options, filter_options, display_options):
    """
    Compute relative residuals and send to plotting function.

    :param df_plot: Pandas dataframe containing only events which are common to ref station and target stations.
    :type df_plot: pandas.DataFrame
    :param ref_stn: Network and station codes for reference network (expected to be just one entry)
    :type ref_stn: dict of corresponding network and station codes under keys 'net' and 'sta'
         (expected to be just one entry)
    :param target_stns: Network and station codes for target network
    :type target_stns: dict of corresponding network and station codes under keys 'net' and 'sta'
    :param batch_options: Runtime options.
    :type batch_options: class BatchOptions
    :param filter_options: Filter options.
    :type filter_options: class FilterOptions
    :param display_options: Display options.
    :type display_options: class DisplayOptions
    """
    # Create column for entire table first
    df_plot['ttResidualRef'] = np.nan

    pbar = tqdm.tqdm(total=len(df_plot), ascii=True)
    pbar.set_description("Broadcasting REF STN residuals")
    # Assuming we only process one ref station at a time.
    df_plot['match_ref'] = (df_plot['net'] == ref['net'][0]) & (df_plot['sta'] == ref['sta'][0])
    for _, grp in df_plot.groupby('#eventID'):
        pbar.update(len(grp))
        ref_mask = grp['match_ref']
        grp_ref = grp[ref_mask]
        if grp_ref.empty:
            continue
        # Choose most favourable channel
        cha = None
        available_cha = grp_ref['cha'].values
        for c in CHANNEL_PREF:
            if c in available_cha:
                cha = c
                break
        # We must find a channel
        if cha is None:
            print("WARNING: Channels {} are not amongst allowed channels {}".format(available_cha, CHANNEL_PREF))
            continue
        cha_mask = (grp_ref['cha'] == cha)
        grp_cha = grp_ref[cha_mask]
        tt_ref_series = grp_cha['ttResidual'].unique()
        if len(tt_ref_series) > 1:
            # print("WARNING: Multiple reference times found for event {}\n{},"
            #       " choosing smallest absolute residual".format(eventid, grp_cha))
            # In this case, choose the smallest reference tt residual
            grp_cha['absTTResidual'] = np.abs(grp_cha['ttResidual'].values)
            grp_cha = grp_cha.sort_values('absTTResidual')
            tt_ref_series = grp_cha['ttResidual'].unique()
        ref_time = tt_ref_series[0]
        df_plot.loc[grp.index, 'ttResidualRef'] = ref_time
    pbar.close()

    # Quality check - each event should have only one unique reference tt residual
    assert np.all([len(df['ttResidualRef'].unique()) == 1 for e, df in df_plot.groupby('#eventID')])

    df_plot['relTtResidual'] = df_plot['ttResidual'] - df_plot['ttResidualRef']

    # Re-order columns
    df_plot = df_plot[['#eventID', 'originTimestamp', 'mag', 'originLon', 'originLat', 'originDepthKm',
                       'net', 'sta', 'cha', 'pickTimestamp', 'phase', 'stationLon', 'stationLat',
                       'distance', 'snr', 'ttResidual', 'ttResidualRef', 'relTtResidual',
                       'qualityMeasureCWT', 'qualityMeasureSlope', 'nSigma']]

    # Sort data by event origin time
    df_plot = df_plot.sort_values(['#eventID', 'originTimestamp'])

    def _plot_decorator(opts):
        if opts.events is not None:
            _add_event_marker_lines(opts.events)
        if opts.show_deployments:
            _add_temporary_deployment_intervals()

    _plot_target_network_rel_residuals(df_plot, target_stns, ref, batch_options, filter_options,
                                       annotator=lambda: _plot_decorator(display_options))


def _add_event_marker_lines(events):
    time_lims = plt.xlim()
    y_lims = plt.ylim()
    for date, event in events.iterrows():
        event_time = pytz.utc.localize(datetime.datetime.strptime(date, "%Y-%m-%d"))
        if event_time < matplotlib.dates.num2date(time_lims[0]) or \
           event_time >= matplotlib.dates.num2date(time_lims[1]):
            continue
        plt.axvline(event_time, linestyle='--', linewidth=1, color='#00800080')
        plt.text(event_time, y_lims[0] + 0.01 * (y_lims[1] - y_lims[0]), event['name'],
                 horizontalalignment='center', verticalalignment='bottom',
                 fontsize=12, fontstyle='italic', color='#008000c0', rotation=90)


def utc_time_string_to_plottable_datetime(utc_timestamp_str):
    """
    Convert a UTC timestamp string to datetime type that is plottable by matplotlib

    :param utc_timestamp_str: ISO-8601 UTC timestamp string
    :type utc_timestamp_str: str
    :return: Plottable datetime value
    :rtype: datetime.datetime with tzinfo
    """
    utc_time = obspy.UTCDateTime(utc_timestamp_str)
    return pytz.utc.localize(datetime.datetime.utcfromtimestamp(float(utc_time)))


def _add_temporary_deployment_intervals():
    """
    Graphical decorator for the TT residual charts to add visual indication of the time intervals
    during which selected temporary network deployments took place.
    """

    def _render_deployment_interval(rec_x, rec_y, rec_width, rec_height, rec_color, text):
        plt.gca().add_patch(plt.Rectangle((rec_x, rec_y), width=rec_width,
                                          height=rec_height, color=rec_color, alpha=0.8, clip_on=True))
        plt.text(rec_x, rec_y + rec_height / 2.0, text, fontsize=12, verticalalignment='center', clip_on=True)

    # Keep these in chronological order of start date to ensure optimal y-position staggering
    deployments = (TEMP_DEPLOYMENTS['7X'], TEMP_DEPLOYMENTS['7D'], TEMP_DEPLOYMENTS['7F'],
                   TEMP_DEPLOYMENTS['7G'], TEMP_DEPLOYMENTS['OA'])
    height = 5
    base_ypos = -40 - height / 2.0
    stagger = height / 2.0
    for d in deployments:
        ypos = base_ypos + stagger
        _render_deployment_interval(d[0], ypos, (d[1] - d[0]), height, d[2], d[3])
        stagger = -stagger


def apply_event_quality_filtering(df, ref_stn, filter_options):
    """
    Apply event quality requirements to pick events.

    :param df: Pandas dataframe loaded from a pick event ensemble
    :type df: pandas.DataFrame
    :param ref_stn: Network and station codes for reference network (expected to be just one entry)
    :type ref_stn: dict of corresponding network and station codes under keys 'net' and 'sta'
         (expected to be just one entry)
    :param filter_options: Filter options.
    :type filter_options: class FilterOptions
    :return: Filtered dataframe of pick events
    :rtype: pandas.DataFrame
    """
    # Remove records where the SNR is too low
    mask_snr = (df['snr'] >= filter_options.min_event_snr)

    # Filter to constrained quality metrics
    mask_cwt = (df['qualityMeasureCWT'] >= filter_options.cwt_cutoff)
    mask_slope = (df['qualityMeasureSlope'] >= filter_options.slope_cutoff)
    mask_sigma = (df['nSigma'] >= filter_options.nsigma_cutoff)

    # For events from ISC catalogs the quality metrics are zero, so we use event magnitude instead.
    mask_zero_quality_stats = (df[['snr', 'qualityMeasureCWT', 'qualityMeasureSlope', 'nSigma']] == 0).all(axis=1)
    mask_origin_mag = (df['mag'] >= filter_options.min_event_mag)

    quality_mask = (mask_snr & mask_cwt & mask_slope & mask_sigma) | (mask_zero_quality_stats & mask_origin_mag)

    mask_ref = df[list(ref_stn)].isin(ref_stn).all(axis=1)
    if filter_options.strict_filtering:
        # But never apply quality mask to ref stations that have all zero quality stats, as we just can't judge quality
        # and don't want to arbitrarily exclude them.
        quality_mask = (mask_zero_quality_stats & mask_ref) | quality_mask
    else:
        # Only apply quality mask to stations that are not the reference station, i.e. use all ref station events
        # regardless of pick quality at the ref station. This gives more results, but possibly more noise.
        quality_mask = mask_ref | (~mask_ref & quality_mask)

    assert np.sum(quality_mask) > 100, 'Not enough points left after quality filtering'
    df_qual = df[quality_mask]
    return df_qual


def analyze_target_relative_to_ref(df_picks, ref_stn, target_stns, batch_options, filter_options, display_options):
    """
    Analyze a single (reference) station's residuals relative to all the other stations
    in a (target) network.

    :param df_picks: Pandas dataframe loaded from a pick event ensemble
    :type df_picks: pandas.DataFrame
    :param ref_stn: Network and station codes for reference network (expected to be just one entry)
    :type ref_stn: dict of corresponding network and station codes under keys 'net' and 'sta'
         (expected to be just one entry)
    :param target_stns: Network and station codes for target network
    :type target_stns: dict of corresponding network and station codes under keys 'net' and 'sta'
    :param batch_options: Runtime options.
    :type batch_options: class BatchOptions
    :param filter_options: Filter options.
    :type filter_options: class FilterOptions
    :param display_options: Display options.
    :type display_options: class DisplayOptions
    """
    # Event quality filtering
    df_qual = apply_event_quality_filtering(df_picks, ref_stn, filter_options)

    # Filter to desired ref and target networks
    mask_ref = df_qual[list(ref_stn)].isin(ref_stn).all(axis=1)
    mask_targ = df_qual[list(target_stns)].isin(target_stns).all(axis=1)
    mask = mask_ref | mask_targ
    np.any(mask)
    df_nets = df_qual.loc[mask]

    # Filter out events in which ref_stn and TARGET stations are not both present
    print("Narrowing dataframe to events common to REF and TARGET networks...")
    df_nets['ref_match'] = df_nets[list(ref_stn)].isin(ref_stn).all(axis=1)
    df_nets['target_match'] = df_nets[list(target_stns)].isin(target_stns).all(axis=1)
    keep_events = [e for e, d in df_nets.groupby('#eventID') if np.any(d['ref_match']) and np.any(d['target_match'])]
    event_mask = df_nets['#eventID'].isin(keep_events)
    df_nets = df_nets[event_mask]
    print("Remaining picks after narrowing to common events: {}".format(len(df_nets)))
    if df_nets.empty:
        print("WARNING: no events left to analyze!")
        return

    # Alias for dataset at the end of all filtering, a static name that can be used from here onwards.
    ds_final = df_nets
    # print(getOverlappingDateRange(ds_final, ref_stn, target_stns))

    plot_network_relative_to_ref_station(ds_final, ref_stn, target_stns, batch_options, filter_options,
                                         display_options)


@click.command()
@click.argument('picks-file', type=click.Path('r'), required=True)
@click.option('--network1', type=str, required=True, help='Network against which to compute residuals')
@click.option('--stations1', type=str, required=False,
              help='Comma separated list of specific network 1 station codes to use.'
              'If empty, use all stations in network 1.')
@click.option('--networks2', type=str, required=True, help='Comma separated list of network(s) '
              'for which to compute residuals. Can be same as network1 to use itself as reference.')
@click.option('--stations2', type=str, required=False, help='Comma separated list of specific network 2 '
              'station codes to use. If empty, use all stations for each network in networks2. Should not '
              'be used unless networks2 has just a single network.')
@click.option('--min-distance', type=float, default=DEFAULT_MIN_DISTANCE, show_default=True,
              help='Minimum teleseismic distance (degrees) of events to plot')
@click.option('--max-distance', type=float, default=DEFAULT_MAX_DISTANCE, show_default=True,
              help='Maximum teleseismic distance (degrees) of events to plot')
@click.option('--min-event-snr', type=float, default=DEFAULT_MIN_EVENT_SNR, show_default=True,
              help='Filter out events with signal-to-noise ratio less than this')
@click.option('--cwt-cutoff', type=float, default=DEFAULT_CWT_CUTOFF, show_default=True,
              help='Filter out events with CWT quality value less than this')
@click.option('--slope-cutoff', type=float, default=DEFAULT_SLOPE_CUTOFF, show_default=True,
              help='Filter out events with slope quality value less than this')
@click.option('--nsigma-cutoff', type=int, default=DEFAULT_NSIGMA_CUTOFF, show_default=True,
              help='Filter out events with nsigma quality value less than this')
@click.option('--min-event-magnitude', type=float, default=DEFAULT_MIN_EVENT_MAG, show_default=True,
              help='Filter out events with magnitude value less than this. '
                   '(Only applies when pick quality stats are all zero.)')
@click.option('--strict-filtering/--no-strict-filtering', default=DEFAULT_STRICT_FILTERING, show_default=True,
              help='Whether to strictly apply quality filters to reference station.')
@click.option('--show-deployments', is_flag=True, default=False, show_default=True,
              help='Show temporary deployments time durations on the plots.')
@click.option('--show-historical/--no-show-historical', default=True, show_default=True,
              help='Show historical events on the plots.')
@click.option('--include-alternate-catalog/--no-include-alternate-catalog', default=True, show_default=True,
              help='Add matching stations from alternate networks in IRIS catalog')
@click.option('--export-path', type=click.Path(file_okay=False),
              help='Folder in which to store raw time series of TT residuals')
def main(picks_file, network1, networks2, stations1=None, stations2=None, 
         min_distance=DEFAULT_MIN_DISTANCE, max_distance=DEFAULT_MAX_DISTANCE,
         min_event_snr=DEFAULT_MIN_EVENT_SNR, cwt_cutoff=DEFAULT_CWT_CUTOFF,
         slope_cutoff=DEFAULT_SLOPE_CUTOFF, nsigma_cutoff=DEFAULT_NSIGMA_CUTOFF,
         min_event_magnitude=DEFAULT_MIN_EVENT_MAG, strict_filtering=DEFAULT_STRICT_FILTERING,
         show_deployments=False, show_historical=True, include_alternate_catalog=True,
         export_path=None):
    """
    Main function for running relative traveltime residual plotting. The picks ensemble file should
    have column headings:
    `#eventID originTimestamp mag originLon originLat originDepthKm net sta cha pickTimestamp phase stationLon stationLat az baz distance ttResidual snr qualityMeasureCWT domFreq qualityMeasureSlope bandIndex nSigma`

    :param picks_file: Picks ensemble file
    :type picks_file: str or path
    :param network1: Reference network code
    :type network1: str
    :param networks2: Comma separated list of network code(s) to analyze
    :type networks2: str
    :param stations1: Comma separated list of reference station codes to use for network1, defaults to None
    :param stations1: str, optional
    :param stations2: Comma separated list of station codes to analyze for network2, if only one network specified
        in networks2. Otherwise should not be used. Defaults to None.
    :param stations2: str, optional
    """

    print("Pandas version: " + pd.__version__)
    print("Matplotlib version: " + matplotlib.__version__)

    dtype = {'#eventID': object,
             'originTimestamp': np.float64,
             'mag': np.float64,
             'originLon': np.float64,
             'originLat': np.float64,
             'originDepthKm': np.float64,
             'net': object,
             'sta': object,
             'cha': object,
             'pickTimestamp': np.float64,
             'phase': object,
             'stationLon': np.float64,
             'stationLat': np.float64,
             'az': np.float64,
             'baz': np.float64,
             'distance': np.float64,
             'ttResidual': np.float64,
             'snr': np.float64,
             'qualityMeasureCWT': np.float64,
             'domFreq': np.float64,
             'qualityMeasureSlope': np.float64,
             'bandIndex': np.int64,
             'nSigma': np.int64}

    # Load in event ensemble
    print("Loading picks file {}".format(picks_file))
    df_raw_picks = pd.read_csv(picks_file, ' ', header=0, dtype=dtype)
    print("Number of raw picks: {}".format(len(df_raw_picks)))

    # Generate catalog of major regional events (mag 8+) for overlays
    if show_historical:
        df_mag8 = df_raw_picks[df_raw_picks['mag'] >= 8.0]
        df_mag8['day'] = df_mag8['originTimestamp'].transform(datetime.datetime.utcfromtimestamp)\
            .transform(lambda x: x.strftime("%Y-%m-%d"))
        df_mag8 = df_mag8.sort_values(['day', 'originTimestamp'])

        day_mag8_count = [(day, len(df_day)) for day, df_day in df_mag8.groupby('day')]
        dates, counts = zip(*day_mag8_count)
        mag8_dict = {'date': dates, 'counts': counts}
        mag8_events_df = pd.DataFrame(mag8_dict, columns=['date', 'counts'])

        event_count_threshold = 400
        significant_events = mag8_events_df[mag8_events_df['counts'] >= event_count_threshold]
        significant_events = significant_events.set_index('date')

        significant_events.loc['2001-06-23', 'name'] = '2001 South Peru Earthquake'
        significant_events.loc['2001-11-14', 'name'] = '2001 Kunlun earthquake'
        significant_events.loc['2002-11-03', 'name'] = '2002 Denali earthquake'
        significant_events.loc['2003-09-25', 'name'] = '2003 Tokachi-Oki earthquake'
        significant_events.loc['2004-12-26', 'name'] = '2004 Indian Ocean earthquake and tsunami'
        significant_events.loc['2005-03-28', 'name'] = '2005 Nias-Simeulue earthquake'
        significant_events.loc['2009-09-29', 'name'] = '2009 Samoa earthquake and tsunami'
        significant_events.loc['2010-02-27', 'name'] = '2010 Chile earthquake'
        significant_events.loc['2011-03-11', 'name'] = '2011 Tohoku earthquake and tsunami'
        significant_events.loc['2012-04-11', 'name'] = '2012 Indian Ocean earthquakes'
        significant_events.loc['2013-02-06', 'name'] = '2013 Solomon Islands earthquakes'
        significant_events.loc['2013-09-24', 'name'] = '2013 Balochistan earthquakes'
        significant_events.loc['2014-04-01', 'name'] = '2014 Iquique earthquake'
        significant_events.loc['2015-09-16', 'name'] = '2015 Illapel earthquake'
        significant_events.loc['2016-08-24', 'name'] = '2016 Myanmar earthquake'
    # end if

    # Populate temporary deployments details
    TEMP_DEPLOYMENTS['7X'] = (utc_time_string_to_plottable_datetime('2009-06-16T03:42:00.000000Z'),
                              utc_time_string_to_plottable_datetime('2011-04-01T23:18:49.000000Z'),
                              'C7', 'Deployment 7X')
    TEMP_DEPLOYMENTS['7D'] = (utc_time_string_to_plottable_datetime('2012-01-01T00:01:36.000000Z'),
                              utc_time_string_to_plottable_datetime('2014-03-27T15:09:51.000000Z'),
                              'C1', 'Deployment 7D')
    TEMP_DEPLOYMENTS['7F'] = (utc_time_string_to_plottable_datetime('2012-12-31T23:59:59.000000Z'),
                              utc_time_string_to_plottable_datetime('2014-11-15T00:43:14.000000Z'),
                              'C3', 'Deployment 7F')
    TEMP_DEPLOYMENTS['7G'] = (utc_time_string_to_plottable_datetime('2014-01-01T00:00:06.000000Z'),
                              utc_time_string_to_plottable_datetime('2016-02-09T21:04:29.000000Z'),
                              'C4', 'Deployment 7G')
    TEMP_DEPLOYMENTS['OA'] = (utc_time_string_to_plottable_datetime('2017-09-13T23:59:13.000000Z'),
                              utc_time_string_to_plottable_datetime('2018-11-28T01:11:14.000000Z'),
                              'C8', 'Deployment OA')

    # Query time period for source dataset
    print("Raw picks date range: {} to {}".format(obspy.UTCDateTime(df_raw_picks['originTimestamp'].min()),
                                                  obspy.UTCDateTime(df_raw_picks['originTimestamp'].max())))

    # Remove unwanted channels as their picks are not considered reliable enough to use
    df_picks = df_raw_picks[df_raw_picks['cha'].isin(CHANNEL_PREF)].reset_index()
    print("Remaining picks after channel filter: {}".format(len(df_picks)))

    # Remove unused columns
    df_picks = df_picks[['#eventID', 'originTimestamp', 'mag', 'originLon', 'originLat', 'originDepthKm',
                         'net', 'sta', 'cha', 'pickTimestamp', 'phase', 'stationLon', 'stationLat', 'az', 'baz',
                         'distance', 'ttResidual', 'snr', 'qualityMeasureCWT', 'qualityMeasureSlope', 'nSigma']]

    if DATE_FILTER_TEMP_DEPLOYMENTS:
        # Some select stations require custom date filters to remove events outside the known valid
        # date range of a network
        DATE_FILTER = (
            ('7D', pd.Timestamp(datetime.datetime(2010, 1, 1))),
            ('7G', pd.Timestamp(datetime.datetime(2010, 1, 1))),
        )
        before = len(df_picks)
        for net, min_date in DATE_FILTER:
            date_mask = (df_picks['net'] == net) & (df_picks['originTimestamp'] < min_date.timestamp())
            df_picks = df_picks[~date_mask]
        after = len(df_picks)
        print('Removed {} events due to known invalid timestamps'.format(before - after))

    # Filter to teleseismic events.
    # 'distance' is angular distance (degrees) between event and station.
    mask_tele = (df_picks['distance'] >= min_distance) & (df_picks['distance'] <= max_distance)
    df_picks = df_picks.loc[mask_tele]
    print("Remaining picks after filtering to teleseismic distances: {}".format(len(df_picks)))

    TARGET_NET = network1
    TARGET_STN = stations1.split(',') if stations1 else get_network_stations(df_picks, TARGET_NET)
    TARGET_STNS = {'net': [TARGET_NET] * len(TARGET_STN), 'sta': [s for s in TARGET_STN]}

    # Find additional network.station codes that match AU network, and add them to target
    if include_alternate_catalog and TARGET_NET in IRIS_ALTERNATE_STATIONS_FILE:
        alt_file = IRIS_ALTERNATE_STATIONS_FILE[TARGET_NET]
        new_nets, new_stas = determine_alternate_matching_codes(df_picks, alt_file, TARGET_NET)
        print("Adding {} more stations from alternate networks file {}".format(len(new_nets), alt_file))
        TARGET_STNS['net'].extend(list(new_nets))
        TARGET_STNS['sta'].extend(list(new_stas))

    REF_NETS = networks2.split(',')
    assert len(REF_NETS) == 1 or not stations2, \
        "Can't specify multiple networks and custom station list at the same time!"

    filter_options = FilterOptions()
    filter_options.strict_filtering = strict_filtering
    filter_options.min_event_snr = min_event_snr
    filter_options.min_event_mag = min_event_magnitude
    filter_options.cwt_cutoff = cwt_cutoff
    filter_options.slope_cutoff = slope_cutoff
    filter_options.nsigma_cutoff = nsigma_cutoff

    display_options = DisplayOptions()
    if show_historical:
        display_options.events = significant_events
    display_options.show_deployments = show_deployments

    batch_options = BatchOptions()
    batch_options.batch_label = '_strict' if strict_filtering else '_no_strict'
    batch_options.export_path = export_path
    for REF_NET in REF_NETS:
        if len(REF_NETS) == 1 and stations2:
            REF_STN = stations2.split(',')
        else:
            REF_STN = get_network_stations(df_picks, REF_NET)
        REF_STNS = {'net': [REF_NET] * len(REF_STN), 'sta': [s for s in REF_STN]}

        if REF_NET == 'AU' and TARGET_NET == 'AU':
            REF_STNS['net'].extend(list(new_nets))
            REF_STNS['sta'].extend(list(new_stas))

        # For certain temporary deployments, force a fixed time range on the x-axis so that all stations
        # in the deployment can be compared on a common time range.
        if REF_NET in TEMP_DEPLOYMENTS:
            batch_options.x_range = (TEMP_DEPLOYMENTS[REF_NET][0], TEMP_DEPLOYMENTS[REF_NET][1])
        else:
            batch_options.x_range = None

        for ref_net, ref_sta in zip(REF_STNS['net'], REF_STNS['sta']):
            print("Plotting against REF: " + ".".join([ref_net, ref_sta]))
            single_ref = {'net': [ref_net], 'sta': [ref_sta]}
            analyze_target_relative_to_ref(df_picks, single_ref, TARGET_STNS, batch_options, filter_options,
                                           display_options)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
