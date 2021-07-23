#!/bin/env python
"""
Bulk analysis script for analysing relative traveltime residuals from a pick ensemble
for the purpose of identifying time periods of GPS clock error in specific stations.

Example usage, which plots 7X.MA11 and 7X.MA12 residuals relative to all common events on
AU network::

    relative_tt_residuals_plotter.py --network1=AU --networks2=7X --stations2="MA11,MA12" /c/data_cache/Picks/20190320/ensemble.p.txt

"""

import os
import sys
import datetime
import logging

import numpy as np
import pandas as pd
from pandas.plotting import register_matplotlib_converters, deregister_matplotlib_converters
import pytz
import matplotlib
import matplotlib.dates
import matplotlib.pyplot as plt

import click
# Progress bar helper to indicate that slow tasks have not stalled
import tqdm
import obspy

if sys.version_info[0] < 3:
    import pathlib2 as pathlib  # pylint: disable=import-error
else:
    import pathlib  # pylint: disable=import-error

from seismic.gps_corrections.picks_reader_utils import (read_picks_ensemble, get_network_stations, compute_matching_network_mask, generate_large_events_catalog)

# pylint: disable=invalid-name, fixme, too-many-locals, too-many-statements
# pylint: disable=attribute-defined-outside-init, logging-format-interpolation, logging-not-lazy

logging.basicConfig()

# Priority order of trusted channels
CHANNEL_PREF_NO_SHZ = ['HHZ', 'HHZ_10', 'H?Z', 'BHZ_00', 'BHZ', 'BHZ_10', 'B?Z', 'EHZ']
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


class FilterOptions:
    """Simple container type for filtering options.
    """
    def __init__(self):
        # Setting strict_filtering to False will generate a lot more events per station chart, sometimes making it
        # easier to spot drift. But it may also add many events with significant non-zero residual (e.g. poor picks).
        self.strict_filtering = DEFAULT_STRICT_FILTERING
        self.min_event_snr = DEFAULT_MIN_EVENT_SNR
        self.cwt_cutoff = DEFAULT_CWT_CUTOFF
        self.slope_cutoff = DEFAULT_SLOPE_CUTOFF
        self.nsigma_cutoff = DEFAULT_NSIGMA_CUTOFF
        self.min_event_mag = DEFAULT_MIN_EVENT_MAG
        self.channel_preference = CHANNEL_PREF


class DisplayOptions:
    """Simple container type for display options.
    """
    def __init__(self):
        self.deployments = None
        # Historical events to add to the plot
        self.events = None


class BatchOptions:
    """
    Simple container type for run time options.
    """
    def __init__(self):  # pragma: no cover
        self.save_file = True
        self.batch_label = ''
        # X-axis time range
        self.x_range = None
        # Path to folder in which to save the clock errors to csv file
        self.export_path = None


def get_iris_station_codes(src_file, original_network):
    """
    Extract the station codes for a given network code from a IRIS query result file.

    :param src_file: IRIS catalog query result file containing network and station information, in stationtxt format
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
    :param iris_file: IRIS catalog query result file containing network and station information in stationtxt format
    :type iris_file: str or path
    :param original_network: Network code whose station codes will be extracted from IRIS file
    :type original_network: str
    :return: Matching sequences of network and station codes
    :rtype: tuple(str), tuple(str)
    """
    log = logging.getLogger(__name__)

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
    log.info("Computing distances to original network station locations...")
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
            plt.gca().xaxis.set_major_formatter(time_formatter)
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
            plt.text(0.01, 0.96, "Channel selection: {}".format(filter_options.channel_preference),
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
                pathlib.Path(subfolder).mkdir(exist_ok=True, parents=True)
                plt_file = os.path.join(subfolder, '_'.join([ref_code, net_code]) + '_' +
                                        ylabel.replace(" ", "").replace(".*", "") + ".png")
                plt.savefig(plt_file, dpi=150)
            else:  # pragma: no cover
                plt.show()
            # end if
            plt.close()
        else:
            log.warning("No values to plot for {}".format(ref_code))
        # end if
    # end plot_dataset

    df_times = pandas_timestamp_to_plottable_datetime(df['originTimestamp'])
    if batch_options.x_range is not None:
        time_range = batch_options.x_range
    else:
        time_range = (df_times.min(), df_times.max())
    ref_code = ".".join([ref['net'][0], ref['sta'][0]])
    log = logging.getLogger(__name__)
    log.info("Plotting time range " + " to ".join([t.strftime("%Y-%m-%d %H:%M:%S") for t in time_range]) +
             " for " + ref_code)
    yaxis = 'relTtResidual'
    xlabel = 'Event Origin Timestamp'

    # Remove reference station from target set before producing composite image.
    # The reference station may not be there, but remove it if it is.
    mask_ref = compute_matching_network_mask(df, ref)
    mask_targ = compute_matching_network_mask(df, target)
    df_agg = df[(mask_targ) & (~mask_ref)]
    _plot_dataset(df_agg, ','.join(np.unique(target['net'])), ref_code)

    # Export median error for each origin event if export path is provided
    if batch_options.export_path is not None:
        pathlib.Path(batch_options.export_path).mkdir(parents=True, exist_ok=True)
        df_agg = df_agg.sort_values('originTimestamp')
        df_export = df_agg[['originTimestamp', 'relTtResidual']]
        median_errors = {'originTimestamp': [], 'rawClockError': []}
        for origin_ts, df_event in df_export.groupby('originTimestamp'):
            median_errors['originTimestamp'].append(origin_ts)
            median_errors['rawClockError'].append(df_event['relTtResidual'].median())
        df_export = pd.DataFrame(median_errors)
        fname = ref_code + ".raw_clock_error.csv"
        fname = os.path.join(batch_options.export_path, fname)
        df_export.to_csv(fname, index=False)


def _plot_network_relative_to_ref_station(df_plot, ref, target_stns, batch_options, filter_options, display_options):
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
    register_matplotlib_converters()
    if df_plot is None:
        return
    df_plot = df_plot.assign(relTtResidual=(df_plot['ttResidual'] - df_plot['ttResidualRef']))

    # Re-order columns
    df_plot = df_plot[['#eventID', 'originTimestamp', 'mag', 'originLon', 'originLat', 'originDepthKm',
                       'net', 'sta', 'cha', 'pickTimestamp', 'stationLon', 'stationLat',
                       'distance', 'snr', 'ttResidual', 'ttResidualRef', 'relTtResidual',
                       'qualityMeasureCWT', 'qualityMeasureSlope', 'nSigma']]

    # Sort data by event origin time
    df_plot = df_plot.sort_values(['#eventID', 'originTimestamp'])

    def _plot_decorator(opts):
        if opts.events is not None:
            _add_event_marker_lines(opts.events)
        if opts.deployments is not None:
            _add_temporary_deployment_intervals(opts.deployments)

    _plot_target_network_rel_residuals(df_plot, target_stns, ref, batch_options, filter_options,
                                       annotator=lambda: _plot_decorator(display_options))
    deregister_matplotlib_converters()


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


def _add_temporary_deployment_intervals(deployments):
    """
    Graphical decorator for the TT residual charts to add visual indication of the time intervals
    during which selected temporary network deployments took place.
    """

    def _render_deployment_interval(rec_x, rec_y, rec_width, rec_height, rec_color, text):
        plt.gca().add_patch(plt.Rectangle((rec_x, rec_y), width=rec_width,
                                          height=rec_height, color=rec_color, alpha=0.8, clip_on=True))
        plt.text(rec_x, rec_y + rec_height / 2.0, text, fontsize=12, verticalalignment='center', clip_on=True)

    # Keep these in chronological order of start date to ensure optimal y-position staggering
    sorted_deployments = sorted(deployments, key=lambda v: v[0])
    height = 5
    base_ypos = -40 - height / 2.0
    stagger = height / 2.0
    for k in sorted_deployments:
        d = deployments[k]
        ypos = base_ypos + stagger
        _render_deployment_interval(d[0], ypos, (d[1] - d[0]), height, d[2], d[3])
        stagger = -stagger


def apply_event_quality_filtering(df, ref_stn, filter_options):
    """
    Apply event quality requirements to pick events. The event magnitude threshold is only used to
    filter rows where the quality metrics ['snr', 'qualityMeasureCWT', 'qualityMeasureSlope', 'nSigma']
    are all zero.

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

    if 'net' not in ref_stn:
        raise KeyError('Expect \'net\' in dict keys')
    if 'sta' not in ref_stn:
        raise KeyError('Expect \'sta\' in dict keys')

    if not len(ref_stn['net']) == len(ref_stn['sta']) == 1:
        raise ValueError('Expect only one network and station code when calling apply_event_quality_filtering')

    # Remove records where the SNR is too low
    mask_snr = (df['snr'] >= filter_options.min_event_snr)

    # Filter to constrained quality metrics
    mask_cwt = (df['qualityMeasureCWT'] >= filter_options.cwt_cutoff)
    mask_slope = (df['qualityMeasureSlope'] >= filter_options.slope_cutoff)
    mask_sigma = (df['nSigma'] >= filter_options.nsigma_cutoff)

    # For events from ISC catalogs the quality metrics are zero, so we use event magnitude instead.
    mask_zero_quality_stats = (df[['snr', 'qualityMeasureCWT', 'qualityMeasureSlope', 'nSigma']] == 0).all(axis=1)
    mask_origin_mag = (df['mag'] >= filter_options.min_event_mag)

    quality_mask = ((~mask_zero_quality_stats & mask_snr & mask_cwt & mask_slope & mask_sigma) |
                    (mask_zero_quality_stats & mask_origin_mag))

    mask_ref = compute_matching_network_mask(df, ref_stn)
    if filter_options.strict_filtering:
        # But never apply quality mask to ref stations that have all zero quality stats, as we just can't judge quality
        # and don't want to arbitrarily exclude them.
        quality_mask = (mask_zero_quality_stats & mask_ref) | quality_mask
    else:
        # Only apply quality mask to stations that are not the reference station, i.e. use all ref station events
        # regardless of pick quality at the ref station. This gives more results, but possibly more noise.
        quality_mask = mask_ref | (~mask_ref & quality_mask)

    df_qual = df[quality_mask]
    return df_qual


def broadcast_ref_residual_per_event(df_plot, ref_netcode, ref_stacode, filter_options):
    """For each event in the dataframe, figure out the best reference residual and add it to new column 'ttResidualRef'

    :param df: Pandas dataframe loaded from a pick event ensemble
    :type df: pandas.DataFrame
    :param ref_netcode: The reference network code
    :type ref_netcode: str
    :param ref_stacode: The reference station code
    :type ref_stacode: str
    :param filter_options: Filter options.
    :type filter_options: class FilterOptions
    :return: Dataframe with a consistent reference residual populated for each event for given reference station.
    :rtype: pandas.DataFrame
    """
    # Create column for entire table first
    df_plot['ttResidualRef'] = np.nan
    log = logging.getLogger(__name__)

    pbar = tqdm.tqdm(total=len(df_plot), ascii=True)
    pbar.set_description("Broadcasting {}.{} residuals".format(ref_netcode, ref_stacode))
    # Assuming we only process one ref station at a time.
    df_plot['match_ref'] = (df_plot['net'] == ref_netcode) & (df_plot['sta'] == ref_stacode)
    for eventid, grp in df_plot.groupby('#eventID'):
        pbar.update(len(grp))
        ref_mask = grp['match_ref']
        grp_ref = grp[ref_mask]
        if grp_ref.empty:
            continue
        # Choose most favourable channel
        cha = None
        available_cha = grp_ref['cha'].values
        for c in filter_options.channel_preference:
            if c in available_cha:
                cha = c
                break
        # We must find a channel
        if cha is None:
            log.warning("Channels {} are not amongst allowed channels {}".format(available_cha,
                                                                                 filter_options.channel_preference))
            continue
        cha_mask = (grp_ref['cha'] == cha)
        grp_cha = grp_ref[cha_mask]
        tt_ref_series = grp_cha['ttResidual'].unique()
        if len(tt_ref_series) > 1:
            log.warning("WARNING: Multiple reference times found for event {}\n{},"
                        " choosing smallest absolute residual".format(eventid, tt_ref_series))
            # In this case, choose the smallest reference tt residual
            grp_cha['absTTResidual'] = np.abs(grp_cha['ttResidual'].values)
            grp_cha = grp_cha.sort_values('absTTResidual')
            tt_ref_series = grp_cha['ttResidual'].unique()
        ref_time = tt_ref_series[0]
        assert not np.isnan(ref_time)
        df_plot.loc[grp.index, 'ttResidualRef'] = ref_time
    pbar.close()

    # Quality check - each event should have only one unique reference tt residual
    assert np.all([len(df['ttResidualRef'].unique()) == 1 for e, df in df_plot.groupby('#eventID')])

    return df_plot


def analyze_target_relative_to_ref(df_picks, ref_stn, target_stns, filter_options):
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
    :param filter_options: Filter options.
    :type filter_options: class FilterOptions
    :return: Dataframe with filtering applied and reference residual populated ready for plotting.
    :rtype: pandas.DataFrame
    """
    assert len(ref_stn['net']) == 1
    assert len(ref_stn['sta']) == 1
    # Event quality filtering
    df_qual = apply_event_quality_filtering(df_picks, ref_stn, filter_options)

    # Filter to desired ref and target networks
    mask_ref = compute_matching_network_mask(df_qual, ref_stn)
    mask_targ = compute_matching_network_mask(df_qual, target_stns)
    mask = mask_ref | mask_targ
    # np.any(mask)
    df_nets = df_qual.loc[mask]

    # Filter out events in which ref_stn and TARGET stations are not both present
    log = logging.getLogger(__name__)
    log.info("Narrowing dataframe to events common to REF and TARGET networks...")
    df_nets = df_nets.assign(ref_match=compute_matching_network_mask(df_nets, ref_stn))
    df_nets = df_nets.assign(target_match=compute_matching_network_mask(df_nets, target_stns))
    keep_events = [e for e, d in df_nets.groupby('#eventID') if np.any(d['ref_match']) and np.any(d['target_match'])]
    event_mask = df_nets['#eventID'].isin(keep_events)
    df_nets = df_nets[event_mask]
    log.info("Remaining picks after narrowing to common events: {}".format(len(df_nets)))
    if df_nets.empty:
        log.warning("No events left to analyze!")
        return
    # print(get_overlapping_date_range(df_nets, ref_stn, target_stns))

    # Compute the reference TT residual for each event and store in new column
    df_final = broadcast_ref_residual_per_event(df_nets, ref_stn['net'][0], ref_stn['sta'][0], filter_options)
    assert 'ttResidualRef' in df_final.columns

    return df_final


def filter_limit_channels(df_picks, channel_pref):
    """Filter picks dataframe to a limited range of preferred channel types.

    :param df_picks: Picks dataframe to filter
    :type df_picks: pandas.DataFrame
    :return: Filtered picks dataframe with only preferred channels
    :rtype: pandas.DataFrame
    """
    df_picks = df_picks[df_picks['cha'].isin(channel_pref)].reset_index()
    return df_picks


def filter_duplicated_network_codes(df_picks):
    """Filter picks dataframe to remove records for selected known duplicate network codes based on date
    of the records we want to keep.

    :param df_picks: Picks dataframe to filter
    :type df_picks: pandas.DataFrame
    :return: Filtered picks dataframe with unwanted records removed.
    :rtype: pandas.DataFrame
    """
    DATE_FILTER = (
        ('7D', pd.Timestamp(datetime.datetime(2010, 1, 1))),
        ('7G', pd.Timestamp(datetime.datetime(2010, 1, 1))),
    )
    for net, min_date in DATE_FILTER:
        date_mask = (df_picks['net'] == net) & (df_picks['originTimestamp'] < min_date.timestamp())
        df_picks = df_picks[~date_mask]
    return df_picks


def filter_to_teleseismic(df_picks, min_dist_deg, max_dist_deg):
    """Filter picks dataframe to limited range of teleseismic distances of the original events.

    :param df_picks: Picks dataframe to filter
    :type df_picks: pandas.DataFrame
    :param min_dist_deg: Minimum teleseismic distance (degrees)
    :type min_dist_deg: float
    :param max_dist_deg: Maximum teleseismic distance (degrees)
    :type max_dist_deg: float
    :return: Filtered picks dataframe containing only events with limited range of teleseismic distance.
    :rtype: pandas.DataFrame
    """
    mask_tele = (df_picks['distance'] >= min_dist_deg) & (df_picks['distance'] <= max_dist_deg)
    df_picks = df_picks.loc[mask_tele]
    return df_picks


def _get_known_temporary_deployments():
    # TODO: Change implementation to read this data from JSON file
    temp_deps = {}
    # Data fields are (start date, end date, plot color, plot label)
    temp_deps['7X'] = (utc_time_string_to_plottable_datetime('2009-06-16T03:42:00.000000Z'),
                       utc_time_string_to_plottable_datetime('2011-04-01T23:18:49.000000Z'),
                       'C7', 'Deployment 7X')
    temp_deps['7D'] = (utc_time_string_to_plottable_datetime('2012-01-01T00:01:36.000000Z'),
                       utc_time_string_to_plottable_datetime('2014-03-27T15:09:51.000000Z'),
                       'C1', 'Deployment 7D')
    temp_deps['7F'] = (utc_time_string_to_plottable_datetime('2012-12-31T23:59:59.000000Z'),
                       utc_time_string_to_plottable_datetime('2014-11-15T00:43:14.000000Z'),
                       'C3', 'Deployment 7F')
    temp_deps['7G'] = (utc_time_string_to_plottable_datetime('2014-01-01T00:00:06.000000Z'),
                       utc_time_string_to_plottable_datetime('2016-02-09T21:04:29.000000Z'),
                       'C4', 'Deployment 7G')
    temp_deps['OA'] = (utc_time_string_to_plottable_datetime('2017-09-13T23:59:13.000000Z'),
                       utc_time_string_to_plottable_datetime('2018-11-28T01:11:14.000000Z'),
                       'C8', 'Deployment OA')
    return temp_deps


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
              help='Whether to apply quality filters to picks of the station being analyzed/plotted.'
                   'If False, quality filters are still applied to the reference station picks, but'
                   'not to the target station, resulting in a larger number of plotted events per station.')
@click.option('--show-deployments', is_flag=True, default=False, show_default=True,
              help='Show temporary deployments time durations on the plots.')
@click.option('--show-historical/--no-show-historical', default=True, show_default=True,
              help='Show historical events on the plots.')
@click.option('--include-alternate-catalog/--no-include-alternate-catalog', default=True, show_default=True,
              help='Add matching stations from alternate networks in IRIS catalog')
@click.option('--export-path', type=click.Path(file_okay=False),
              help='Folder in which to store raw time series of TT residuals')
@click.option('--interactive', is_flag=True, default=False, show_default=True,
              help='If True, plots will be displayed as popup windows instead of saving to file. '
                   'Use this option to interact with the data.')
def main(picks_file, network1, networks2, stations1=None, stations2=None,
         min_distance=DEFAULT_MIN_DISTANCE, max_distance=DEFAULT_MAX_DISTANCE,
         min_event_snr=DEFAULT_MIN_EVENT_SNR, cwt_cutoff=DEFAULT_CWT_CUTOFF,
         slope_cutoff=DEFAULT_SLOPE_CUTOFF, nsigma_cutoff=DEFAULT_NSIGMA_CUTOFF,
         min_event_magnitude=DEFAULT_MIN_EVENT_MAG, strict_filtering=DEFAULT_STRICT_FILTERING,
         show_deployments=False, show_historical=True, include_alternate_catalog=True,
         export_path=None, interactive=False):  # pragma: no cover
    """
    Main function for running relative traveltime residual plotting. The picks ensemble file should
    have column headings in accordance with picks_reader_utils.PICKS_TABLE_SCHEMA.

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
    # TODO: invert naming of "target" and "reference" in the code. The main help documentation is correct,
    #       but internal naming is actually swapped.

    log = logging.getLogger(__name__)
    # log.setLevel(logging.INFO)
    log.setLevel(logging.DEBUG)

    log.info("Pandas version: " + pd.__version__)
    log.info("Matplotlib version: " + matplotlib.__version__)

    # Load in event ensemble
    log.info("Loading picks file {}".format(picks_file))
    df_raw_picks = read_picks_ensemble(picks_file)
    # Remove unused columns
    df_raw_picks = df_raw_picks[['#eventID', 'originTimestamp', 'mag', 'originLon', 'originLat', 'originDepthKm',
                                 'net', 'sta', 'cha', 'pickTimestamp', 'stationLon', 'stationLat', 'az',
                                 'baz', 'distance', 'ttResidual', 'snr', 'qualityMeasureCWT', 'qualityMeasureSlope',
                                 'nSigma']]
    log.debug("Number of raw picks: {}".format(len(df_raw_picks)))
    # Log time period for source dataset
    log.debug("Raw picks date range: {} to {}".format(obspy.UTCDateTime(df_raw_picks['originTimestamp'].min()),
                                                      obspy.UTCDateTime(df_raw_picks['originTimestamp'].max())))

    filter_options = FilterOptions()
    filter_options.strict_filtering = strict_filtering
    filter_options.min_event_snr = min_event_snr
    filter_options.min_event_mag = min_event_magnitude
    filter_options.cwt_cutoff = cwt_cutoff
    filter_options.slope_cutoff = slope_cutoff
    filter_options.nsigma_cutoff = nsigma_cutoff

    # Perform series of global filters on picks
    df_picks = df_raw_picks
    log.debug("Picks before global filter: {}".format(len(df_picks)))
    filter_chain = [lambda _df: filter_limit_channels(_df, filter_options.channel_preference),
                    filter_duplicated_network_codes,
                    lambda _df: filter_to_teleseismic(_df, min_distance, max_distance)]
    for f in filter_chain:
        df_picks = f(df_picks)
    log.debug("Picks after global filter: {}".format(len(df_picks)))

    # Populate temporary deployments details
    temporary_deployments = _get_known_temporary_deployments()

    TARGET_NET = network1
    TARGET_STN = stations1.split(',') if stations1 else get_network_stations(df_picks, TARGET_NET)
    TARGET_STNS = {'net': [TARGET_NET] * len(TARGET_STN), 'sta': [s for s in TARGET_STN]}

    # Find additional network.station codes that match AU network, and add them to target
    if include_alternate_catalog and TARGET_NET in IRIS_ALTERNATE_STATIONS_FILE:
        alt_file = IRIS_ALTERNATE_STATIONS_FILE[TARGET_NET]
        alt_file = os.path.join(os.path.split(__file__)[0], alt_file)
        new_nets, new_stas = determine_alternate_matching_codes(df_picks, alt_file, TARGET_NET)
        log.info("Adding {} more stations from alternate networks file {}".format(len(new_nets), alt_file))
        TARGET_STNS['net'].extend(list(new_nets))
        TARGET_STNS['sta'].extend(list(new_stas))

    REF_NETS = networks2.split(',')
    assert len(REF_NETS) == 1 or not stations2, \
        "Can't specify multiple networks and custom station list at the same time!"

    display_options = DisplayOptions()
    if show_historical:
        # Generate catalog of major regional events (mag 8+) for overlays
        display_options.events = generate_large_events_catalog(df_raw_picks, 8.0)
    if show_deployments:
        display_options.deployments = temporary_deployments

    batch_options = BatchOptions()
    batch_options.save_file = not interactive
    batch_options.batch_label = '_strict' if strict_filtering else '_no_strict'
    batch_options.export_path = export_path
    assert len(TARGET_STNS['net']) == len(TARGET_STNS['sta'])
    for ref_net in REF_NETS:
        if len(REF_NETS) == 1 and stations2:
            REF_STN = stations2.split(',')
        else:
            REF_STN = get_network_stations(df_picks, ref_net)
        # end if
        REF_STNS = {'net': [ref_net] * len(REF_STN), 'sta': [s for s in REF_STN]}

        if ref_net == 'AU' and TARGET_NET == 'AU':
            REF_STNS['net'].extend(list(new_nets))
            REF_STNS['sta'].extend(list(new_stas))

        # For certain temporary deployments, force a fixed time range on the x-axis so that all stations
        # in the deployment can be compared on a common time range.
        if ref_net in temporary_deployments:
            batch_options.x_range = (temporary_deployments[ref_net][0], temporary_deployments[ref_net][1])
        else:
            batch_options.x_range = None

        assert len(REF_STNS['net']) == len(REF_STNS['sta'])
        for net, sta in zip(REF_STNS['net'], REF_STNS['sta']):
            log.info("Plotting against REF: " + ".".join([net, sta]))
            single_ref = {'net': [net], 'sta': [sta]}
            # Peform analysis
            df_final = analyze_target_relative_to_ref(df_picks, single_ref, TARGET_STNS, filter_options)
            # Plot results
            _plot_network_relative_to_ref_station(df_final, single_ref, TARGET_STNS, batch_options,
                                                  filter_options, display_options)
        # end for
    # end for


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
