#!/bin/env python
"""
Bulk analysis script for analysing relative traveltime residuals from a pick ensemble
for the purpose of identifying time periods of GPS clock error in specific stations.
"""

import sys
import os
import datetime
import numpy as np
import pandas as pd

import pytz
import matplotlib
import matplotlib.dates
import matplotlib.pyplot as plt

# Progress bar helper to indicate that slow tasks have not stalled
import tqdm

import obspy

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()


# Priority order of trusted channels
# CHANNEL_PREF = ['BHZ_00', 'BHZ', 'BHZ_10', 'B?Z', 'S?Z', 'SHZ', '???', '?']
CHANNEL_PREF = ['BHZ_00', 'BHZ', 'BHZ_10', 'B?Z', 'S?Z', 'SHZ']
# CHANNEL_PREF = ['BHZ_00', 'BHZ', 'BHZ_10', 'B?Z']

# TODO: Bundle these filter settings into a NamedTuple or class
MIN_REF_SNR = 10
# MIN_REF_SNR = 0

CWT_CUTOFF = 15
SLOPE_CUTOFF = 3
NSIGMA_CUTOFF = 4
# CWT_CUTOFF = 0
# slope_cutoff = 0
# nsigma_cutoff = 0

# This only applies when pick quality stats are all zero.
MIN_EVENT_MAG = 5.5
# MIN_EVENT_MAG = 4.0
# MIN_EVENT_MAG = 0.0


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


def plot_target_network_rel_residuals(df, target, ref, tt_scale=50, snr_scale=(0, 60), save_file=False, file_label='', annotator=None):

    def plot_dataset(ds, net_code, ref_code):
        # Sort ds rows by SNR, so that the weakest SNR points are drawn first and the high SNR point last,
        # to make sure high SNR point are in the top rendering layer.
        ds = ds.sort_values('snr')
        times = pandas_timestamp_to_plottable_datetime(ds['originTimestamp'])
        vals = ds[yaxis].values
        qual = ds['snr'].values
        min_mag = 4.0
        mag = ds['mag'].values - min_mag
        ylabel = 'Relative TT residual (sec)'
        title = r"Network {} TT residual relative to {} (filtering: ref SNR$\geq${}, CWT$\geq${}, slope$\geq${}, $n\sigma\geq{}$)".format(
                net_code, ref_code, str(MIN_REF_SNR), str(CWT_CUTOFF), str(SLOPE_CUTOFF), str(NSIGMA_CUTOFF))
        if len(vals) > 0:
            # print(time_range[1] - time_range[0])
            plt.figure(figsize=(32, 9))
            sc = plt.scatter(times, vals, c=qual, alpha=0.5, cmap='gnuplot_r', s=np.maximum(50 * mag, 10), edgecolors=None, linewidths=0)
            time_formatter = matplotlib.dates.DateFormatter("%Y-%m-%d")
            plt.axes().xaxis.set_major_formatter(time_formatter)
            cb = plt.colorbar(sc, drawedges=False)
            cb.set_label('Signal to noise ratio', fontsize=12)
            plt.grid(color='#80808080', linestyle=':')
            plt.xlabel(xlabel, fontsize=14)
            plt.ylabel(ylabel, fontsize=14)
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            plt.xlim(time_range)
            plt.ylim((-tt_scale, tt_scale))
            plt.clim(snr_scale)
            plt.title(title, fontsize=18)
            plt.legend(['Point size = Mag - {}, Color = SNR'.format(min_mag)], fontsize=12, loc=1)
            plt.text(0.01, 0.96, "Channel selection: {}".format(CHANNEL_PREF), transform=plt.gca().transAxes, fontsize=12)
            plt.text(0.01, 0.92, "Start date: {}".format(str(time_range[0])), transform=plt.gca().transAxes, fontsize=12)
            plt.text(0.01, 0.88, "  End date: {}".format(str(time_range[1])), transform=plt.gca().transAxes, fontsize=12)
            plt.tight_layout(pad=1.05)
            if annotator is not None:
                annotator()
            if save_file:
                # subfolder = os.path.join(net_code, ref_code)
                subfolder = net_code
                os.makedirs(subfolder, exist_ok=True)
                plt_file = os.path.join(subfolder, '_'.join([net_code, ref_code]) + '_' + ylabel.replace(" ", "").replace(".*", "") + file_label + ".png")
                plt.savefig(plt_file, dpi=150)
                plt.close()
            else:
                plt.show()
    # end plotDataset

    df_times = pandas_timestamp_to_plottable_datetime(df['originTimestamp'])
    time_range = (df_times.min(), df_times.max())
    ref_code = ".".join([ref['net'][0], ref['sta'][0]])
    print("Plotting time range " + " to ".join([t.strftime("%Y-%m-%d %H:%M:%S") for t in time_range]) + " for " + ref_code)
    yaxis = 'relTtResidual'
    xlabel = 'Event Origin Timestamp'

    # Remove reference station from target set before producing composite image.
    # The reference station may not be there, but remove it if it is.
    mask_ref = df[list(ref)].isin(ref).all(axis=1)
    mask_targ = df[list(target)].isin(target).all(axis=1)
    df_agg = df[(mask_targ) & (~mask_ref)]
    plot_dataset(df_agg, ','.join(np.unique(target['net'])), ref_code)


def plot_network_relative_to_ref_station(df_plot, ref, target_stns, events=None):
    # For each event, create column for reference traveltime residual

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
    df_plot = df_plot[['#eventID', 'originTimestamp', 'mag', 'originLon', 'originLat', 'originDepthKm', 'net', 'sta', 'cha',
                       'pickTimestamp', 'phase', 'stationLon', 'stationLat', 'distance', 'snr', 'ttResidual', 'ttResidualRef',
                       'relTtResidual', 'qualityMeasureCWT', 'qualityMeasureSlope', 'nSigma']]

    # Sort data by event origin time
    df_plot = df_plot.sort_values(['#eventID', 'originTimestamp'])

    save_file = True
    if events is not None:
        plot_target_network_rel_residuals(df_plot, target_stns, ref, save_file=save_file, annotator=lambda: add_event_marker_lines(events))
    else:
        plot_target_network_rel_residuals(df_plot, target_stns, ref, save_file=save_file)


def add_event_marker_lines(events):
    time_lims = plt.xlim()
    y_lims = plt.ylim()
    for date, event in events.iterrows():
        event_time = pytz.utc.localize(datetime.datetime.strptime(date, "%Y-%m-%d"))
        if event_time < matplotlib.dates.num2date(time_lims[0]) or event_time >= matplotlib.dates.num2date(time_lims[1]):
            continue
        plt.axvline(event_time, linestyle='--', linewidth=1, color='#00800080')
        plt.text(event_time, y_lims[0] + 0.01 * (y_lims[1] - y_lims[0]), event['name'], horizontalalignment='center', verticalalignment='bottom',
                 fontsize=12, fontstyle='italic', color='#008000c0', rotation=90)


def apply_event_quality_filtering(df, ref_stn, apply_quality_to_ref=True):
    """
    Apply event quality requirements to pick events.

    :param ref_network: Network and station codes for reference network
    :type ref_network: dict of corresponding network and station codes under keys 'net' and 'sta'

    :param df: Pandas dataframe loaded from a pick event ensemble
    :type df: pandas.DataFrame
    :param ref_stn: Network and station codes for reference network (expected to be just one entry)
    :type ref_stn: dict of corresponding network and station codes under keys 'net' and 'sta'
         (expected to be just one entry)
    :param apply_quality_to_ref: Whether to apply quality criteria to the reference station events, defaults to True
    :param apply_quality_to_ref: bool, optional
    :return: Filtered dataframe of pick events
    :rtype: pandas.DataFrame
    """
    # Remove records where the SNR is too low
    mask_snr = (df['snr'] >= MIN_REF_SNR)

    # Filter to constrained quality metrics
    mask_cwt = (df['qualityMeasureCWT'] >= CWT_CUTOFF)
    mask_slope = (df['qualityMeasureSlope'] >= SLOPE_CUTOFF)
    mask_sigma = (df['nSigma'] >= NSIGMA_CUTOFF)

    # For events from ISC catalogs the quality metrics are zero, so we use event magnitude instead.
    mask_zero_quality_stats = (df[['snr', 'qualityMeasureCWT', 'qualityMeasureSlope', 'nSigma']] == 0).all(axis=1)
    mask_origin_mag = (df['mag'] >= MIN_EVENT_MAG)

    quality_mask = (mask_snr & mask_cwt & mask_slope & mask_sigma) | (mask_zero_quality_stats & mask_origin_mag)

    mask_ref = df[list(ref_stn)].isin(ref_stn).all(axis=1)
    if apply_quality_to_ref:
        # But never apply quality mask to ref stations that have all zero quality stats, as we just can't judge quality
        # and don't want to arbitrarily exclude them.
        quality_mask = (mask_zero_quality_stats & mask_ref) | quality_mask
    else:
        # Only apply quality mask stations that are not the reference station, i.e. use all ref station events
        # regardless of pick quality at the ref station. This gives more results, but possibly more noise.
        quality_mask = mask_ref | (~mask_ref & quality_mask)

    assert np.sum(quality_mask) > 100, 'Not enough points left after quality filtering'
    df_qual = df[quality_mask]
    return df_qual


def analyze_target_relative_to_ref(df_picks, ref_stn, target_stns, significant_events):
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
    :param significant_events: Pandas dataframe of historically significant seismic events, with a
        'name' and indexed by date string.
    :type significant_events: pandas.DataFrame
    """
    # Setting this to False will generate a lot more events per station chart, sometimes making it
    # easier to spot drift. But it may also add many events with significant non-zero residual.
    APPLY_QUALITY_TO_REF = False

    # Event quality filtering
    df_qual = apply_event_quality_filtering(df_picks, ref_stn, apply_quality_to_ref=APPLY_QUALITY_TO_REF)

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
    if len(df_nets) == 0:
        print("WARNING: no events left to analyze!")
        return

    # Alias for dataset at the end of all filtering, a static name that can be used from here onwards.
    ds_final = df_nets
    # print(getOverlappingDateRange(ds_final, ref_stn, target_stns))

    plot_network_relative_to_ref_station(ds_final, ref_stn, target_stns, significant_events)


def main(input_file):

    print(pd.__version__)
    print(matplotlib.__version__)

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
    print("Loading picks file {}".format(input_file))
    df_raw_picks = pd.read_csv(input_file, ' ', header=0, dtype=dtype)
    print("Number of raw picks: {}".format(len(df_raw_picks)))

    # Generate catalog of major regional events (mag 8+) for overlays
    if True:  # TODO: Make this an option
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
        significant_events.loc['2005-03-28', 'name'] = '2005 Nias–Simeulue earthquake'
        significant_events.loc['2009-09-29', 'name'] = '2009 Samoa earthquake and tsunami'
        significant_events.loc['2010-02-27', 'name'] = '2010 Chile earthquake'
        significant_events.loc['2011-03-11', 'name'] = '2011 Tohoku earthquake and tsunami'
        significant_events.loc['2012-04-11', 'name'] = '2012 Indian Ocean earthquakes'
        significant_events.loc['2013-02-06', 'name'] = '2013 Solomon Islands earthquakes'
        significant_events.loc['2013-09-24', 'name'] = '2013 Balochistan earthquakes'
        significant_events.loc['2014-04-01', 'name'] = '2014 Iquique earthquake'
        significant_events.loc['2015-09-16', 'name'] = '2015 Illapel earthquake'
        significant_events.loc['2016-08-24', 'name'] = '2016 Myanmar earthquake'

        print(significant_events)

    # Query time period for source dataset
    print("Raw picks date range: {} to {}".format(obspy.UTCDateTime(df_raw_picks['originTimestamp'].min()),
                                                  obspy.UTCDateTime(df_raw_picks['originTimestamp'].max())))

    # Remove non-BHZ channels as their picks are not considered reliable enough to use
    df_picks = df_raw_picks[df_raw_picks['cha'].isin(CHANNEL_PREF)].reset_index()
    print("Remaining picks after channel filter: {}".format(len(df_picks)))

    # Remove unused columns for readability
    df_picks = df_picks[['#eventID', 'originTimestamp', 'mag', 'originLon', 'originLat', 'originDepthKm', 'net', 'sta', 'cha', 'pickTimestamp', 'phase',
                        'stationLon', 'stationLat', 'az', 'baz', 'distance', 'ttResidual', 'snr', 'qualityMeasureCWT', 'qualityMeasureSlope', 'nSigma']]

    if True:  # TODO: Make this an option
        # Some select stations require custom date filters to remove events outside the known valid date range of a network
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

    # ## Filter to teleseismic events
    ANG_DIST = 'distance'  # angular distance (degrees) between event and station
    mask_tele = (df_picks[ANG_DIST] >= 30.0) & (df_picks[ANG_DIST] <= 90.0)
    df_picks = df_picks.loc[mask_tele]
    print("Remaining picks after filtering to teleseismic distances: {}".format(len(df_picks)))

    TARGET_NET = 'AU'
    TARGET_STN = get_network_stations(df_picks, TARGET_NET)
    TARGET_STNS = {'net': [TARGET_NET] * len(TARGET_STN), 'sta': [s for s in TARGET_STN]}
    # getNetworkDateRange(df_picks, TARGET_NET)

    # Find additional network.station codes that match AU network, and add them to target
    # TODO: Make this a command line argument.
    IRIS_AU_STATIONS_FILE = "AU_irisws-fedcatalog_20190305T012747Z.txt"
    new_nets, new_stas = determine_alternate_matching_codes(df_picks, IRIS_AU_STATIONS_FILE, TARGET_NET)
    print("Adding {} more stations from alternate international networks".format(len(new_nets)))
    TARGET_STNS['net'].extend(list(new_nets))
    TARGET_STNS['sta'].extend(list(new_stas))
    # getNetworkDateRange(df_picks, TARGET_NET)

    REF_NET = 'AU'
    REF_STN = get_network_stations(df_picks, REF_NET)
    REF_STNS = {'net': [REF_NET] * len(REF_STN), 'sta': [s for s in REF_STN]}
    # REF_STN = REF_STN[0:10]
    # custom_stns = ['ARMA', 'WB0', 'WB1', 'WB2', 'WB3', 'WB4', 'WB6', 'WB7', 'WB8', 'WB9', 'WC1', 'WC2', 'WC3', 'WC4', 'WR1', 'WR2', 'WR3', 'WR4', 'WR5', 'WR6', 'WR7', 'WR8', 'WR9',
    #                'PSA00', 'PSAA1', 'PSAA2', 'PSAA3', 'PSAB1', 'PSAB2', 'PSAB3', 'PSAC1', 'PSAC2', 'PSAC3', 'PSAD1', 'PSAD2', 'PSAD3', 'QIS', 'QLP', 'RKGY', 'STKA']
    # Hand-analyze stations with identified possible clock issues:
    # custom_stns = ['ARMA', 'HTT', 'KAKA', 'KMBL', 'MEEK', 'MOO', 'MUN', 'RKGY', 'RMQ', 'WC1', 'YNG']
    # REF_STNS = {'net': ['AU'] * len(custom_stns), 'sta': custom_stns}
    REF_STNS['net'].extend(list(new_nets))
    REF_STNS['sta'].extend(list(new_stas))

    for ref_net, ref_sta in zip(REF_STNS['net'], REF_STNS['sta']):
        print("Plotting against REF: " + ".".join([ref_net, ref_sta]))
        single_ref = {'net': [ref_net], 'sta': [ref_sta]}
        analyze_target_relative_to_ref(df_picks, single_ref, TARGET_STNS, significant_events)

    # REF_NET = 'AU'
    # REF_STN = getNetworkStations(df_picks, REF_NET)
    # # REF_STN = REF_STN[0:10]
    # REF_STNS = {'net': [REF_NET] * len(REF_STN), 'sta': [s for s in REF_STN]}

    # # TARGET_NETWORKS = ['AU', '7B', '7D', '7G', '7X']
    # TARGET_NETWORKS = ['7B', '7D', '7G', '7X']
    # for TARGET_NET in TARGET_NETWORKS:
    #     TARGET_STN = getNetworkStations(df_picks, TARGET_NET)
    #     TARGET_STNS = {'net': [TARGET_NET] * len(TARGET_STN), 'sta': [s for s in TARGET_STN]}
    #     # getNetworkDateRange(df_picks, TARGET_NET)

    #     for ref_net, ref_sta in zip(REF_STNS['net'], REF_STNS['sta']):
    #         print("Plotting against REF: " + ".".join([ref_net, ref_sta]))
    #         single_ref = {'net': [ref_net], 'sta': [ref_sta]}
    #         analyzeTargetRelativeToRef(df_picks, single_ref, TARGET_STNS, significant_events)


if __name__ == "__main__":
    # Example usage: ``python relative_tt_residuals_plotter.py "C:\data_cache\Picks\20190219\ensemble.p.txt"``
    input_file = sys.argv[1]
    main(input_file)
