#!/bin/env python

import os
import sys
import datetime
import numpy as np
import pandas as pd
# pd.set_option('display.max_rows', 200)
# pd.set_option('display.max_columns', None)
# pd.set_option('display.max_colwidth', -1)
# pd.set_option('display.width', 240)

import datetime
import pytz
import matplotlib
import matplotlib.dates
import matplotlib.pyplot as plt
# matplotlib.rcParams['figure.figsize'] = (16.0, 9.0)
# matplotlib.rcParams['figure.max_open_warning'] = 100

# Progress bar helper to indicate that slow tasks have not stalled
import tqdm

import obspy

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()


pd.__version__
matplotlib.__version__

PICKS_PATH = r"C:\data_cache\Picks\20190219\ensemble.p.txt"
# PICKS_PATH = r"C:\data_cache\Picks\20190219\ensemble_small.p.txt"

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

# Priority order of trusted channels
# channel_pref = ['BHZ_00', 'BHZ', 'BHZ_10', 'B?Z', 'S?Z', 'SHZ', '???', '?']
channel_pref = ['BHZ_00', 'BHZ', 'BHZ_10', 'B?Z']

min_ref_snr = 10
# min_ref_snr = 0

# cwt_cutoff = 10
# slope_cutoff = 2
# nsigma_cutoff = 4
cwt_cutoff = 15
slope_cutoff = 3
nsigma_cutoff = 4
# cwt_cutoff = 0
# slope_cutoff = 0
# nsigma_cutoff = 0

def getNetworkStations(df, netcode):
    return sorted(df[df['net'] == netcode]['sta'].unique().tolist())


def getNetworkMean(df, netcode):
    mean_lat = df[df['net'] == netcode]['stationLat'].mean()
    mean_lon = df[df['net'] == netcode]['stationLon'].mean()
    return (mean_lat, mean_lon)


def getNetworkDateRange(df, netcode):
    mask = (df['net'] == netcode)
    df_net = df.loc[mask]
    min_date = df_net['originTimestamp'].min()
    max_date = df_net['originTimestamp'].max()
    return (obspy.UTCDateTime(min_date), obspy.UTCDateTime(max_date))


def getStationDateRange(df, netcode, statcode):
    mask = (df['net'] == netcode)
    df_net = df.loc[mask]
    mask = (df_net['sta'] == statcode)
    df_sta = df_net.loc[mask]
    min_date = df_sta['originTimestamp'].min()
    max_date = df_sta['originTimestamp'].max()
    return (obspy.UTCDateTime(min_date), obspy.UTCDateTime(max_date))


def getOverlappingDateRange(df, ref_station, target_network):
    mask_ref = df[list(ref_station)].isin(ref_station).all(axis=1)
    mask_targ = df[list(target_network)].isin(target_network).all(axis=1)
    mask = mask_ref | mask_targ
    if not np.any(mask):
        return (None, None)
    df_nets = df.loc[mask]
    keep_events = [e for e, d in df_nets.groupby('#eventID') if np.any(d[list(ref_station)].isin(ref_station).all(axis=1)) and np.any(d[list(target_network)].isin(target_network).all(axis=1))]
    event_mask = df_nets['#eventID'].isin(keep_events)
    df_nets = df_nets[event_mask]
    return (obspy.UTCDateTime(df_nets['originTimestamp'].min()), obspy.UTCDateTime(df_nets['originTimestamp'].max()))


def display_styled_table(df):
    # Display table with blocks of same event ID highlighted
    df['lastEventID'] = df['#eventID'].shift(1)
    df['lastEventID'].iloc[0] = df['#eventID'].iloc[0]
    cols = ['#ffffff', '#e0e0ff']
    def block_highlighter(r):
        if r['lastEventID'] != r['#eventID']:
            block_highlighter.current_col = (block_highlighter.current_col + 1) % len(cols)
        return ['background-color: ' + cols[block_highlighter.current_col]]*len(r)
    block_highlighter.current_col = 0
    return df.style.apply(block_highlighter, axis=1)


def pandasTimestampToPlottableDatetime(data):
    return data.transform(datetime.datetime.utcfromtimestamp).astype('datetime64[ms]').dt.to_pydatetime()


def plotTargetNetworkRelResiduals(df, target, ref, tt_scale=50, snr_scale=(0, 60), save_file=False, file_label='', annotator=None):

    def plotDataset(ds, net_code, ref_code):
        # Sort ds rows by SNR, so that the weakest SNR points are drawn first and the high SNR point last,
        # to make sure high SNR point are in the top rendering layer.
        ds = ds.sort_values('snr')
        times = pandasTimestampToPlottableDatetime(ds['originTimestamp'])
        vals = ds[yaxis].values
        qual = ds['snr'].values
        min_mag = 4.0
        mag = ds['mag'].values - min_mag
        ylabel = 'Relative TT residual (sec)'
        title = r"Network {} TT residual relative to {} (filtering: ref SNR$\geq${}, CWT$\geq${}, slope$\geq${}, $n\sigma\geq{}$)".format(
                net_code, ref_code, str(min_ref_snr), str(cwt_cutoff), str(slope_cutoff), str(nsigma_cutoff))
        if len(vals) > 0:
            plt.figure(figsize=(32,9))
            sc = plt.scatter(times, vals, c=qual, alpha=0.5, cmap='gnuplot_r', s=np.maximum(50*mag, 10), edgecolors=None, linewidths=0)
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
            lgd = plt.legend(['Point size = Mag - {}, Color = SNR'.format(min_mag)], fontsize=12, loc=1)
            plt.text(0.01, 0.96, "Channel selection: {}".format(channel_pref), transform=plt.gca().transAxes, fontsize=12)
            plt.text(0.01, 0.92, "Start date: {}".format(str(time_range[0])), transform=plt.gca().transAxes, fontsize=12)
            plt.text(0.01, 0.88, "  End date: {}".format(str(time_range[1])), transform=plt.gca().transAxes, fontsize=12)
            if annotator is not None:
                annotator()
            if save_file:
                subfolder = os.path.join(net_code, ref_code)
                os.makedirs(subfolder, exist_ok=True)
                plt_file = os.path.join(subfolder, '_'.join([net_code, ref_code]) + '_' + ylabel.replace(" ", "").replace(".*", "") + file_label + ".png")
                plt.savefig(plt_file, dpi=150)
            else:
                plt.show()
    # end plotDataset

    df_times = pandasTimestampToPlottableDatetime(df['originTimestamp'])
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
    plotDataset(df_agg, ','.join(np.unique(target['net'])), ref_code)


def plotNetworkRelativeToRefStation(df_plot, ref, target_stns, events=None):
    # ## For each event, create column for reference traveltime residual

    # Create column for entire table first
    df_plot['ttResidualRef'] = np.nan

    pbar = tqdm.tqdm(total=len(df_plot), ascii=True)
    pbar.set_description("Broadcasting REF STN residuals")
    for eventid, grp in df_plot.groupby('#eventID'):
        pbar.update(len(grp))
        ref_mask = (grp['net'] == ref['net'][0]) & (grp['sta'] == ref['sta'][0])  # TODO: remove direct addressing of [0]
        grp_ref = grp[ref_mask]
        if grp_ref.empty:
            continue
        # Choose most favourable channel
        cha = None
        available_cha = grp_ref['cha'].values
        for c in channel_pref:
            if c in available_cha:
                cha = c
                break
        # We must find a channel
        if cha is None:
            print("WARNING: Channels {} are not amongst allowed channels {}".format(available_cha, channel_pref))
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

    if events is not None:
        plotTargetNetworkRelResiduals(df_plot, target_stns, ref, save_file=True, annotator=lambda: addEventMarkerLines(events))
    else:
        plotTargetNetworkRelResiduals(df_plot, target_stns, ref, save_file=True)


def addEventMarkerLines(events):
    time_lims = plt.xlim()
    y_lims = plt.ylim()
    for date, event in events.iterrows():
        event_time = pytz.utc.localize(datetime.datetime.strptime(date, "%Y-%m-%d"))
        if event_time < matplotlib.dates.num2date(time_lims[0]) or event_time >= matplotlib.dates.num2date(time_lims[1]):
            continue
        plt.axvline(event_time, linestyle='--', linewidth=1, color='#00800080')
        plt.text(event_time, y_lims[0] + 0.01 * (y_lims[1] - y_lims[0]), event['name'], horizontalalignment='center', verticalalignment='bottom',
                 fontsize=12, fontstyle='italic', color='#008000c0', rotation=90)


def main():

    df_raw_picks = pd.read_csv(PICKS_PATH, ' ', header=0, dtype=dtype)
    print("Number of raw picks: {}".format(len(df_raw_picks)))

    # ## Generate catalog of major regional events (mag 8+) for overlays
    if True:
        df_mag8 = df_raw_picks[df_raw_picks['mag'] >= 8.2]
        df_mag8['day'] = df_mag8['originTimestamp'].transform(datetime.datetime.utcfromtimestamp).transform(lambda x: x.strftime("%Y-%m-%d"))
        df_mag8 = df_mag8.sort_values(['day', 'originTimestamp'])

        day_mag8_count = [(day, len(df_day)) for day, df_day in df_mag8.groupby('day')]
        dates, counts = zip(*day_mag8_count)
        mag8_dict = {'date': dates, 'counts': counts}
        mag8_events_df = pd.DataFrame(mag8_dict, columns=['date', 'counts'])

        event_count_threshold = 400
        significant_events = mag8_events_df[mag8_events_df['counts'] >= event_count_threshold]
        significant_events = significant_events.set_index('date')

        significant_events.loc['2001-06-23', 'name'] = '2001 South Peru Earthquake'
        significant_events.loc['2004-12-26', 'name'] = '2004 Indian Ocean earthquake and tsunami'
        significant_events.loc['2005-03-28', 'name'] = '2005 Nias–Simeulue earthquake'
        significant_events.loc['2010-02-27', 'name'] = '2010 Chile earthquake'
        significant_events.loc['2011-03-11', 'name'] = '2011 Tohoku earthquake and tsunami'
        significant_events.loc['2012-04-11', 'name'] = '2012 Indian Ocean earthquakes'

        print(significant_events)

    # Query time period for source dataset
    print("Raw picks date range: {} to {}".format(obspy.UTCDateTime(df_raw_picks['originTimestamp'].min()), 
                                                  obspy.UTCDateTime(df_raw_picks['originTimestamp'].max())))

    # Remove non-BHZ channels as their picks are not considered reliable enough to use
    df_picks = df_raw_picks[df_raw_picks['cha'].isin(channel_pref)].reset_index()
    print("Remaining picks after channel filter: {}".format(len(df_picks)))

    # Remove unused columns for readability
    df_picks = df_picks[['#eventID', 'originTimestamp', 'mag', 'originLon', 'originLat', 'originDepthKm', 'net', 'sta', 'cha', 'pickTimestamp', 'phase', 
                        'stationLon', 'stationLat', 'az', 'baz', 'distance', 'ttResidual', 'snr', 'qualityMeasureCWT', 'qualityMeasureSlope', 'nSigma']]

    if True:
        # Some select stations require custom date filters to remove singular events outside the date range of the rest of the network
        DATE_FILTER = (
            ('7D', 'DG41', pd.Timestamp(datetime.datetime(2012, 7, 1))), 
            ('7D', 'DF43', pd.Timestamp(datetime.datetime(2012, 7, 1)))
        )
        before = len(df_picks)
        for net, sta, min_date in DATE_FILTER:
            date_mask = (df_picks['net'] == net) & (df_picks['sta'] == sta) & (df_picks['originTimestamp'] < min_date.timestamp())
            df_picks = df_picks[~date_mask]
        after = len(df_picks)
        print('Removed {} events due to timestamps'.format(before - after))

    # ## Filter to teleseismic events
    ANG_DIST = 'distance'  # angular distance (degrees) between event and station
    mask_tele = (df_picks[ANG_DIST] >= 30.0) & (df_picks[ANG_DIST] <= 90.0)
    df_picks = df_picks.loc[mask_tele]
    print("Remaining picks after filtering to teleseismic distances: {}".format(len(df_picks)))

    TARGET_NET = 'AU'
    STN_LIST = getNetworkStations(df_picks, TARGET_NET)
    TARGET_STNS = {'net': [TARGET_NET] * len(STN_LIST), 'sta': [s for s in STN_LIST]}
    # getNetworkDateRange(df_picks, TARGET_NET)

    REF_NET = 'AU'
    REF_STN = 'HTT'  # QIS, CTA, QLP, TOO, WB2, WR0, WR2, HTT, ARMA, CMSA
    REF = {'net': [REF_NET], 'sta': [REF_STN]}
    # getStationDateRange(df_picks, REF_NET, REF_STN)

    # ## Remove reference station records where the SNR is too low
    mask_ref = df_picks[list(REF)].isin(REF).all(axis=1)
    mask_ref_snr = ~mask_ref | (mask_ref & (df_picks['snr'] >= min_ref_snr))
    df_good_ref_snr = df_picks.loc[mask_ref_snr]
    print("Remaining picks after filtering to SNR >= {}: {}".format(min_ref_snr, len(df_good_ref_snr)))

    # ## Filter to constrained quality metrics
    mask_cwt = (df_good_ref_snr['qualityMeasureCWT'] >= cwt_cutoff)
    mask_slope = (df_good_ref_snr['qualityMeasureSlope'] >= slope_cutoff)
    mask_sigma = (df_good_ref_snr['nSigma'] >= nsigma_cutoff)
    # Make sure we DON'T filter out the reference station, which may have zero quality values
    mask_ref = df_good_ref_snr[list(REF)].isin(REF).all(axis=1)
    quality_mask = (mask_cwt & mask_slope & mask_sigma) | mask_ref
    assert np.sum(quality_mask) > 100, 'Not enough points left after quality filtering'
    df_qual = df_good_ref_snr[quality_mask]
    print("Remaining picks after filtering to minimum quality metrics: {}".format(len(df_qual)))

    # ## Filter to desired ref and target networks
    mask_ref = df_qual[list(REF)].isin(REF).all(axis=1)
    mask_targ = df_qual[list(TARGET_STNS)].isin(TARGET_STNS).all(axis=1)
    mask = mask_ref | mask_targ
    np.any(mask)
    df_nets = df_qual.loc[mask]
    # len(df_nets)
    # Filter out events in which REF and TARGET stations are not both present
    print("Narrowing dataframe to events common to REF and TARGET networks...")
    keep_events = [e for e, d in df_nets.groupby('#eventID') if np.any(d[list(REF)].isin(REF).all(axis=1)) and np.any(d[list(TARGET_STNS)].isin(TARGET_STNS).all(axis=1))]
    # len(keep_events)
    event_mask = df_nets['#eventID'].isin(keep_events)
    df_nets = df_nets[event_mask]
    print("Remaining picks after narrowing to common events: {}".format(len(df_nets)))
    if len(df_nets) == 0:
        print("WARNING: no events left to analyze!")
        return

    # Alias for dataset at the end of all filtering, a static name that can be used from here onwards.
    ds_final = df_nets
    # print(getOverlappingDateRange(ds_final, REF, TARGET_STNS))

    # plotNetworkRelativeToRefStation(ds_final, REF, TARGET_STNS)
    plotNetworkRelativeToRefStation(ds_final, REF, TARGET_STNS, significant_events)


if __name__ == "__main__":
    main()
