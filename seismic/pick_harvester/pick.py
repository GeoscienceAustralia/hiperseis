#!/bin/env python
"""
Description:
    Harvests picks from ASDF data-sets in parallel
References:

CreationDate:   13/09/18
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     13/09/18   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

from mpi4py import MPI
import glob, os, sys

import numpy as np
from obspy import Stream, Trace, UTCDateTime
from datetime import datetime
from ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet

import click
import traceback
from obspy import UTCDateTime
from obspy.taup import TauPyModel
from obspy.signal.rotate import rotate_ne_rt
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees
from PhasePApy.phasepapy.phasepicker import fbpicker
from PhasePApy.phasepapy.phasepicker import ktpicker
from PhasePApy.phasepapy.phasepicker import aicdpicker

from utils import EventParser, Catalog, CatalogCSV, ProgressTracker
import psutil
import gc

from quality import compute_quality_measures

def extract_p(taupy_model, pickerlist, event, station_longitude, station_latitude,
              st, win_start=-50, win_end=50, resample_hz=20,
              bp_freqmins = [0.5, 2., 5.],
              bp_freqmaxs = [5., 10., 10.],
              margin=None,
              max_amplitude=1e8,
              plot_output_folder=None):

    po = event.preferred_origin
    if(not po): return None

    atimes = []
    try:
        atimes = taupy_model.get_travel_times_geo(po.depthkm, po.lat,
                                                  po.lon, station_latitude,
                                                  station_longitude,
                                                  phase_list=('P',))
    except:
        return None
    # end try

    if (len(atimes) == 0): return None
    tat = atimes[0].time # theoretical arrival time

    try:
        buffer_start = -5
        buffer_end = 5

        snrst = st.slice(po.utctime + tat + win_start + buffer_start, po.utctime + tat + win_end + buffer_end)
        snrst = snrst.copy()
        snrst.resample(resample_hz)
        snrst.detrend('linear')

        st = snrst.slice(po.utctime + tat + win_start, po.utctime + tat + win_end)
    except:
        return None
    # end try

    if(len(st) == 0 or len(snrst) == 0): return None

    tr = st[0]
    snrtr = snrst[0]
    if(type(tr.data) == np.ndarray):
        if(np.max(tr.data) > max_amplitude): return None

        pickslist = []
        snrlist = []
        residuallist = []
        bandindex = -1
        pickerindex = -1

        foundpicks = False
        for i in range(len(bp_freqmins)):
            trc = tr.copy()
            trc.filter('bandpass', freqmin=bp_freqmins[i],
                       freqmax=bp_freqmaxs[i], corners=4,
                       zerophase=True)

            for ipicker, picker in enumerate(pickerlist):
                try:
                    scnl, picks, polarity, snr, uncert = picker.picks(trc)

                    for ipick, pick in enumerate(picks):
                        actualArrival = pick - po.utctime
                        residual = actualArrival - tat

                        if ((margin and np.fabs(residual) < margin) or (margin == None)):
                            pickslist.append(pick)

                            plotinfo = None
                            if (plot_output_folder):
                                plotinfo = {'eventid': event.public_id,
                                            'origintime': po.utctime,
                                            'mag': event.preferred_magnitude.magnitude_value,
                                            'net': trc.stats.network,
                                            'sta': trc.stats.station,
                                            'phase': 'p',
                                            'ppsnr': snr[ipick],
                                            'pickid': ipick,
                                            'outputfolder': plot_output_folder}
                            # end if

                            wab = snrtr.slice(pick - 3, pick + 3)
                            scales = np.logspace(0.15, 1.5, 30)
                            cwtsnr, dom_freq, slope_ratio = compute_quality_measures(wab, scales, plotinfo)
                            snrlist.append([snr[ipick], cwtsnr, dom_freq, slope_ratio])

                            residuallist.append(residual)
                            bandindex = i
                            pickerindex = ipicker

                            foundpicks = True
                        # end if
                    # end for
                except:
                    continue
                # end try
                if (foundpicks): break
            # end for
            if (foundpicks): break
        # end for

        if (len(pickslist)):
            return pickslist, residuallist, \
                   np.array(snrlist), bandindex, pickerindex
        # end if
    # end if

    return None
# end func

def extract_s(taupy_model, pickerlist, event, station_longitude, station_latitude,
              stn, ste, ba, win_start=-50, win_end=50, resample_hz=20,
              bp_freqmins = [0.01, 0.01, 0.5],
              bp_freqmaxs = [   1,   2., 5.],
              margin=None,
              max_amplitude=1e8,
              plot_output_folder=None):

    po = event.preferred_origin
    if(not po): return None

    atimes = []
    try:
        atimes = taupy_model.get_travel_times_geo(po.depthkm, po.lat,
                                                  po.lon, station_latitude,
                                                  station_longitude,
                                                  phase_list=('S',))
    except:
        return None
    # end try

    if (len(atimes) == 0): return None
    tat = atimes[0].time # theoretical arrival time

    tr = None
    snrtr = None
    try:
        buffer_start = -5
        buffer_end = 5
        stn = stn.slice(po.utctime + tat + win_start + buffer_start, po.utctime + tat + win_end + buffer_end)
        stn = stn.copy()
        stn.resample(resample_hz)
        stn.detrend('linear')

        if(ste):
            ste = ste.slice(po.utctime + tat + win_start + buffer_start, po.utctime + tat + win_end + buffer_end)
            ste = ste.copy()
            ste.resample(resample_hz)
            ste.detrend('linear')
        # end if

        if(ste):
            if(type(stn[0].data) == np.ndarray and type(ste[0].data) == np.ndarray):
                rc, tc = rotate_ne_rt(stn[0].data, ste[0].data, ba)
                snrtr = Trace(data=tc, header=stn[0].stats)
                tr = snrtr.slice(po.utctime + tat + win_start, po.utctime + tat + win_end)
                #tr = Trace(data=np.sqrt(np.power(rc,2) + np.power(tc,2)), header=stn[0].stats)
            # end if
        else:
            if(type(stn[0].data) == np.ndarray):
                snrtr = stn[0]
                tr = snrtr.slice(po.utctime + tat + win_start, po.utctime + tat + win_end)
            # end if
        # end if
    except Exception as e:
        return None
    # end try

    if(tr):
        if(np.max(tr.data) > max_amplitude): return None

        pickslist = []
        snrlist = []
        residuallist = []
        bandindex = -1
        pickerindex = -1

        foundpicks = False
        for i in range(len(bp_freqmins)):
            trc = tr.copy()
            trc.filter('bandpass', freqmin=bp_freqmins[i],
                       freqmax=bp_freqmaxs[i], corners=4,
                       zerophase=True)

            for ipicker, picker in enumerate(pickerlist):
                try:
                    scnl, picks, polarity, snr, uncert = picker.picks(trc)

                    for ipick, pick in enumerate(picks):
                        actualArrival = pick - po.utctime
                        residual = actualArrival - tat

                        if ((margin and np.fabs(residual) < margin) or (margin == None)):
                            pickslist.append(pick)

                            plotinfo = None
                            if (plot_output_folder):
                                plotinfo = {'eventid': event.public_id,
                                            'origintime': po.utctime,
                                            'mag': event.preferred_magnitude.magnitude_value,
                                            'net': trc.stats.network,
                                            'sta': trc.stats.station,
                                            'phase': 's',
                                            'ppsnr': snr[ipick],
                                            'pickid': ipick,
                                            'outputfolder': plot_output_folder}
                            # end if

                            wab = snrtr.slice(pick - 3, pick + 3)
                            scales = np.logspace(0.5, 4, 30)
                            cwtsnr, dom_freq, slope_ratio = compute_quality_measures(wab, scales, plotinfo)
                            snrlist.append([snr[ipick], cwtsnr, dom_freq, slope_ratio])

                            residuallist.append(residual)
                            bandindex = i
                            pickerindex = ipicker

                            foundpicks = True
                        # end if
                    # end for
                except:
                    continue
                # end try
                if (foundpicks): break
            # end for
            if (foundpicks): break
        # end for

        if (len(pickslist)):
            return pickslist, residuallist, \
                   np.array(snrlist), bandindex, pickerindex
        # end if
    # end if

    return None
# end func

def getWorkloadEstimate(fds, originTimestamps):
    totalTraceCount = 0
    for nc, sc, start_time, end_time in fds.local_net_sta_list():

        day = 24 * 3600
        curr = start_time
        step = day
        while (curr < end_time):
            if (curr + step > end_time):
                step = end_time - curr
            # end if

            eventIndices = (np.where((originTimestamps >= curr.timestamp) & \
                                     (originTimestamps <= (curr + day).timestamp)))[0]

            if(eventIndices.shape[0]>0): totalTraceCount += 1
            curr += step
        # wend
    # end for
    return totalTraceCount
# end func

def dropBogusTraces(st, sampling_rate_cutoff=5):
    badTraces = [tr for tr in st if tr.stats.sampling_rate < sampling_rate_cutoff]

    for tr in badTraces: st.remove(tr)
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('asdf-source', required=True,
                type=click.Path(exists=True))
@click.argument('event-folder', required=True,
                type=click.Path(exists=True))
@click.argument('output-path', required=True,
                type=click.Path(exists=True))
@click.option('--min-magnitude', default=4.0, help='Minimum magnitude of event')
@click.option('--restart', default=False, is_flag=True, help='Restart job')
@click.option('--save-quality-plots', default=False, is_flag=True, help='Save plots of quality estimates')
def process(asdf_source, event_folder, output_path, min_magnitude, restart, save_quality_plots):
    """
    ASDF_SOURCE: Text file containing a list of paths to ASDF files
    EVENT_FOLDER: Path to folder containing event files\n
    OUTPUT_PATH: Output folder \n
    """

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()
    proc_workload = None

    if(rank == 0):
        def outputConfigParameters():
            # output config parameters
            fn = 'pick.%s.cfg' % (datetime.now().strftime('%Y-%m-%d-%H-%M-%S'))
            fn = os.path.join(output_path, fn)

            f = open(fn, 'w+')
            f.write('Parameter Values:\n\n')
            f.write('%25s\t\t: %s\n' % ('ASDF_SOURCE', asdf_source))
            f.write('%25s\t\t: %s\n' % ('EVENT_FOLDER', event_folder))
            f.write('%25s\t\t: %s\n' % ('OUTPUT_PATH', output_path))
            f.write('%25s\t\t: %s\n' % ('MIN_MAGNITUDE', min_magnitude))
            f.write('%25s\t\t: %s\n' % ('RESTART_MODE', 'TRUE' if restart else 'FALSE'))
            f.write('%25s\t\t: %s\n' % ('SAVE_PLOTS', 'TRUE' if save_quality_plots else 'FALSE'))
            f.close()
        # end func

        outputConfigParameters()
    # end if

    # ==================================================
    # Create output-folder for snr-plots
    # ==================================================
    plot_output_folder = None
    if(save_quality_plots):
        plot_output_folder = os.path.join(output_path, 'plots')
        if(rank == 0):
            if(not os.path.exists(plot_output_folder)):
                os.mkdir(plot_output_folder)
        # end if
        comm.Barrier()
    # end if

    # ==================================================
    # Read catalogue and retrieve origin times
    # ==================================================
    cat = CatalogCSV(event_folder)
    events = cat.get_events()
    originTimestamps = cat.get_preferred_origin_timestamps()

    # ==================================================
    # Create lists of pickers for both p- and s-arrivals
    # ==================================================
    sigmalist = np.arange(8,3,-1)
    pickerlist_p = []
    pickerlist_s = []
    for sigma in sigmalist:
        picker_p = aicdpicker.AICDPicker(t_ma=5, nsigma=sigma, t_up=1, nr_len=5,
                                       nr_coeff=2, pol_len=10, pol_coeff=10, uncert_coeff=3)
        picker_s = aicdpicker.AICDPicker(t_ma=15, nsigma=sigma, t_up=1, nr_len=5,
                                       nr_coeff=2, pol_len=10, pol_coeff=10, uncert_coeff=3)

        pickerlist_p.append(picker_p)
        pickerlist_s.append(picker_s)
    # end for

    # ==================================================
    # Define theoretical model
    # Instantiate data-access object
    # Retrieve estimated workload
    # ==================================================
    taupyModel = TauPyModel(model='iasp91')
    fds = FederatedASDFDataSet(asdf_source, use_json_db=False, logger=None)
    workload = getWorkloadEstimate(fds, originTimestamps)

    # ==================================================
    # Define output header and open output files
    # depending on the mode of operation (fresh/restart)
    # ==================================================
    header = '#eventID originTimestamp mag originLon originLat originDepthKm net sta cha pickTimestamp stationLon stationLat az baz distance ttResidual snr qualityMeasureCWT domFreq qualityMeasureSlope bandIndex nSigma\n'
    ofnp = os.path.join(output_path, 'p_arrivals.%d.txt' % (rank))
    ofns = os.path.join(output_path, 's_arrivals.%d.txt' % (rank))
    ofp = None
    ofs = None
    if(restart == False):
        ofp = open(ofnp, 'w+')
        ofs = open(ofns, 'w+')
        ofp.write(header)
        ofs.write(header)
    else:
        ofp = open(ofnp, 'a+')
        ofs = open(ofns, 'a+')
    # end if

    progTracker = ProgressTracker(output_folder=output_path, restart_mode=restart)
    totalTraceCount = 0
    for nc, sc, start_time, end_time in fds.local_net_sta_list():
        day = 24 * 3600
        dayCount = 0
        curr = start_time
        traceCountP = 0
        pickCountP = 0
        traceCountS = 0
        pickCountS = 0
        sw_start = datetime.now()
        step = day
        while (curr < end_time):
            if (curr + step > end_time):
                step = end_time - curr
            # end if

            eventIndices = (np.where((originTimestamps >= curr.timestamp) & \
                                     (originTimestamps <= (curr + day).timestamp)))[0]

            if(eventIndices.shape[0]>0):
                totalTraceCount += 1
                stations = fds.get_stations(curr, curr + day, network=nc, station=sc)
                stations_zch = [s for s in stations if 'Z' in s[3]]  # only Z channels
                stations_nch = [s for s in stations if 'N' in s[3] or '1' in s[3]]  # only N channels
                stations_ech = [s for s in stations if 'E' in s[3] or '2' in s[3]]  # only E channels

                for codes in stations_zch:
                    if(progTracker.increment()): pass
                    else: continue

                    st = fds.get_waveforms(codes[0], codes[1], codes[2], codes[3],
                                           curr,
                                           curr + step,
                                           automerge=True,
                                           trace_count_threshold=200)

                    if (len(st) == 0): continue
                    dropBogusTraces(st)

                    slon, slat = codes[4], codes[5]
                    for ei in eventIndices:
                        event = events[ei]
                        po = event.preferred_origin
                        da = gps2dist_azimuth(po.lat, po.lon, slat, slon)
                        mag = None
                        if(event.preferred_magnitude): mag = event.preferred_magnitude.magnitude_value
                        elif(len(po.magnitude_list)): mag = po.magnitude_list[0].magnitude_value
                        if(mag == None): mag = np.NaN

                        if(np.isnan(mag) or mag < min_magnitude): continue

                        result = extract_p(taupyModel, pickerlist_p, event, slon, slat, st,
                                           plot_output_folder = plot_output_folder)
                        if(result):
                            picklist, residuallist, snrlist, bandindex, pickerindex = result

                            arcdistance = kilometers2degrees(da[0]/1e3)
                            for ip, pick in enumerate(picklist):
                                line = '%s %f %f %f %f %f ' \
                                       '%s %s %s %f %f %f ' \
                                       '%f %f %f ' \
                                       '%f %f %f %f %f '\
                                       '%d %d\n' % (event.public_id, po.utctime.timestamp, mag, po.lon, po.lat, po.depthkm,
                                                    codes[0], codes[1], codes[3], pick.timestamp, slon, slat,
                                                    da[1], da[2], arcdistance,
                                                    residuallist[ip], snrlist[ip, 0], snrlist[ip, 1], snrlist[ip, 2], snrlist[ip, 3],
                                                    bandindex, sigmalist[pickerindex])
                                ofp.write(line)
                            # end for
                            ofp.flush()
                            pickCountP += 1
                        # end if

                        if (len(stations_nch) == 0 and len(stations_ech) == 0):
                            result = extract_s(taupyModel, pickerlist_s, event, slon, slat, st, None, da[2],
                                               plot_output_folder=plot_output_folder)
                            if (result):
                                picklist, residuallist, snrlist, bandindex, pickerindex = result

                                arcdistance = kilometers2degrees(da[0] / 1e3)
                                for ip, pick in enumerate(picklist):
                                    line = '%s %f %f %f %f %f ' \
                                           '%s %s %s %f %f %f ' \
                                           '%f %f %f ' \
                                           '%f %f %f %f %f ' \
                                           '%d %d\n' % (event.public_id, po.utctime.timestamp, mag, po.lon, po.lat, po.depthkm,
                                                        codes[0], codes[1], codes[3], pick.timestamp, slon, slat,
                                                        da[1], da[2], arcdistance,
                                                        residuallist[ip], snrlist[ip, 0], snrlist[ip, 1], snrlist[ip, 2], snrlist[ip, 3],
                                                        bandindex, sigmalist[pickerindex])
                                    ofs.write(line)
                                # end for
                                ofs.flush()
                                pickCountS += 1
                            # end if
                        # end if
                    # end for

                    traceCountP += len(st)
                # end for

                if(len(stations_nch)>0 and len(stations_nch) == len(stations_ech)):
                    for codesn, codese in zip(stations_nch, stations_ech):
                        if (progTracker.increment()): pass
                        else: continue

                        stn = fds.get_waveforms(codesn[0], codesn[1], codesn[2], codesn[3],
                                               curr,
                                               curr + step,
                                               automerge=True,
                                               trace_count_threshold=200)
                        ste = fds.get_waveforms(codese[0], codese[1], codese[2], codese[3],
                                               curr,
                                               curr + step,
                                               automerge=True,
                                               trace_count_threshold=200)

                        dropBogusTraces(stn)
                        dropBogusTraces(ste)

                        if (len(stn) == 0): continue
                        if (len(ste) == 0): continue

                        slon, slat = codesn[4], codesn[5]

                        for ei in eventIndices:
                            event = events[ei]
                            po = event.preferred_origin
                            da = gps2dist_azimuth(po.lat, po.lon, slat, slon)

                            mag = None
                            if (event.preferred_magnitude):
                                mag = event.preferred_magnitude.magnitude_value
                            elif (len(po.magnitude_list)):
                                mag = po.magnitude_list[0].magnitude_value
                            if (mag == None): mag = np.NaN

                            if (np.isnan(mag) or mag < min_magnitude): continue

                            result = extract_s(taupyModel, pickerlist_s, event, slon, slat, stn, ste, da[2],
                                               plot_output_folder=plot_output_folder)
                            if (result):
                                picklist, residuallist, snrlist, bandindex, pickerindex = result

                                arcdistance = kilometers2degrees(da[0] / 1e3)
                                for ip, pick in enumerate(picklist):
                                    line = '%s %f %f %f %f %f ' \
                                           '%s %s %s %f %f %f ' \
                                           '%f %f %f ' \
                                           '%f %f %f %f %f ' \
                                           '%d %d\n' % (event.public_id, po.utctime.timestamp, mag, po.lon, po.lat, po.depthkm,
                                                        codesn[0], codesn[1], '00T', pick.timestamp, slon, slat,
                                                        da[1], da[2], arcdistance,
                                                        residuallist[ip], snrlist[ip, 0], snrlist[ip, 1], snrlist[ip, 2], snrlist[ip, 3],
                                                        bandindex, sigmalist[pickerindex])
                                    ofs.write(line)
                                # end for
                                ofs.flush()
                                pickCountS += 1
                            # end if
                        # end for

                        traceCountS += (len(stn) + len(ste))
                    # end for
                # end if
            # end if
            curr += step
            dayCount += 1
        # wend
        sw_stop = datetime.now()
        totalTime = (sw_stop - sw_start).total_seconds()

        gc.collect()
        print '(Rank %d: %5.2f%%, %d/%d) Processed %d traces and found %d p-arrivals and %d s-arrivals for ' \
              'network %s station %s in %f s. Memory usage: %5.2f MB.' % \
              (rank, (float(totalTraceCount) / float(workload) * 100) if workload > 0 else 100, totalTraceCount, workload,
               traceCountP + traceCountS, pickCountP, pickCountS, nc, sc, totalTime,
               round(psutil.Process().memory_info().rss / 1024. / 1024., 2))
    # end for
    ofp.close()
    ofs.close()

    print 'Processing complete on rank %d'%(rank)

    del fds
# end func

if (__name__ == '__main__'):
    process()
# end if
