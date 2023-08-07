#!/bin/env python
"""
Description:
    Uses EarthquakeTransformer to harvest picks from ASDF data-sets in parallel

References:

CreationDate:   18/05/21
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     18/05/21   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

from mpi4py import MPI
import os
import logging

from ordered_set import OrderedSet as set
import numpy as np
from datetime import datetime
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet

import click
from obspy.core import UTCDateTime, Stats
from scipy import signal
from obspy.signal.filter import bandpass, highpass, lowpass

from seismic.pick_harvester.utils import recursive_glob, split_list
from seismic.xcorqc.utils import get_stream
from seismic.xcorqc.xcorqc import taper
from seismic.misc_p import ProgressTracker
import matplotlib.pyplot as plt

from keras.models import load_model
from tensorflow.keras.optimizers import Adam
from EQTransformer.core.EqT_utils import f1, SeqSelfAttention, FeedForward, LayerNormalization
from EQTransformer.core.mseed_predictor import _picker
from collections import defaultdict

logging.basicConfig()
DAY_SECONDS = 24 * 3600

def setup_logger(name, log_file, level=logging.INFO):
    """
    Function to setup a logger; adapted from stackoverflow
    """
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    handler = logging.FileHandler(log_file, mode='w')
    handler.setFormatter(formatter)

    logger = logging.getLogger(name+log_file)
    logger.setLevel(level)
    logger.addHandler(handler)
    logger.propagate = False
    return logger
# end func

def getWorkLoad(fds:FederatedASDFDataSet, netsta_list:str,
                start_time:UTCDateTime, end_time:UTCDateTime):
    """
    @param fds:
    @param netsta_list:
    @param start_time:
    @param end_time:
    @return: list with each row containing:
             network.station, start_time and end_time for each <=24 hr period
             with data available
    """
    if(os.path.exists(netsta_list)):
        netsta_list = ' '.join(open(netsta_list).readlines()).replace('\n', ' ').strip()
    # end if

    # Gather station metadata
    netsta_list_subset = set(netsta_list.split(' ')) if netsta_list != '*' else netsta_list
    netsta_list = []

    for netsta in list(fds.unique_coordinates.keys()):
        if(netsta_list_subset != '*'):
            if netsta not in netsta_list_subset:
                continue

        netsta_list.append(netsta)
    # end for

    if(len(netsta_list) == 0): print('Warning: no stations found to process..')

    result_list = []
    for netsta in netsta_list:
        st = UTCDateTime(start_time)
        et = UTCDateTime(end_time)

        nc, sc = netsta.split('.')
        gSt, gEt = fds.get_global_time_range(nc, sc)

        if(st < gSt): st = gSt
        if(et > gEt): et = gEt

        cTime = st
        while(cTime < et):
            cStep = DAY_SECONDS
            if (cTime + cStep) > et:
                cStep = et - cTime
            # end if

            r = fds.get_stations(cTime, cTime + cStep, network=nc, station=sc)
            if(len(r)): # has data
                result_list.append([netsta, cTime, cTime + cStep])
            # end if
            cTime += cStep
        # wend
    # end for

    return result_list
# end func

def processData(ztrc, ntrc, etrc, model, picking_args,
                output_path, ofh_p, ofh_s, save_plots,
                logger=None, window_seconds=60,
                buffer_seconds=10, overlap=0.3):
    """
    @param ztrc: z-trace
    @param ntrc: n-trace
    @param etrc: e-trace
    @param model: ML-model
    @param picking_args: picking arguments
    @param output_path: output path
    @param p_ofh: file handle for outputting p-arrivals
    @param s_ofh: file handle for outputting s-arrivals
    @param save_plots flag for saving plots
    @param window_seconds: length of window in seconds
    @param buffer_seconds: data buffer length around data windows to be able to
                           exclude preprocessing artefacts
    @param overlap: window overlap
    """

    assert (ztrc.stats.sampling_rate == ntrc.stats.sampling_rate == etrc.stats.sampling_rate), \
           "Discrepant sampling rates found among channel data. Aborting.."

    REQ_SAMPLING_RATE = 100
    NSAMPLES = 6000
    sr = ztrc.stats.sampling_rate

    max_st = UTCDateTime(0)
    min_et = UTCDateTime(1e10)
    for trc in [ztrc, ntrc, etrc]:
        if(trc.stats.starttime > max_st): max_st = trc.stats.starttime
        if(trc.stats.endtime < min_et): min_et = trc.stats.endtime
    # end for

    net, sta, loc, cha = ztrc.stats.network, ztrc.stats.station, \
                         ztrc.stats.location, ztrc.stats.channel

    ml_input = np.zeros((1, NSAMPLES, 3))
    stime = max_st
    etime = min_et
    ctime = stime
    step = window_seconds - int(window_seconds * overlap)
    wcount = 0
    while(ctime < etime):

        zdata = ndata = edata = None
        sidx_start = int((wcount * step - buffer_seconds/2.) * sr)
        sidx_end   = int(sidx_start + (window_seconds + buffer_seconds/2.) * sr)

        if((sidx_start>=0) and
           (sidx_end < ztrc.data.shape[0]) and \
           (sidx_end < ntrc.data.shape[0]) and \
           (sidx_end < etrc.data.shape[0])):

            zdata = np.array(ztrc.data[sidx_start:sidx_end])
            ndata = np.array(ntrc.data[sidx_start:sidx_end])
            edata = np.array(etrc.data[sidx_start:sidx_end])

            if(type(zdata)==np.ndarray and \
               type(ndata)==np.ndarray and \
               type(edata)==np.ndarray):

                if(np.sum(zdata) != 0 and
                   np.sum(ndata) != 0 and
                   np.sum(edata) != 0): # traces cannot be all zeros

                    for i, data in enumerate([edata, ndata, zdata]):
                        data = signal.detrend(data)
                        data -= np.mean(data)
                        taper(data, int(buffer_seconds * sr / 2.))
                        data = bandpass(data, 1.0, 45, sr, corners=2, zerophase=True)

                        if(sr != REQ_SAMPLING_RATE):
                            data = signal.resample(data,
                                                   (window_seconds + buffer_seconds) * REQ_SAMPLING_RATE)
                        # end if

                        ml_input[0, :, i] = data[int(buffer_seconds/2.*REQ_SAMPLING_RATE) :
                                                 int(buffer_seconds/2.*REQ_SAMPLING_RATE)+NSAMPLES]
                    # end for

                    eprob, pprob, sprob = model(ml_input)
                    matches, pick_errors, sprob = _picker(picking_args,
                                                          np.squeeze(eprob),
                                                          np.squeeze(pprob),
                                                          np.squeeze(sprob))

                    #=============================================================
                    # matches: { detection start-time:[ 0 detection-end-time,
                    #                                   1 detection-probability,
                    #                                   2 detection-uncertainty,
                    #                                   3 P-arrival,
                    #                                   4 P-probability,
                    #                                   5 P-uncertainty,
                    #                                   6 S-arrival,
                    #                                   7 S-probability,
                    #                                   8 S-uncertainty] }
                    #=============================================================

                    def generate_plots():
                        """
                        Note that file names are reused, across ranks.
                        """
                        fig, axes = plt.subplots(3,1)

                        axes[0].plot(ml_input[0, :, 0], c='k', lw=0.3)
                        axes[0].set_xlim(0, NSAMPLES)
                        axes[0].grid()

                        for k, v in matches.items():
                            if(v[3]):
                                axes[0].axvline(v[3], c='blue', alpha=0.75)
                                axes[0].text(v[3], 0, '{}'.format(v[3]), c='b', fontsize=6, zorder=10)
                            if(v[6]):
                                axes[0].axvline(v[6], c='red', alpha=0.75)
                                axes[0].text(v[6], 0, '{}'.format(v[6]), c='r', fontsize=6, zorder=10)
                        # end for

                        axes[1].plot(np.squeeze(eprob), '--', c='g', label='Event')
                        axes[1].plot(np.squeeze(pprob), '--', c='blue', label='P')
                        axes[1].plot(np.squeeze(sprob), '--', c='red', label='S')
                        axes[1].set_xlim(0, NSAMPLES)
                        axes[1].grid()
                        axes[1].set_xlabel('Samples')
                        axes[1].set_ylabel('Probability')
                        axes[1].legend()

                        f, t, ps = signal.stft(ml_input[0, :, 2], fs=REQ_SAMPLING_RATE, nperseg=80)
                        ps = np.abs(ps)
                        ps[ps<np.std(ps)/2.] = 1
                        axes[2].pcolormesh(t, f, ps, cmap='jet', rasterized=True)
                        axes[2].set_ylabel('Frequency [Hz]')
                        axes[2].set_xlabel('Time [s]')

                        plt.suptitle('%s : %s - %s' %(ztrc.stats.station, ctime.strftime('%Y-%m-%dT%H:%M:%S'),
                                                 (ctime+window_seconds).strftime('%Y-%m-%dT%H:%M:%S')),
                                     y=1.0)

                        plt.tight_layout()
                        plt.savefig('%s/plots/%s.%d.png'%(output_path, ztrc.stats.station, wcount), dpi=300)
                        plt.close(fig)
                    # end func

                    if(len(matches)):
                        logger.info('Time window: %s - %s: found %d arrivals..' %
                                    (ctime.strftime('%Y-%m-%dT%H:%M:%S.%f'),
                                     ((ctime + window_seconds).strftime('%Y-%m-%dT%H:%M:%S.%f')),
                                     len(matches)))
                        if(save_plots): generate_plots()

                        for k, v in matches.items():
                            if(v[3]):
                                # output p-arrivals
                                ts = ctime.timestamp + v[3] / REQ_SAMPLING_RATE
                                line = '{}, {}, {}, {}, {}, {}\n'.format(net, sta, loc, cha, ts, v[4])
                                ofh_p.write(line)
                            # end if

                            if(v[6]):
                                # output s-arrivals
                                ts = ctime.timestamp + v[6] / REQ_SAMPLING_RATE
                                line = '{}, {}, {}, {}, {}, {}\n'.format(net, sta, loc, cha, ts, v[7])
                                ofh_s.write(line)
                            # end if
                        # end for
                    # end if
                # end if
            # end if
        # end if


        ctime += step
        wcount += 1
    # wend
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('asdf-source', required=True,
                type=click.Path(exists=True))
@click.argument('ml-model-path', required=True,
                type=click.Path(exists=True))
@click.argument('output-path', required=True,
                type=click.Path(exists=True))
@click.option('--picking-params', type=(float, float, float),
              default=(0.5, 0.3, 0.3),
              show_default=True,
              help="Detection threshold, P-threshold and S-threshold, specified as three "
                   "space-separated floating point values. The default is (0.5, 0.3, 0.3), "
                   "which are the recommended values for the original EQTransformer model. "
                   "The corresponding recommended values for the conservative EQTransformer "
                   "model are: (0.2, 0.03, 0.03)")
@click.option('--station-names', default='*', type=str,
              show_default=True,
              help="Either station name(s) (space-delimited) or a text file containing NET.STA entries in each line to "
                   "process; default is '*', which processes all available stations.")
@click.option('--start-time', default='1970-01-01T00:00:00',
              type=str,
              show_default=True,
              help="Date and time (in UTC format) to start from; default is year 1900.")
@click.option('--end-time', default='2100-01-01T00:00:00',
              type=str,
              show_default=True,
              help="Date and time (in UTC format) to stop at; default is year 2100.")
@click.option('--zchan', default='BHZ',
              type=str,
              show_default=True,
              help="Name of z-channel. This parameter and the three following are required to "
                   "specify channel names for the stations being processed. Simple wildcards, e.g. '*Z', are "
                   "also supported.")
@click.option('--nchan', default='BHN',
              type=str,
              show_default=True,
              help="Name of n-channel")
@click.option('--echan', default='BHE',
              type=str,
              show_default=True,
              help="Name of e-channel")
@click.option('--restart', default=False, is_flag=True,
              help='Restart job')
@click.option('--save-plots', default=False, is_flag=True,
              help='Save plots in output folder. Note that file-names are recycled and '
                   'pertain to the latest block of 24 hr data being processed.')
def process(asdf_source, ml_model_path, output_path, picking_params, station_names,
            start_time, end_time, zchan, nchan, echan, restart, save_plots):
    """
    ASDF_SOURCE: Text file containing a list of paths to ASDF files
    ML_MODEL_PATH: Path to EQT Model in H5 format
    OUTPUT_PATH: Output folder \n
    """
    # initialize MPI
    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()
    proc_workload = None

    plot_output_folder = None
    if (save_plots):
        plot_output_folder = os.path.join(output_path, 'plots')
        if (rank == 0):
            if (not os.path.exists(plot_output_folder)):
                os.mkdir(plot_output_folder)
        # end if
        comm.Barrier()
    # end if

    # load ML model
    loss_weights=[0.03, 0.40, 0.58],
    loss_types=['binary_crossentropy', 'binary_crossentropy', 'binary_crossentropy'],

    model = load_model(ml_model_path,
                       custom_objects={'SeqSelfAttention': SeqSelfAttention,
                                       'FeedForward': FeedForward,
                                       'LayerNormalization': LayerNormalization,
                                       'f1': f1})
    model.compile(loss = loss_types,
                  loss_weights = loss_weights,
                  optimizer = Adam(lr = 0.001),
                  metrics = [f1])
    # define picking parameters
    picking_args = {
        "detection_threshold": picking_params[0],
        "P_threshold": picking_params[1],
        "S_threshold": picking_params[2],
    }

    if (rank == 0):
        def outputConfigParameters():
            # output config parameters
            fn = 'pick_eqt.cfg'
            fn = os.path.join(output_path, fn)

            f = open(fn, 'w+')
            f.write('Parameter Values:\n\n')
            f.write('%25s\t\t: %s\n' % ('ASDF_SOURCE', asdf_source))
            f.write('%25s\t\t: %s\n' % ('ML_MODEL', ml_model_path))
            f.write('%25s\t\t: %s\n' % ('OUTPUT_PATH', output_path))
            f.write('%25s\t\t: %s\n' % ('RESTART_MODE', 'TRUE' if restart else 'FALSE'))
            f.close()

        # end func

        outputConfigParameters()
    # end if

    startTime = UTCDateTime(start_time)
    endTime = UTCDateTime(end_time)
    fds = FederatedASDFDataSet(asdf_source, logger=None)

    if(rank == 0):
        workload = getWorkLoad(fds, station_names, startTime, endTime)
        proc_workload = split_list(workload, npartitions=nproc)
    # end if

    # broadcast workload to all procs
    proc_workload = comm.bcast(proc_workload, root=0)

    # ==================================================
    # Open output files depending on the mode of
    # operation (fresh/restart)
    # ==================================================
    ofnp = os.path.join(output_path, 'p_arrivals.%d.txt' % (rank))
    ofns = os.path.join(output_path, 's_arrivals.%d.txt' % (rank))
    ofhp = None
    ofhs = None
    if (restart == False):
        ofhp = open(ofnp, 'w+')
        ofhs = open(ofns, 'w+')
    else:
        ofhp = open(ofnp, 'a+')
        ofhs = open(ofns, 'a+')
    # end if

    # Progress tracker
    progTracker = ProgressTracker(output_folder=output_path, restart_mode=restart)

    loggerCache = defaultdict(list)
    # main loop
    for netsta, st, et in proc_workload[rank]:
        nc, sc = netsta.split('.')
        loc_pref_dict = defaultdict(lambda: None) # ignoring location codes

        if (progTracker.increment()):
            pass
        else:
            continue
        # end if

        if(netsta in loggerCache.keys()):
            logger = loggerCache[netsta]
        else:
            # set up logger for current station
            fn = os.path.join(output_path, '%s.log' % (netsta))
            logger = setup_logger('%s' % (netsta), fn)
            loggerCache[netsta] = logger
        # end if

        # get streams
        stz, stn, ste = [], [], []
        try:
            stz = get_stream(fds, nc, sc, zchan, st, et, loc_pref_dict, logger=logger)
            if(len(stz)): stn = get_stream(fds, nc, sc, nchan, st, et, loc_pref_dict, logger=logger)
            if(len(stn)): ste = get_stream(fds, nc, sc, echan, st, et, loc_pref_dict, logger=logger)
        except Exception as e:
            logger.error('\t' + str(e))
            logger.warning('\tError encountered while fetching data. Skipping along..')
        # end try

        if(len(ste)): # we have data in all three streams
            processData(stz.traces[0], stn.traces[0], ste.traces[0], model, picking_args,
                        output_path, ofhp, ofhs, save_plots, logger=logger)
        # end if
    # end for

    ofhp.close()
    ofhs.close()
    print(('Processing complete on rank %d' % (rank)))

    del fds

    # Ensuring all processes have completed before merging results
    comm.Barrier()

    # Merge results on proc 0
    if (rank == 0):
        header = '#net, sta, loc, cha, timestamp, probability \n'
        merge_results(output_path, header)
    # end if
# end func

def merge_results(output_path, header):
    search_strings = ['p_arrivals*', 's_arrivals*']
    output_fns = ['p_combined.txt', 's_combined.txt']

    for ss, ofn in zip(search_strings, output_fns):
        files = recursive_glob(output_path, ss)
        ofn = open('%s/%s' % (output_path, ofn), 'w+')
        ofn.write(header)

        data = set()
        for i, fn in enumerate(files):
            lines = open(fn, 'r').readlines()

            for j in np.arange(len(lines)):
                data.add(lines[j])
            # end for

            os.system('rm %s' % (fn))
        # end for

        for l in data:
            ofn.write(l)
        # end for

        ofn.close()
    # end for
# end func

if (__name__ == '__main__'):
    process()
# end if
