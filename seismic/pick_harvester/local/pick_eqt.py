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
from obspy import Trace
from datetime import datetime
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet

import click
from obspy.signal.rotate import rotate_ne_rt
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees
from obspy.core import UTCDateTime, Stats
from scipy import signal
from obspy.signal.filter import bandpass, highpass, lowpass

from seismic.pick_harvester.utils import CatalogCSV, ProgressTracker, recursive_glob, split_list
from seismic.xcorqc.utils import get_stream
from seismic.xcorqc.xcorqc import taper

import matplotlib.pyplot as plt

import keras
from keras import backend as K
from keras.models import load_model
from keras.optimizers import Adam
import tensorflow as tf
from EQTransformer.core.EqT_utils import f1, SeqSelfAttention, FeedForward, LayerNormalization
from EQTransformer.core.mseed_predictor import _picker

logging.basicConfig()
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

def dropBogusTraces(st, sampling_rate_cutoff=5):
    badTraces = [tr for tr in st if tr.stats.sampling_rate < sampling_rate_cutoff]

    for tr in badTraces: st.remove(tr)
# end func

def stationsToProcess(fds, netsta_list):
    if(os.path.exists(netsta_list)):
        netsta_list = ' '.join(open(netsta_list).readlines()).replace('\n', ' ').strip()
    # end if

    # Gather station metadata
    netsta_list_subset = set(netsta_list.split(' ')) if netsta_list != '*' else netsta_list
    netsta_list_result = []

    for netsta in list(fds.unique_coordinates.keys()):
        if(netsta_list_subset != '*'):
            if netsta not in netsta_list_subset:
                continue

        netsta_list_result.append(netsta)
    # end for

    return netsta_list_result
# end func

def processData(ztrc, ntrc, etrc, model, picking_args, window_seconds=60, buffer_seconds=10, overlap=0.5):
    assert(ztrc.stats.sampling_rate == ntrc.stats.sampling_rate == etrc.stats.sampling_rate)

    sr = ztrc.stats.sampling_rate

    max_st = UTCDateTime(0)
    min_et = UTCDateTime(1e10)
    for trc in [ztrc, ntrc, etrc]:
        if(trc.stats.starttime > max_st): max_st = trc.stats.starttime
        if(trc.stats.endtime < min_et): min_et = trc.stats.endtime
    # end for

    ml_input = np.zeros((1, 6000, 3))
    stime = max_st
    etime = min_et
    ctime = stime
    step = window_seconds + buffer_seconds - (window_seconds + buffer_seconds) * overlap
    wcount = 0
    while(ctime < etime):

        zdata = ndata = edata = None
        sidx_start = int(wcount * step * sr)
        sidx_end   = sidx_start + int((window_seconds + buffer_seconds) * sr)

        if((sidx_end < ztrc.data.shape[0]) and \
           (sidx_end < ntrc.data.shape[0]) and \
           (sidx_end < etrc.data.shape[0])):

            zdata = np.array(ztrc.data[sidx_start:sidx_end])
            ndata = np.array(ntrc.data[sidx_start:sidx_end])
            edata = np.array(etrc.data[sidx_start:sidx_end])

            if(type(zdata)==np.ndarray and \
               type(ndata)==np.ndarray and \
               type(edata)==np.ndarray):

                for i, data in enumerate([edata, ndata, zdata]):
                    data = signal.detrend(data)
                    data -= np.mean(data)
                    taper(data, int(buffer_seconds * sr / 2))
                    data = bandpass(data, 1.0, 45, sr, corners=2, zerophase=True)

                    if(sr != 100):
                        data = signal.resample(data, (window_seconds + buffer_seconds) * 100)
                    # end if

                    ml_input[0, :, i] = data[int(buffer_seconds/2*100) : int(buffer_seconds/2*100)+6000]
                    #plt.plot(ml_input[:, i], lw=0.5, alpha=0.2)
                # end for

                #plt.savefig('/g/data/ha3/rakib/seismic/pst/tests/plots/%d.pdf'%(wcount))
                #plt.close()

                eprob, pprob, sprob = model(ml_input)
                matches, pick_errors, sprob = _picker(picking_args,
                                                      np.squeeze(eprob),
                                                      np.squeeze(pprob),
                                                      np.squeeze(sprob))
                print((ctime, matches))
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
@click.option('--station-names', default='*', type=str,
              help="Either station name(s) (space-delimited) or a text file containing NET.STA entries in each line to "
                   "process; default is '*', which processes all available stations.")
@click.option('--start-time', default='1970-01-01T00:00:00',
              type=str,
              help="Date and time (in UTC format) to start from; default is year 1900.")
@click.option('--end-time', default='2100-01-01T00:00:00',
              type=str,
              help="Date and time (in UTC format) to stop at; default is year 2100.")
@click.option('--zchan', default='BHZ',
              type=str,
              help="Name of z-channel. This parameter and the three following are required to "
                   "specify channel names for the stations being processed. Simple wildcards, e.g. '*Z', are "
                   "also supported.")
@click.option('--nchan', default='BHN',
              type=str,
              help="Name of n-channel")
@click.option('--echan', default='BHE',
              type=str,
              help="Name of e-channel")
@click.option('--restart', default=False, is_flag=True, help='Restart job')
@click.option('--save-quality-plots', default=False, is_flag=True, help='Save plots of quality estimates')
def process(asdf_source, ml_model_path, output_path, station_names, start_time, end_time, zchan, nchan, echan, restart, save_quality_plots):
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
        "detection_threshold": 0.8,
        "P_threshold": 0.8,
        "S_threshold": 0.6,
    }

    if (rank == 0):
        def outputConfigParameters():
            # output config parameters
            fn = 'pick_eqt.%s.cfg' % (datetime.now().strftime('%Y-%m-%d-%H-%M-%S'))
            fn = os.path.join(output_path, fn)

            f = open(fn, 'w+')
            f.write('Parameter Values:\n\n')
            f.write('%25s\t\t: %s\n' % ('ASDF_SOURCE', asdf_source))
            f.write('%25s\t\t: %s\n' % ('ML_MODEL', ml_model_path))
            f.write('%25s\t\t: %s\n' % ('OUTPUT_PATH', output_path))
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
    if (save_quality_plots):
        plot_output_folder = os.path.join(output_path, 'plots')
        if (rank == 0):
            if (not os.path.exists(plot_output_folder)):
                os.mkdir(plot_output_folder)
        # end if
        comm.Barrier()
    # end if

    fds = FederatedASDFDataSet(asdf_source, logger=None)

    if(rank == 0):
        station_list = stationsToProcess(fds, station_names)
        proc_workload = split_list(station_list, npartitions=nproc)
    # end if

    # broadcast workload to all procs
    proc_workload = comm.bcast(proc_workload, root=0)

    # ==================================================
    # Define output header and open output files
    # depending on the mode of operation (fresh/restart)
    # ==================================================
    #header = '#eventID originTimestamp mag originLon originLat originDepthKm net sta cha pickTimestamp stationLon stationLat az baz distance ttResidual snr qualityMeasureCWT domFreq qualityMeasureSlope bandIndex nSigma\n'
    header = 'net sta cha pickTimestamp stationLon stationLat\n'
    ofnp = os.path.join(output_path, 'p_arrivals.%d.txt' % (rank))
    ofns = os.path.join(output_path, 's_arrivals.%d.txt' % (rank))
    ofp = None
    ofs = None
    if (restart == False):
        ofp = open(ofnp, 'w+')
        ofs = open(ofns, 'w+')
        ofp.write(header)
        ofs.write(header)
    else:
        ofp = open(ofnp, 'a+')
        ofs = open(ofns, 'a+')
    # end if

    # Progress tracker
    progTracker = ProgressTracker(output_folder=output_path, restart_mode=restart)

    # main loop
    startTime = UTCDateTime(start_time)
    endTime = UTCDateTime(end_time)
    step = 24 * 3600 # day
    for netsta in proc_workload[rank]:
        # set up logger for current station
        fn = os.path.join(output_path, '%s.log' % (netsta))
        logger = setup_logger('%s' % (netsta), fn)

        startTime = UTCDateTime(start_time)
        endTime = UTCDateTime(end_time)

        cTime = startTime
        while cTime < endTime:
            cStep = step

            if (cTime + cStep) > endTime:
                cStep = endTime - cTime
            # end if

            if (progTracker.increment()):
                pass
            else:
                continue
            # end if

            # get streams
            nc, sc = netsta.split('.')
            stz, stn, ste = [], [], []
            stz = get_stream(fds, nc, sc, zchan, cTime, cTime + cStep, logger=logger)
            if(len(stz)): stn = get_stream(fds, nc, sc, nchan, cTime, cTime + cStep, logger=logger)
            if(len(stn)): ste = get_stream(fds, nc, sc, echan, cTime, cTime + cStep, logger=logger)

            if(len(ste)): # we have data in all three streams
                processData(stz.traces[0], stn.traces[0], ste.traces[0], model, picking_args)
            # end if

            cTime += step
        # wend
    # end for

    ofp.close()
    ofs.close()
    print(('Processing complete on rank %d' % (rank)))

    del fds

    # Ensuring all processes have completed before merging results
    comm.Barrier()

    # Merge results on proc 0
    if (rank == 0):
        merge_results(output_path)
    # end if
# end func

def merge_results(output_path):
    search_strings = ['p_arrivals*', 's_arrivals*']
    output_fns = ['p_combined.txt', 's_combined.txt']

    for ss, ofn in zip(search_strings, output_fns):
        files = recursive_glob(output_path, ss)
        ofn = open('%s/%s' % (output_path, ofn), 'w+')

        data = set()
        for i, fn in enumerate(files):
            lines = open(fn, 'r').readlines()

            if (i == 0):
                ofn.write(lines[0])
            # end if

            for j in range(1, len(lines)):
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
