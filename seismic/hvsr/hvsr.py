from seismic.hvsr.batch import create_HVSR
from mpi4py import MPI
import os
from obspy.core import UTCDateTime
import numpy as np
import matplotlib
from collections import defaultdict
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import click
from seismic.hvsr.utils import waveform_iterator3c
from seismic.xcorqc.utils import SpooledMatrix
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
import re
from copy import deepcopy
import uuid
from matplotlib.backends.backend_pdf import PdfPages
from seismic.receiver_fn.rf_plot_utils import pdf_merge
from tqdm import tqdm

def generate_master_curve(station:str, output_path:str,
                          sm:SpooledMatrix, hvsr_freq:np.ndarray,
                          master_method:str, cutoff_value:float,
                          clip_freq:bool, lowest_freq:float,
                          highest_freq:float, window_length:float):
    """
    :param master_curve_method: string
        How to determine the master curve. Available are 'mean', 'geometric
        average', 'full gaussian' , and 'median'
    :param cutoff_value: float
        If given, than this value determines which part of the bottom and top
        frequencies are discarded for the calculation of the master HVSR curve.
        e.g. a value of 0.1 will throw away the bottom and top 10% for each
        frequency.
    """
    def find_nearest_idx(array, value): return (np.abs(array - value)).argmin()

    if((sm.nrows == 0) or (hvsr_freq is None)): return None

    hvsr_matrix = sm.get_matrix()
    nwindows = len(hvsr_matrix)
    days_processed = nwindows * window_length / 86400.

    # Copy once to be able to calculate standard deviations.
    original_matrix = deepcopy(hvsr_matrix)

    # Sort it for quantile operations.
    hvsr_matrix.sort(axis=0)

    # Only senseful for mean calculations. Omitted for the median.
    if cutoff_value != 0.0 and master_method != 'median':
        hvsr_matrix = hvsr_matrix[int(nwindows * cutoff_value):
                                  int(np.ceil(nwindows * (1 - cutoff_value))), :]
    # end if

    # Mean.
    if master_method == 'mean':
        master_curve = hvsr_matrix.mean(axis=0)
    # Geometric average.
    elif master_method == 'geometric average':
        master_curve = hvsr_matrix.prod(axis=0) ** (1.0 / nwindows)
    # Median.
    elif master_method == 'median':
        master_curve = np.empty(len(hvsr_matrix[0, :]))
        error = np.empty((len(master_curve), 2))
        for _i in range(len(master_curve)):
            cur_row = hvsr_matrix[:, _i]
            master_curve[_i] = np.quantile(cur_row, 0.5)
            error[_i, 0] = np.quantile(cur_row, 0.25)
            error[_i, 1] = np.quantile(cur_row, 0.75)
    # end if

    std = (np.log1p(hvsr_matrix[:][:]) - np.log1p(master_curve))
    errormag = np.zeros(nwindows)
    for i in range(nwindows):
        errormag[i] = np.dot(std[i, :], std[i, :].T)
    error = np.dot(std.T, std)
    error /= float(nwindows - 1)

    if clip_freq:
        lclip = find_nearest_idx(hvsr_freq, lowest_freq)
        uclip = find_nearest_idx(hvsr_freq, highest_freq)
        hvsr_matrix = hvsr_matrix[:, lclip:uclip]
        master_curve = master_curve[lclip:uclip]
        hvsr_freq = hvsr_freq[lclip:uclip]
        error = error[lclip:uclip, lclip:uclip]
    # end if

    #========================================================
    # find peaks
    #========================================================
    peak_indices = []
    for iw in np.arange(len(hvsr_matrix)):
        peak_indices.append(np.argmax(hvsr_matrix[iw, :]))
    # end for
    peak_indices = np.array(peak_indices, dtype='i4')
    f0_mean = f0_std = None
    if(len(peak_indices)):
        valid_indices = np.nonzero(hvsr_freq[peak_indices])[0]
        if(len(valid_indices)):
            #print(hvsr_freq[peak_indices[valid_indices]], np.log(hvsr_freq[peak_indices[valid_indices]]))
            f0_mean = np.exp(np.mean(np.log(hvsr_freq[peak_indices[valid_indices]])))
            f0_std = np.std(hvsr_freq[peak_indices[valid_indices]])
            f0_std = f0_mean * np.sqrt(np.exp(f0_std*f0_std) - 1)
        # end if
    # end if

    #print("Master curve shape: " + str(master_curve.shape))
    #print("Frequencies shape: " + str(hvsr_freq.shape))
    #print("Error shape: " + str(error.shape))

    diagerr = np.sqrt(np.diag(error))
    lerr = np.exp(np.log(master_curve) - diagerr)
    uerr = np.exp(np.log(master_curve) + diagerr)

    output_file = os.path.join(output_path, '{}.pdf'.format(station))
    with PdfPages(output_file) as pdf:
        fig, ax = plt.subplots()
        fig.set_size_inches(7, 7)

        # Add amount of data processed
        ax.plot([], [], ' ', label='#win: {} (~{:0.2f} days)'.format(nwindows, days_processed))

        # Plot master curves and +/- 1 sigma
        ax.plot(hvsr_freq, master_curve, 'k-', lw=1, label='Master-curve')
        ax.fill_between(hvsr_freq, lerr, uerr, color='lightgray', label='Master-curve ± 1 STD')

        if(f0_mean):
            ymin = np.min(hvsr_matrix)
            ymax = np.max(hvsr_matrix)
            ax.vlines(f0_mean, ymin, ymax,
                      colors='r', linestyles='solid', label='Mean f0: {:0.2f} Hz'.format(f0_mean))
            ax.fill_betweenx([ymin, ymax],
                             f0_mean-f0_std, f0_mean+f0_std, color='salmon',
                             label='f0 ± 1 STD', alpha=0.3)
        # end if
        ax.legend()

        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.grid(which='major', linestyle='-', linewidth=0.5, color='k')
        ax.grid(which='minor', linestyle=':', linewidth=0.5, color='grey')
        plt.suptitle(station)
        pdf.savefig(dpi=300)
        plt.close()
    # end with
    return output_file, hvsr_freq, master_curve, lerr, uerr
# end func

def get_stations_to_process(fds:FederatedASDFDataSet, network:str, station_list:list):
    netsta_list= fds.unique_coordinates.keys()
    result = []
    net_dict = defaultdict(list)

    for netsta in netsta_list: net_dict[netsta.split('.')[0]].append(netsta.split('.')[1])

    if(type(station_list) == list):
        for station in station_list:
            if(station not in net_dict[network]):
                raise ValueError('Station {} not found. Aborting..'.format(station))
            else:
                result.append('{}.{}'.format(network, station))
            # end if
        # end for
    elif(station_list == '*'):
        for station in net_dict[network]: result.append('{}.{}'.format(network, station))
    else:
        raise ValueError('Invalid station-list: {}. Aborting..'.format(station_list))
    # end if

    return result
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('asdf-source',
                type=click.Path(exists=True))
@click.argument('network', type=str)
@click.argument('spec-method',
                required=True,
                type=click.Choice(['single-taper', 'multitaper', 'st', 'cwt2']))
@click.argument('output-path', required=True,
                type=click.Path(exists=True))
@click.option('--win-length', default=200, help="Window length in seconds")
@click.option('--trigger-method', default='zdetect', type=click.Choice(['zdetect', 'stalta']),
              help="Triggering method to use")
@click.option('--trigger-wlen', default=0.5, type=float, help="Triggering window length in seconds if method='zdetect', otherwise"
                                                              " this is window size for short time average in 'stalta'")
@click.option('--trigger-wlen-long', default=30, type=float, help="Window length in seconds for long time average if method='stalta'; "
              "this parameter has no effect if method='zdetect'")
@click.option('--trigger-threshold', default=0.95, type=float, help="Threshold, as a percentile, for the characteristic function to find quiet areas")
@click.option('--trigger-lowpass-value', default=None, type=float, help="Lowpass filter value (Hz) to use for triggering only and not for data "
                                                            "data processing")
@click.option('--trigger-highpass-value', default=None, type=float, help="Highpass filter value (Hz) to use for triggering only and not for data "
                                                             "data processing")
@click.option('--resample-rate', default=None, type=int, help="Traces are resampled to this sampling rate as a pre-processing step."
                                                                "Default is None.")
@click.option('--nfreq', default=50, help="Number of frequency bins")
@click.option('--fmin', default=0.1, help="Lowest frequency")
@click.option('--fmax', default=40., help="Highest frequency, which is clipped to the Nyquist value if larger")
@click.option('--lowpass-value', default=None, type=float, help="Lowpass filter value (Hz)")
@click.option('--highpass-value', default=None, type=float, help="Highpass filter value (Hz)")
@click.option('--freq-sampling', default='log',
              type=click.Choice(['linear', 'log']),
              help="Sampling method for frequency bins")
@click.option('--resample-log-freq', is_flag=True,
              help="Resample frequency bins in log-space. Only applicable if --freq-sampling is 'linear' and --spec-method is 'cwt2")
@click.option('--smooth-spectra-method', default='konno-ohmachi',
              type=click.Choice(['konno-ohmachi', 'none']),
              help="Smooth spectra using the Konno & Ohmachi method; 'none' to skip smoothing")
@click.option('--clip-fmin', default=0.3,
              help="Minimum clip frequency for master HVSR-curve. Only applicable if --clip-freq is given")
@click.option('--clip-fmax', default=50.,
              help="Maximum clip frequency for master HVSR-curve. Only applicable if --clip-freq is given")
@click.option('--clip-freq', is_flag=True,
              help="Clip master HVSR-curve to range specified by --clip-fmin and --clip-fmax")
@click.option('--master-curve-method', default='mean',
              type=click.Choice(['mean', 'median', 'geometric-average']),
              help="Method for computing master HVSR curve")
@click.option('--station-list', default='*', type=str,
              help="'A space-separated list of stations (within quotes) to process; default is '*', which "
                   "processes all available stations.",
              show_default=True)
@click.option('--start-time', default='1970-01-01T00:00:00',
              type=str,
              help="Date and time (in UTC format) to start from; default is year 1900.")
@click.option('--end-time', default='2100-01-01T00:00:00',
              type=str,
              help="Date and time (in UTC format) to stop at; default is year 2100.")
@click.option('--output-prefix', default='$network.$spec_method',
              type=str,
              help="Output file prefix; default is composed of network name and spectra method")
def process(asdf_source, network, spec_method, output_path, win_length,
            trigger_method, trigger_wlen, trigger_wlen_long, 
            trigger_threshold, trigger_lowpass_value, trigger_highpass_value,
            resample_rate, nfreq, fmin, fmax, lowpass_value, highpass_value,
            freq_sampling, resample_log_freq, smooth_spectra_method, clip_fmin,
            clip_fmax, clip_freq, master_curve_method, station_list, start_time,
            end_time, output_prefix):
    """
    ASDF_SOURCE: Path to text file containing paths to ASDF files\n
    NETWORK: Name of the network to process\n
    SPEC_METHOD: Method for computing spectra; ['single-taper', 'st', 'cwt2']. \n
    OUTPUT_PATH: Output folder \n
    """

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()
    proc_stations = defaultdict(list)

    if(output_prefix=='$network.$spec_method'):
        output_prefix = '{}.{}'.format(network, spec_method)
    # end if

    if(rank == 0):
        print('\n=== RunBatch Parameters ===\n')
        print(('ASDF-suource:            %s' % asdf_source))
        print(('Network:                 %s' % network))
        print(('Spec. Method:            %s' % spec_method))
        print(('Output-path:             %s' % output_path))
        print(('Win. Length:             %d (seconds)' % win_length))
        
        print(('Triggering method:       %s' % trigger_method))
        if(trigger_method=='zdetect'):
            print(('Trigger Win. Length:     %f (seconds)' % trigger_wlen))
        else:
            print(('Trigger Win. Length sta: %f (seconds)' % trigger_wlen))
            print(('Trigger Win. Length lta: %f (seconds)' % trigger_wlen_long))

        print(('Trigger threshold:       %3.2f' % trigger_threshold))
        if(trigger_lowpass_value):
            print(('Trigger lowpass value:   %3.2f' % trigger_lowpass_value))
        if(trigger_highpass_value):
            print(('Trigger highpass value:  %3.2f' % trigger_highpass_value))
        
        print(('nfreq:                   %d' % nfreq))
        print(('fmin:                    %f' % fmin))
        print(('fmax:                    %f' % fmax))
        if(lowpass_value):
            print(('lowpass_value:           %f (Hz)' % lowpass_value))
        if(highpass_value):
            print(('highpass_value:          %f (Hz)' % highpass_value))
        if(resample_rate):
            print(('resample_rate:          %f (Hz)' % resample_rate))
        print(('freq_sampling:           %s' % freq_sampling))
        print(('resample_log_freq:       %d' % resample_log_freq))
        print(('smooth_spectra_method:   %s' % smooth_spectra_method))
        print(('clip_freq:               %d' % clip_freq))
        if(clip_freq):
            print(('\tclip_fmin:             %d' % clip_fmin))
            print(('\tclip_fmax:             %d' % clip_fmax))
        print(('Start date and time:     %s' % start_time))
        print(('End date and time:       %s' % end_time))
        print('\n===========================\n')
    # end if

    # Removing '-' in options which are meant to avoid having to pass
    # a quoted string at the commandline
    if (spec_method == 'single-taper'): spec_method = 'single taper'
    if (master_curve_method == 'geometric-average'):
        master_curve_method = 'geometric average'

    if(smooth_spectra_method=='none'): smooth_spectra_method = None

    # Check input consistency
    if((resample_log_freq and spec_method != 'cwt2') or
       (resample_log_freq and freq_sampling != 'linear')):
        msg = "Error: Log-space frequency bin resampling can only be applied if --freq-sampling is 'linear' and --spec-method is 'cwt2'. Aborting..\n"
        sys.exit(msg)
    #end if

    # Triggering parameters
    triggering_options = {'method':trigger_method,
                          'trigger_wlen': trigger_wlen,
                          'trigger_wlen_long': trigger_wlen_long,
                          'trigger_threshold': trigger_threshold,
                          'trigger_lowpass_value': trigger_lowpass_value,
                          'trigger_highpass_value': trigger_highpass_value}

    # Get start and end times
    try:
        start_time = UTCDateTime(start_time)
        end_time = UTCDateTime(end_time)
    except:
        raise NameError('Failed to convert start or end time to UTCDateTime')
    # end try

    # check station names
    try:
        if(station_list != '*'):
            station_list = re.findall('\S+', station_list)
        # end if
    except:
        raise NameError('Invalid station-list..')
    # end try


    nfrequencies    = nfreq
    initialfreq     = fmin
    finalfreq       = fmax

    spectra_method  = spec_method
    lowest_freq     = clip_fmin
    highest_freq    = clip_fmax

    fds = FederatedASDFDataSet(asdf_source)
    if(rank == 0):
        stations = get_stations_to_process(fds, network, station_list)
        print("")
        print('Stations to process:')
        print(stations)
        print("")

        def split_list(lst, npartitions):
            k, m = divmod(len(lst), npartitions)
            return [lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(npartitions)]
        # end func

        proc_stations = split_list(stations, nproc)
    # end if

    # broadcast workload to all procs
    proc_stations = comm.scatter(proc_stations, root=0)

    # create temporary folder
    tempdir = None
    if(rank == 0):
        tempdir = os.path.join(output_path, str(uuid.uuid4()))
        os.makedirs(tempdir, exist_ok=True)
    # end if
    tempdir = comm.bcast(tempdir, root=0)

    spooled_results = []
    pbar = tqdm(total=len(proc_stations))
    for station in proc_stations:
        spooled_mat = None
        good_freqs = None

        net, sta = station.split('.')
        for st in waveform_iterator3c(fds, net, sta, start_time, end_time):
            if(not len(st)): continue # no data found

            pbar.set_description("Rank {}: {}: [{} - {}]".format(rank, station,
                                                                 st[0].stats.starttime.strftime('%Y-%m-%d'),
                                                                 st[0].stats.endtime.strftime('%Y-%m-%d')))

            freqs, hvsr_matrix = create_HVSR( st,
                                              spectra_method=spectra_method,
                                              spectra_options={'time_bandwidth': 3.5,
                                                               'number_of_tapers': None,
                                                               'quadratic': False,
                                                               'adaptive': True, 'nfft': None,
                                                               'taper': 'blackman'},
                                              window_length=win_length,
                                              bin_samples=nfrequencies,
                                              bin_sampling=freq_sampling,
                                              f_min=initialfreq,
                                              f_max=np.min([finalfreq, (1./st[0].stats.delta)*0.5]),
                                              triggering_options=triggering_options,
                                              lowpass_value = lowpass_value,
                                              highpass_value = highpass_value,
                                              new_sample_rate=resample_rate,
                                              resample_log_freq=resample_log_freq,
                                              smoothing=smooth_spectra_method )

            if(freqs is not None):
                good_freqs = freqs

                if (spooled_mat is None):
                    spooled_mat = SpooledMatrix(hvsr_matrix.shape[1],
                                                dtype=np.float32, max_size_mb=1048,
                                                prefix=station, dir=tempdir)
                # end if

                # append results to spooled_mat
                if(spooled_mat.ncols == hvsr_matrix.shape[1]):
                    for i in np.arange(len(hvsr_matrix)): spooled_mat.write_row(hvsr_matrix[i, :])
                else:
                    print('Warning: Discrepant hvsr_matrix shape found for {} at time {}. '
                          'Expected {} columns, found {}. Ignoring and moving along..'
                          .format(station, st[0].stats.starttime.strftime('%Y-%m-%d'),
                                  spooled_mat.ncols, hvsr_matrix.shape[1]))
                # end if
                #break
            # end if
        # end for
        spooled_results.append((station, spooled_mat, good_freqs))
        pbar.update()
    # end for
    pbar.close()

    # serialize processing of hvsr results
    pdf_files = []
    for irank in np.arange(nproc):
        if(rank==irank):
            for station, sm, freqs in spooled_results:
                result = generate_master_curve(station, tempdir, sm, freqs, master_curve_method, 0,
                                            clip_freq, lowest_freq, highest_freq,
                                            win_length)
                if(result):
                    net, sta = station.split('.')
                    pdf_fn, hvsr_freq, master_curve, lerr, uerr = result
                    pdf_files.append(pdf_fn)

                    # output hv
                    ofn = os.path.join(output_path, '{}.{}.hv.txt'.format(output_prefix, sta))
                    np.savetxt(ofn, np.column_stack((hvsr_freq, master_curve,
                                                     lerr, uerr)),
                               header='Freq[Hz] Amplitude Amp-1STD Amp+1STD')
                # end if
            # end for
        # end if
        comm.Barrier()
    # end for

    # gather pdf files
    comm.Barrier()
    pdf_files = comm.gather(pdf_files, root=0)

    if(rank == 0):
        # merge pdf files
        ofn = os.path.join(output_path, '{}.pdf'.format(output_prefix))
        pdf_files = [item for items in pdf_files for item in items]
        pdf_merge(pdf_files, ofn)

        os.removedirs(tempdir)

        print("Finishing...")
        print("HVSR:runbatch SUCCESS!")
    # end if
    comm.Barrier()
    del fds
# end func

if (__name__ == '__main__'):
    process()
# end if
