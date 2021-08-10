#!/usr/bin/env python
# coding: utf-8
"""Bulk plotting helper functions based on NetworkEventDataset
"""

# pylint: disable-all

import os
import math
import logging
import warnings

import click
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from tqdm.auto import tqdm

from seismic.network_event_dataset import NetworkEventDataset
from seismic.inversion.wavefield_decomp.runners import load_mcmc_solution
from seismic.stream_quality_filter import curate_seismograms
from seismic.stream_io import get_obspyh5_index


@click.command()
@click.option('--src-file', type=click.Path(exists=True, dir_okay=False), required=True,
              help='Event waveform source file for seismograms, generated using extract_event_traces.py script')
@click.option('--output-file', type=click.Path(dir_okay=False), required=True,
              help='Name of the output PDF file')
@click.option('--soln-file', type=click.Path(exists=True, dir_okay=False), required=True,
              help='Solution file that contains energy-per-event and filtering settings.')
def main(src_file, output_file, soln_file):
    """
    Example usage:

        python plot_network_events.py  --src-file OA_event_waveforms_for_rf_20170911T000036-20181128T230620_rev8.h5 \
         --output-file OA_event_seismograms.pdf
         --soln-file /g/data/ha3/am7399/shared/OA_wfdecomp_inversion/OA_wfd_out.h5

    :param src_file: File containing event waveforms
    :param output_file: Output pdf file name. Will be extended with a job id string.
    :param soln_file: Solution file that contains energy-per-event and filtering settings
    """

    if os.path.splitext(output_file)[1].lower() != '.pdf':
        output_file += '.pdf'
    # end if

    index = get_obspyh5_index(src_file, seeds_only=True)

    # Per page
    ncols = 2
    nrows = 8
    pagesize = (8.3, 8.3*1.414)  # A4
    channel_order = 'ZRT'

    bbstyle = dict(boxstyle="round", fc="w", alpha=0.5, linewidth=0.5)
    annot_fontsize = 5
    axes_fontsize = 6

    logger = logging.getLogger(__name__)

    soln_configs, job_id = load_mcmc_solution(soln_file)
    station_esu = {}
    for soln, config in soln_configs:
        station = config['station_id']
        station_esu[station] = (soln.esu, config)
    # end if
    out_basename, ext = os.path.splitext(output_file)
    out_basename += '_' + job_id
    output_file = out_basename + ext

    with PdfPages(output_file) as pdf:

        r_on_z = []
        t_on_z = []
        t_on_r = []
        z_cov_r = []
        z_cov_t = []
        r_cov_t = []
        energy_category = []

        pb = tqdm(index)
        for seedid in pb:
            if seedid not in station_esu:
                continue
            # end if
            soln_esu, sta_config = station_esu[seedid]

            net, sta, loc = seedid.split('.')
            pb.set_description(seedid)
            pb.write('Loading ' + seedid)

            with warnings.catch_warnings():
                warnings.simplefilter('ignore', category=FutureWarning)
                ned = NetworkEventDataset(src_file, net, sta, loc)
                # BEGIN PREPARATION & CURATION
                # Trim streams to time window
                su_energy_opts = sta_config["su_energy_opts"]
                time_window = su_energy_opts["time_window"]
                ned.apply(lambda stream: stream.trim(stream[0].stats.onset + time_window[0],
                                                     stream[0].stats.onset + time_window[1]))
                # Curate
                curation_opts = sta_config["curation_opts"]
                curate_seismograms(ned, curation_opts, logger)
                fs = su_energy_opts["sampling_rate"]
                # Downsample
                ned.apply(lambda stream: stream.filter('lowpass', freq=fs/2.0, corners=2, zerophase=True)
                          .interpolate(fs, method='lanczos', a=10))
                # END PREPARATION & CURATION
            # end with

            pb.write('Computing stats ' + seedid)
            db_evid = ned.station(sta)
            esu_evid = sta_config['event_ids']
            for esu in soln_esu:
                assert len(esu) == len(esu_evid)
            # end for
            esu_mean = np.array(soln_esu).mean(axis=0)
            assert len(esu_evid) == len(db_evid) == len(esu_mean)
            assert np.all(np.array(db_evid.keys()) == np.array(esu_evid))
            for i, stream in enumerate(db_evid.values()):
                tr_z = stream.select(component='Z')[0]
                tr_r = stream.select(component='R')[0]
                tr_t = stream.select(component='T')[0]
                # Ratio of RMS amplitudes
                rms_z = np.sqrt(np.nanmean(np.square(tr_z.data)))
                rms_r = np.sqrt(np.nanmean(np.square(tr_r.data)))
                rms_t = np.sqrt(np.nanmean(np.square(tr_t.data)))
                r_on_z.append(rms_r / rms_z)
                t_on_z.append(rms_t / rms_z)
                t_on_r.append(rms_t / rms_r)
                # Correlation coefficients
                corr_c = np.corrcoef([tr_z, tr_r, tr_t])
                z_cov_r.append(corr_c[0, 1])
                z_cov_t.append(corr_c[0, 2])
                r_cov_t.append(corr_c[1, 2])
                # Quality category based on energy value
                e = esu_mean[i]
                energy_category.append(0 if e < 2 else 1 if e < 4 else 2)
            # end for

            pb.write('Rendering ' + seedid)
            num_events = len(db_evid)
            num_per_page = nrows*ncols
            ev_ids = list(db_evid.keys())
            num_pages = math.ceil(num_events/num_per_page)
            for pagenum in range(num_pages):
                f = plt.figure(constrained_layout=False, figsize=pagesize)
                f.suptitle('{}.{} event seismograms (pg {}/{})'.format(net, sta, pagenum + 1, num_pages), y=0.98,
                           fontsize=11, va='top')
                pgspec = gridspec.GridSpec(ncols=ncols, nrows=nrows, figure=f, left=0.1, right=0.9, bottom=0.05,
                                           top=0.95, hspace=0.3, wspace=0.2)
                for i, evid in enumerate(ev_ids[pagenum*num_per_page:(pagenum + 1)*num_per_page]):
                    stream = db_evid[evid]
                    # Find max range so we can use a common scale
                    max_half_range = 0.0
                    for tr in stream:
                        max_half_range = max(max_half_range, 0.5*(np.nanmax(tr.data) - np.nanmin(tr.data)))
                    # end for
                    max_half_range *= 1.1

                    # Do plots
                    gs = pgspec[i % nrows, i//nrows].subgridspec(3, 1, hspace=0.0)
                    for tr in stream:
                        cha = tr.stats.channel[-1]
                        idx = channel_order.index(cha)
                        axn = f.add_subplot(gs[idx])
                        t = tr.times() - (tr.stats.onset - tr.stats.starttime)
                        tr_mean = np.nanmean(tr.data)
                        # Rasterize so that size is fairly constant, irrespective
                        # of the number of sample points per trace
                        axn.plot(t, tr.data, 'k', rasterized=True)
                        axn.set_ylim(tr_mean - max_half_range, tr_mean + max_half_range)
                        axn.tick_params(labelsize=axes_fontsize)
                        tag = '.'.join([tr.stats.network, tr.stats.station, tr.stats.location, tr.stats.channel])
                        axn.text(0.02, 0.90, tag, fontsize=annot_fontsize, ha='left', va='top', transform=axn.transAxes,
                                 bbox=bbstyle)
                        if idx == 0:
                            id_label = 'Event: {}'.format(evid)
                            axn.text(0.98, 0.90, id_label, fontsize=annot_fontsize, ha='right', va='top',
                                     transform=axn.transAxes, bbox=bbstyle)
                        # end if
                        if idx != 2:
                            axn.set_xticks([])
                            axn.xaxis.set_visible(False)
                            axn.xaxis.label.set_visible(False)
                        else:
                            axn.xaxis.label.set_size(axes_fontsize)
                            axn.text(0.5, 0.02, 'Onset: {}'.format(str(tr.stats.onset)), fontsize=annot_fontsize,
                                     ha='center', va='bottom', transform=axn.transAxes)
                        # end if
                        axn.yaxis.label.set_size(axes_fontsize)
                        # f.add_subplot(axn)
                    # end for
                # end for

                pdf.savefig(dpi=300, papertype='a4', orientation='portrait')
                plt.close()
            # end for
        # end for
        pb.close()

        # Cluster plot statistics.
        # Colour points by rank on energy spectrum for Tao method, to see
        # if waveforms with high energy (ones we want to remove in advance if
        # possible) are clustered in a way that could be filtered a priori.
        # See https://seaborn.pydata.org/generated/seaborn.pairplot.html
        r_on_z = np.array(r_on_z)
        t_on_z = np.array(t_on_z)
        t_on_r = np.array(t_on_r)
        qual = np.array(energy_category, dtype=int)

        rms_array = {'R/Z': r_on_z, 'T/Z': t_on_z, 'T/R': t_on_r, 'Energy': qual}
        df = pd.DataFrame.from_dict(rms_array)
        _p = sns.pairplot(df, hue='Energy', palette=sns.color_palette("Set2"),
                          plot_kws=dict(s=5, alpha=0.5, rasterized=True))
        _p.set(xlim=(0, 3), ylim=(0, 3))
        plt.suptitle('Ratios of channel RMS amplitudes')
        plt.tight_layout()
        pdf.savefig(dpi=300, papertype='a4', orientation='portrait')
        plt.close()

        z_cov_r = np.array(z_cov_r)
        z_cov_t = np.array(z_cov_t)
        r_cov_t = np.array(r_cov_t)

        cov_array = {'cov(Z,R)': z_cov_r, 'cov(Z,T)': z_cov_t, 'cov(R,T)': r_cov_t, 'Energy': qual}
        df = pd.DataFrame.from_dict(cov_array)
        _p = sns.pairplot(df, hue='Energy', palette=sns.color_palette("Set2"),
                          plot_kws=dict(s=5, alpha=0.5, rasterized=True))
        plt.suptitle('Cross-correlation coefficients between channels')
        plt.tight_layout()
        pdf.savefig(dpi=300, papertype='a4', orientation='portrait')
        plt.close()

    # end with
# end func


if __name__ == "__main__":
    main()
# end if
