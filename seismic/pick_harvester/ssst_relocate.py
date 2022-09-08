#!/bin/env python
"""
Description:
    Relocates events and redefines arrival phases

References:

CreationDate:   10/08/22
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     10/08/22   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import matplotlib.pyplot as plt

from seismic.pick_harvester.ssst_relocator import SSSTRelocator
from ordered_set import OrderedSet as set
from matplotlib.backends.backend_pdf import PdfPages
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
import numpy as np
import os
import click
from seismic.receiver_fn.rf_plot_utils import pdf_merge
from seismic.pick_harvester.ssst_utils import get_iters, h5_to_named_array

AU_SHEAR_DISPLACEMENTS_SHP_FILE = os.path.dirname(os.path.abspath(__file__)) + \
                                  '/notebooks/data/ShearDisplacementLines2_5M.shp'

def plot_ssst_results_for_source(sr, src, h5_output_fn, pdf_output_fn,
                                 p_residual_cutoff, s_residual_cutoff):
    pdf = PdfPages(pdf_output_fn)
    iters = get_iters(h5_output_fn)

    event_attr = 'is_{}_event'.format(src.decode())
    arrival_attr = 'is_{}_arrival'.format(src.decode())

    is_src_event = getattr(sr, event_attr)
    is_src_arrival = getattr(sr, arrival_attr)
    is_relocated = np.sum(is_src_event) == np.sum(sr.events_imask[is_src_event])
    phase_types = {'P': sr.p_phases, 'S': sr.s_phases}

    # Maps and delta distributions
    if (1):
        if (is_relocated):
            reloc_events = h5_to_named_array(h5_output_fn, '{}/events'.format(iters[-1]))[is_src_event]
            orig_events = sr.orig_events[getattr(sr, event_attr)]
            converged_imask = reloc_events['event_quality']

            # maps and magnitude distributions
            crs = ccrs.PlateCarree()
            fig, axes = plt.subplots(1, 2, subplot_kw={'projection': crs})
            fig.set_size_inches(20, 20)

            for ax in axes:
                ax.coastlines('50m')

                if (src in [b'GA', b'GG']):
                    shape_feature = ShapelyFeature(Reader(AU_SHEAR_DISPLACEMENTS_SHP_FILE).geometries(),
                                                   crs, facecolor='none', edgecolor='k', lw=0.05, zorder=1)

                    ax.add_feature(shape_feature)
                    ax.set_extent([110, 160, -10, -45], crs)
                # end if
            # end for
            lons, lats, mags = orig_events['lon'], orig_events['lat'], orig_events['mag']

            axes[0].scatter(lons, lats, s=10, transform=crs, marker='o', c='r', edgecolor='none',
                            alpha=1, zorder=10)
            axes[1].scatter(reloc_events['lon'][converged_imask],
                            reloc_events['lat'][converged_imask], s=10, transform=crs,
                            marker='o', c='g', edgecolor='none', alpha=1, zorder=10)
            axes[0].set_title('Original Events')
            axes[1].set_title('Relocated Events')

            ax0 = fig.add_subplot(331)
            _ = ax0.hist(orig_events['mag'])
            ax0.text(0.7, 0.7, 'N: {}'.format(len(orig_events['mag'])), transform=ax0.transAxes)

            ax1 = fig.add_subplot(333)
            _ = ax1.hist(reloc_events['mag'][converged_imask])
            ax1.text(0.7, 0.7, 'N: {}'.format(len(reloc_events['mag'][converged_imask])),
                     transform=ax1.transAxes)

            ax0.set_xlabel('Magnitude'); ax1.set_xlabel('Magnitude')
            ax0.set_ylabel('Frequency'); ax1.set_ylabel('Frequency')

            pdf.savefig(dpi=300)
            plt.close()

            # delta
            fig, axes = plt.subplots(2, 2)
            fig.set_size_inches(10, 10)
            for ivar, var in enumerate(['lon', 'lat', 'depth_km', 'origin_ts']):
                ax = axes.flatten()[ivar]

                data = orig_events[var][converged_imask] - reloc_events[var][converged_imask]
                _ = ax.hist(data, bins=20)
                ax.set_xlabel('$\Delta$ %s' % (var.upper()))
                ax.set_ylabel('Frequency')

                ax.text(0.7, 0.7, 'N: {}'.format(np.sum(converged_imask)), transform=ax.transAxes)
                ax.text(0.7, 0.65, 'std: {:0.3f}'.format(np.std(data)), transform=ax.transAxes)
            # end for
            pdf.savefig(dpi=300)
            plt.close()
        # end if
    # end if

    # residual distributions for input data
    if (1):
        fig, axes = plt.subplots(2, 2)
        fig.set_size_inches(10, 10)

        phase = sr.orig_phases
        phase_imasks = [sr.is_P(phase),
                        sr.is_S(phase)]

        phc1 = {}
        phc2 = {}
        for p in phase_types['P']:
            phc1[p] = np.sum(p.encode() == phase[is_src_arrival])
        # end for
        for p in phase_types['S']:
            phc2[p] = np.sum(p.encode() == phase[is_src_arrival])
        # end for

        fig.suptitle('Distributions of residuals for input data\n\n'
                     'Phase counts P: {}\n'.format(phc1) +
                     'Phase counts S: {}\n'.format(phc2), size=18)

        PREEX_CUTOFF = 100
        k = 0
        for pt, pi in zip(phase_types.keys(), phase_imasks):
            # preexisting P/S arrivals
            p_imask = pi

            r1 = sr.r0[~sr.is_AUTO_arrival & is_src_arrival & p_imask]
            r2 = sr.r0[sr.is_AUTO_arrival & is_src_arrival & p_imask]

            r1 = r1[np.fabs(r1) < PREEX_CUTOFF]
            r2 = r2[np.fabs(r2) < PREEX_CUTOFF]

            _ = axes[k, 0].hist(r1, bins=50)
            _ = axes[k, 1].hist(r2, bins=50)
            axes[k, 0].text(0.7, 0.7, 'N: {}'.format(len(r1)), transform=axes[k, 0].transAxes)
            axes[k, 0].text(0.7, 0.65, 'std: {:0.3f}'.format(np.std(r1)), transform=axes[k, 0].transAxes)
            axes[k, 1].text(0.7, 0.7, 'N: {}'.format(len(r2)), transform=axes[k, 1].transAxes)
            axes[k, 1].text(0.7, 0.65, 'std: {:0.3f}'.format(np.std(r2)), transform=axes[k, 1].transAxes)

            axes[k, 0].set_xlabel('Residual [s]'); axes[k, 0].set_ylabel('Frequency');
            axes[k, 1].set_xlabel('Residual [s]'); axes[k, 1].set_ylabel('Frequency');
            axes[k, 0].set_title('Preexisting {} Arrivals'.format(pt))
            axes[k, 1].set_title('Automatic {} Arrivals'.format(pt))
            k += 1
        # end for
        plt.tight_layout()
        pdf.savefig(dpi=300)
        plt.close()
    # end if

    # residual and tcorr distribitions through ssst iterations
    if (1):
        for i in iters:
            tcorr_label = 'SSST Time-correction [s]' if i > 0 else 'SST Time-correction [s]'

            # if(i>0): continue

            res = h5_to_named_array(h5_output_fn, '{}/arrivals/residual'.format(i))
            tcorr = h5_to_named_array(h5_output_fn, '{}/arrivals/tcorr'.format(i))
            eq = h5_to_named_array(h5_output_fn, '{}/events/event_quality'.format(i))
            q = h5_to_named_array(h5_output_fn, '{}/arrivals/quality_measure_slope'.format(i))
            phase = h5_to_named_array(h5_output_fn, '{}/arrivals/phase'.format(i))

            cutoff_times = [p_residual_cutoff, s_residual_cutoff]
            phase_imasks = [sr.is_P(h5_to_named_array(h5_output_fn, '{}/arrivals/phase'.format(i))),
                            sr.is_S(h5_to_named_array(h5_output_fn, '{}/arrivals/phase'.format(i)))]

            for pt, ct, pi in zip(phase_types.keys(), cutoff_times, phase_imasks):
                # preexisting P/S arrivals
                imask = pi
                imask &= eq[sr.event_id_to_idx[sr.arrivals['event_id']]]
                imask &= np.fabs(res) < ct

                r = res[~sr.is_AUTO_arrival & is_src_arrival & imask]
                tc = tcorr[~sr.is_AUTO_arrival & is_src_arrival & imask]

                phc = {}
                for p in phase_types[pt]:
                    phc[p] = np.sum(p.encode() == phase[~sr.is_AUTO_arrival & is_src_arrival])
                # end for

                fig, axes = plt.subplots(1, 2)
                fig.set_size_inches(10, 5)
                fig.suptitle('Distributions of Preexisting {} residuals/time-corrections: iter {}\n\n'
                             'Phase counts: {}'.format(pt, i, phc), size=18)

                _ = axes[0].hist(r, bins=50)
                _ = axes[1].hist(tc, bins=50)

                axes[0].text(0.7, 0.7, 'N: {}'.format(len(r)), transform=axes[0].transAxes)
                axes[0].text(0.7, 0.65, 'std: {:0.3f}'.format(np.std(r)), transform=axes[0].transAxes)
                axes[1].text(0.7, 0.7, 'N: {}'.format(len(tc)), transform=axes[1].transAxes)
                axes[1].text(0.7, 0.65, 'std: {:0.3f}'.format(np.std(tc)), transform=axes[1].transAxes)
                axes[0].set_xlabel('Residual [s]');
                axes[0].set_ylabel('Frequency')
                axes[1].set_xlabel(tcorr_label);
                axes[1].set_ylabel('Frequency')

                fig.subplots_adjust(top=0.8)
                plt.tight_layout()
                pdf.savefig(dpi=300)
                plt.close()

                # automatic P/S arrivals
                phc = {}
                for p in phase_types[pt]:
                    phc[p] = np.sum(p.encode() == phase[sr.is_AUTO_arrival & is_src_arrival])
                # end for

                fig, axes = plt.subplots(2, 6)
                fig.set_size_inches(15, 5)
                fig.suptitle('Distributions of automatic {} residuals/time-corrections: iter {}\n\n'
                             'Phase counts: {}'.format(pt, i, phc), size=18)

                for j, cutoff in enumerate([-1, 1, 2, 3, 4, 5]):
                    title = 'All' if j == 0 else 'Slope-ratio > {}'.format(cutoff)
                    r = res[sr.is_AUTO_arrival & is_src_arrival & imask & (q > cutoff)]
                    tc = tcorr[sr.is_AUTO_arrival & is_src_arrival & imask & (q > cutoff)]

                    _ = axes[0, j].hist(r, bins=20)
                    _ = axes[1, j].hist(tc, bins=20)

                    axes[0, j].set_title(title)
                    axes[0, j].text(0.7, 0.7, 'N: {}'.format(len(r)), transform=axes[0, j].transAxes)
                    axes[0, j].text(0.7, 0.6, 'std: {:0.3f}'.format(np.std(r)), transform=axes[0, j].transAxes)
                    axes[1, j].text(0.7, 0.7, 'N: {}'.format(len(tc)), transform=axes[1, j].transAxes)
                    axes[1, j].text(0.7, 0.6, 'std: {:0.3f}'.format(np.std(tc)), transform=axes[1, j].transAxes)
                    axes[0, j].set_xlabel('Residual [s]');
                    axes[0, j].set_ylabel('Frequency')
                    axes[1, j].set_xlabel(tcorr_label);
                    axes[1, j].set_ylabel('Frequency')
                # end for
                plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
                fig.subplots_adjust(top=0.8)
                plt.tight_layout()
                pdf.savefig(dpi=300)
                plt.close()
                # end for
        # end for
    # end if
    pdf.close()
# end func

def plot_ssst_results(sr:SSSTRelocator, h5_output_fn:str, pdf_output_fn:str,
                      p_residual_cutoff=5, s_residual_cutoff=10):
    """
    @param sr:
    @param h5_output_fn:
    @param pdf_output_fn:
    @param p_residual_cutoff
    @param s_residual_cutoff
    @return:
    """
    sources = set(sr.events['source'])

    pdf_files = []
    for src in sources:
        pdf_ofn = os.path.join(sr._temp_dir, 'ssst.{}.pdf'.format(src.decode()))
        try:
            plot_ssst_results_for_source(sr, src, h5_output_fn, pdf_ofn,
                                         p_residual_cutoff, s_residual_cutoff)
            pdf_files.append(pdf_ofn)
        except Exception as e:
            print(str(e))
        # end try
    # end for

    pdf_merge(pdf_files, pdf_output_fn)
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('catalog_csv', type=click.Path(exists=True, dir_okay=False),
                required=True)
@click.argument('config', type=click.Path(exists=True, dir_okay=False),
                required=True)
@click.argument('output_file_name', type=click.Path(exists=False, dir_okay=False),
                required=True)
@click.option('--automatic-picks-p', type=click.Path(dir_okay=False), default=None,
              help='Automatic P picks in txt format (output of pick.py).', show_default=True)
@click.option('--automatic-picks-s', type=click.Path(dir_okay=False), default=None,
              help='Automatic S picks in txt format (output of pick.py).', show_default=True)
@click.option('--phases', default='P Pg Pb Pn S Sg Sb Sn',
              help='A space-separated list of phases (case-sensitive and within quotes) to process. '
                   'Note that phases not in the list are dropped.',
              show_default=True)
@click.option('--p-residual-cutoff', default=5., help='P-arrivals with abs(residual) exceeding this cutoff are dropped.',
              show_default=True)
@click.option('--s-residual-cutoff', default=10., help='S-arrivals with abs(residual) exceeding this cutoff are dropped.',
              show_default=True)
@click.option('--min-slope-ratio', default=5, help='Automatic arrivals with quality_measure_slope less than this '
                                                   'value are excluded from SST and SSST computations and the event '
                                                   'relocation procedure -- the phases of such arrivals are redefined '
                                                   'nonetheless.',
              show_default=True)
@click.option('--ssst-niter', type=int, default=5,
              help='Number of ssst iterations', show_default=True)
@click.option('--ball-radius-km', type=int, default=55,
              help='Radius of sphere in km used for calculating ssst corrections', show_default=True)
def process(catalog_csv, config, output_file_name, automatic_picks_p, automatic_picks_s, phases,
            p_residual_cutoff, s_residual_cutoff, min_slope_ratio, ssst_niter, ball_radius_km):
    """
    CATALOG_CSV: catalog in csv format
    CONFIG: config file in json format
    OUTPUT_FILE_NAME: name of output file
    """

    if('H5' not in output_file_name.upper()): output_file_name = output_file_name + '.h5'
    if(os.path.exists(output_file_name)):
        assert 0, 'Output file exists; Aborting..'
    # end if

    auto_pick_files = []
    auto_pick_phases = []
    if(automatic_picks_p):
        auto_pick_files.append(automatic_picks_p)
        auto_pick_phases.append('P')
    # end if
    if(automatic_picks_s):
        auto_pick_files.append(automatic_picks_s)
        auto_pick_phases.append('S')
    # end if

    sr = SSSTRelocator(catalog_csv,
                       auto_pick_files=auto_pick_files,
                       auto_pick_phases=auto_pick_phases,
                       config_fn=config,
                       phase_list=phases,
                       p_residual_cutoff=p_residual_cutoff,
                       s_residual_cutoff=s_residual_cutoff)

    sr.ssst_relocate(ssst_niter=ssst_niter, ball_radius_km=ball_radius_km,
                     min_slope_ratio=min_slope_ratio, output_fn=output_file_name)

    # generate plots on master rank
    if(sr.rank == 0):
        print('SSST-relocation SUCCESS..')
        print('Generating pdf report..')
        pdf_output_file_name, _ = os.path.splitext(output_file_name)
        pdf_output_file_name += '.pdf'

        plot_ssst_results(sr, output_file_name, pdf_output_file_name,
                          p_residual_cutoff, s_residual_cutoff)
    # end if
# end func

if __name__ == "__main__":
    process()
# end if
