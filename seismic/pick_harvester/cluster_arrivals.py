#!/bin/env python
"""
Description:
    Exports events and arrivals from the HDF5 output of ssst_relocate.py

References:

CreationDate:   11/08/22
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     11/08/22   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from seismic.pick_harvester.cluster_utils import NestedGrid
from seismic.pick_harvester.ssst_utils import SSST_Result
from seismic.pick_harvester.cluster_plot_utils import plot_before_cluster, plot_after_cluster
from collections import defaultdict
from tqdm import tqdm
import os
import click

def get_source_block_indices(ng:NestedGrid, sr:SSST_Result, imask):
    source_block = np.zeros(np.sum(imask), dtype='i4')
    bidx_dict = defaultdict(lambda: -1)
    for i, (elon, elat, edepth_km)  in enumerate(tqdm(zip(sr.elons[imask], sr.elats[imask],
                                                          sr.edepths_km[imask]),
                                                      desc='Finding source-block IDs')):
        if(elon < 0): elon += 360 # convert to being within [0, 360]
        ecolat = 90 - elat

        if((elon, ecolat, edepth_km) in bidx_dict.keys()):
            source_block[i] = bidx_dict[(elon, ecolat, edepth_km)]
        else:
            bidx = ng.get_block_index(elon, ecolat, edepth_km)
            source_block[i] = bidx
            bidx_dict[(elon, ecolat, edepth_km)] = bidx
        # end if
    # end for
    return source_block
# end func

def get_station_block_indices(ng:NestedGrid, sr:SSST_Result, imask):
    station_block = np.zeros(np.sum(imask), dtype='i4')
    bidx_dict = defaultdict(lambda: -1)
    for i, (slon, slat)  in enumerate(tqdm(zip(sr.slons[imask], sr.slats[imask]),
                                           desc='Finding station-block IDs')):
        if(slon < 0): slon += 360 # convert to being within [0, 360]
        scolat = 90 - slat

        if((slon, scolat) in bidx_dict.keys()):
            station_block[i] = bidx_dict[(slon, scolat)]
        else:
            bidx = ng.get_block_index(slon, scolat, 0)
            station_block[i] = bidx
            bidx_dict[(slon, scolat)] = bidx
        # end if
    # end for
    return station_block
# end func

def cluster(ng:NestedGrid, sr:SSST_Result, phases,
            p_residual_cutoff=5., s_residual_cutoff=10.,
            min_slope_ratio=5, only_paired_s=False):

    def cluster_helper(c_imask):
        source_block = get_source_block_indices(ng, sr, c_imask)
        station_block = get_station_block_indices(ng, sr, c_imask)

        cdict = defaultdict(list)
        indices = np.arange(len(sr.arrivals))
        for i, j, idx in tqdm(zip(source_block, station_block, indices[c_imask])):
            cdict[(i, j)].append([i, j, idx])
        # end for

        # convert dict-entries to numpy arrays
        for k in tqdm(cdict.keys()): cdict[k] = np.array(cdict[k])

        result = []
        for k, rows in tqdm(cdict.items()):
            ott = sr.ott[rows[:, 2]]

            med_idx = np.argwhere(ott == np.quantile(ott, 0.5, interpolation='nearest'))[0][0]

            src_block, sta_block, g_med_idx = rows[med_idx, :]
            tup = (src_block, sta_block, sr.residual[g_med_idx], sr.eorigin_ts[g_med_idx],
                   sr.elons[g_med_idx], sr.elats[g_med_idx], sr.edepths_km[g_med_idx],
                   sr.slons[g_med_idx], sr.slats[g_med_idx], sr.selevs_km[g_med_idx],
                   (sr.arrivals['arrival_ts'][g_med_idx] - sr.eorigin_ts[g_med_idx] - sr.tcorr[g_med_idx]),
                   sr.ecdists[g_med_idx], sr.arrivals['phase'][g_med_idx],
                   1 if sr.is_P[g_med_idx] else 2)

            result.append(tup)
        # end for

        fields = {'names': ['source_block', 'station_block', 'residual', 'eorigin_ts', 'elon', 'elat', 'edepth_km',
                            'slon', 'slat', 'selev_km', 'observed_tt', 'ecdist', 'phase', 'phase_type'],
                  'formats': ['i4', 'i4', 'f4', 'f8', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'U10', 'i4']}
        result = np.array(result, dtype=fields)
        return result
    # end func

    paired_S = None
    if(only_paired_s):
        print('Finding paired S arrivals..')
        paired_S = sr.find_paired_S()
    # end if

    cutoff = np.zeros(len(sr.arrivals))
    cutoff[sr.is_P] = p_residual_cutoff
    cutoff[sr.is_S] = s_residual_cutoff

    # Create an imask for phases to be clustered
    imask_phase = np.zeros(len(sr.arrivals), dtype='?')
    for ph in phases.split():
        ph = ph.strip().encode()

        imask_phase |= sr.phase == ph
    # end for

    # =======================================================
    # Create a combined imask:
    # 1. arrival is part of an event marked as being of
    #    good quality
    # 2. its residual lies within cutoff value
    # 3. if an automatic arrival, it should meed the
    #    slope-ratio criteria
    # 4. its phase must belong to the set of phases being
    #    used for clustering
    # =======================================================
    imask = (sr.equality) & \
            (np.fabs(sr.residual) < cutoff) & \
            (~sr.is_AUTO_arrival | (sr.is_AUTO_arrival & (sr.slope_ratio > min_slope_ratio))) & \
            (imask_phase)

    # cluster p-phases
    print('Clustering P phases..')
    cp = cluster_helper(imask & sr.is_P)

    print('Clustering S phases..')
    s_imask = paired_S if paired_S else sr.is_S
    cs = cluster_helper(imask & s_imask)

    return cp, cs
# end func

def save_clustered(ofn, p_clustered:np.ndarray, s_clustered:np.ndarray):
    fh = open(ofn, 'w')

    fmt = ''
    for name, dtype in p_clustered.dtype.descr:
        if ('i4' in dtype):
            fmt += '%10.0f '
        elif ('f8' in dtype):
            fmt += '%16.0f '
        elif ('f4' in dtype):
            fmt += '%16.7f '
        elif ('U' in dtype):
            fmt += '%s '
        # end if
    # end for

    np.savetxt(fh, p_clustered, fmt=fmt, header='   '.join(p_clustered.dtype.names), delimiter=' ')
    np.savetxt(fh, s_clustered, fmt=fmt, delimiter=' ')
    fh.close()
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('input_h5', type=click.Path(exists=True, dir_okay=False),
                required=True)
@click.argument('param_fn', type=click.Path(exists=True, dir_okay=False),
                required=True)
@click.argument('output_file_name', type=click.Path(exists=False, dir_okay=False),
                required=True)
@click.option('--phases', default='P Pg Pb S Sg Sb', help='A space-separated list of phases (case-sensitive and within '
                                                          'quotes) to cluster. Note that phases not in the list are '
                                                          'dropped.',
              show_default=True)
@click.option('--min-slope-ratio', default=5, help='Automatic arrivals with quality_measure_slope less than this '
                                                   'value are discarded prior to the clustering step',
              show_default=True)
@click.option('--p-residual-cutoff', default=5., help='P-arrivals with abs(residual) exceeding this cutoff are dropped.',
              show_default=True)
@click.option('--s-residual-cutoff', default=10., help='S-arrivals with abs(residual) exceeding this cutoff are dropped.',
              show_default=True)
@click.option('--outer-depth-refinement-factor', default=1, help='Depth-refinement factor for outer grid. '
              'Must be a power of 2.',
              show_default=True)
@click.option('--inner-depth-refinement-factor', default=1, help='Depth-refinement factor for inner grid. '
                                                                 'Must be a power of 2.',
              show_default=True)
@click.option('--only-paired-s', default=False, is_flag=True, help='S arrivals not accompanied by corresponding P '
                                                                   'arrivals are dropped before clustering.',
              show_default=True)
def process(input_h5, param_fn, output_file_name, phases, min_slope_ratio,
            p_residual_cutoff, s_residual_cutoff,
            outer_depth_refinement_factor, inner_depth_refinement_factor, only_paired_s):
    """
    INPUT_H5: hdf5 input (output of ssst_relocate.py)
    PARAM_FN: Grid parameterization file
    OUTPUT_FILE_NAME: name of output file
    """

    if('TXT' not in output_file_name.upper()): output_file_name = output_file_name + '.txt'
    pdf_output_fn, _ = os.path.splitext(output_file_name)
    pdf_output_fn += '.pdf'

    ng = NestedGrid(param_fn, outer_depth_refinement_factor=outer_depth_refinement_factor,
                    inner_depth_refinement_factor=inner_depth_refinement_factor)
    print('Instantiated clustering grid: {}'.format(ng))

    print('Loading SSST results..')
    sr = SSST_Result(input_h5)

    pdf = PdfPages(pdf_output_fn)

    print('Generating diagnostic plots before clustering..')
    plot_before_cluster(sr, ng, pdf, min_slope_ratio=min_slope_ratio)

    print('Clustering arrivals..')
    cp, cs = cluster(ng, sr, phases, p_residual_cutoff=p_residual_cutoff,
                     s_residual_cutoff=s_residual_cutoff, min_slope_ratio=min_slope_ratio,
                     only_paired_s=only_paired_s)

    print('Generating diagnostic plots after clustering..')
    plot_after_cluster(cp, cs, ng, pdf)

    pdf.close()

    print('Saving results..')
    save_clustered(output_file_name, cp, cs)
# end func

if __name__ == "__main__":
    process()
# end if