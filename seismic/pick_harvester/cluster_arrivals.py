#!/bin/env python
"""
Description:
    Clusters arrivals from the HDF5 output of ssst_relocate.py

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
    for i, (elon, elat, edepth_km) in enumerate(tqdm(zip(sr.elons[imask], sr.elats[imask],
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
    for i, (slon, slat) in enumerate(tqdm(zip(sr.slons[imask], sr.slats[imask]),
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
            min_slope_ratio=5,
            outer_grid_min_rays=10,
            inner_grid_min_rays=3,
            cluster_only_automatic=False,
            match_p_s=False):

    def cluster_helper(c_imask):
        source_blocks = get_source_block_indices(ng, sr, c_imask)
        station_blocks = get_station_block_indices(ng, sr, c_imask)

        cdict = defaultdict(list)
        indices = np.arange(len(sr.arrivals))
        for i, j, idx in zip(source_blocks, station_blocks, indices[c_imask]):
            cdict[(i, j)].append(idx)
        # end for

        # convert dict-entries to numpy arrays
        for k in cdict.keys(): cdict[k] = np.array(cdict[k])

        result = []
        both_inner = 0
        one_inner = 0
        both_outer = 0
        for k, indices in tqdm(cdict.items()):
            src_block, sta_block = k
            ctt = sr.corrected_travel_time[indices]

            is_src_block_inner = ng.is_inner_block(src_block)
            is_sta_block_inner = ng.is_inner_block(sta_block)

            if(is_src_block_inner and is_sta_block_inner): both_inner += 1
            elif(is_src_block_inner or is_sta_block_inner): one_inner += 1
            else: both_outer += 1

            min_rays = inner_grid_min_rays
            if((not is_src_block_inner) and (not is_sta_block_inner)):
                min_rays = outer_grid_min_rays
            # end if

            a_imask = None
            if(cluster_only_automatic and (is_src_block_inner or is_sta_block_inner)):
                # only these arrivals are to be clustered
                a_imask = sr.is_AUTO_arrival[indices]

                # gather preexisting arrivals that are not to be clustered
                for idx in indices[~a_imask]:
                    tup = (src_block, sta_block, sr.residual[idx], sr.eorigin_ts[idx],
                           sr.elons[idx], sr.elats[idx], sr.edepths_km[idx],
                           sr.slons[idx], sr.slats[idx], sr.selevs_km[idx],
                           sr.corrected_travel_time[idx],
                           sr.ecdists[idx], sr.arrivals['phase'][idx],
                           1 if sr.is_P[idx] else 2,
                           sr.arrivals['event_id'][idx],
                           sr.arrivals['net'][idx],
                           sr.arrivals['sta'][idx],
                           sr.is_AUTO_arrival[idx])
                    result.append(tup)
                # end for
            else:
                a_imask = np.bool_(np.ones(len(indices)))
            # end if

            if(len(ctt[a_imask])):
                if (len(ctt[a_imask]) < min_rays): continue  # must have at least minimum number of rays

                med_idx = np.argwhere(ctt[a_imask] == np.quantile(ctt[a_imask], 0.5,
                                                                  interpolation='nearest'))[0][0]

                g_med_idx = indices[a_imask][med_idx]
                tup = (src_block, sta_block, sr.residual[g_med_idx], sr.eorigin_ts[g_med_idx],
                       sr.elons[g_med_idx], sr.elats[g_med_idx], sr.edepths_km[g_med_idx],
                       sr.slons[g_med_idx], sr.slats[g_med_idx], sr.selevs_km[g_med_idx],
                       sr.corrected_travel_time[g_med_idx],
                       sr.ecdists[g_med_idx], sr.arrivals['phase'][g_med_idx],
                       1 if sr.is_P[g_med_idx] else 2,
                       sr.arrivals['event_id'][g_med_idx],
                       sr.arrivals['net'][g_med_idx],
                       sr.arrivals['sta'][g_med_idx],
                       sr.is_AUTO_arrival[g_med_idx])
                result.append(tup)
            # end if

        # end for

        #print('\n **both inner: {}, one inner: {}, both outer: {} **\n'.format(both_inner,
        #                                                                       one_inner,
        #                                                                       both_outer))

        fields = {'names': ['source_block', 'station_block', 'residual', 'eorigin_ts', 'elon', 'elat', 'edepth_km',
                            'slon', 'slat', 'selev_km', 'observed_tt', 'ecdist', 'phase', 'phase_type', 'event_id',
                            'net', 'sta', 'is_automatic'],
                  'formats': ['i4', 'i4', 'f4', 'f8', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4',
                              'f4', 'f4', 'U10', 'i4', 'i4', 'U10', 'U10', 'i4']}
        result = np.array(result, dtype=fields)
        return result
    # end func

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
    # 3. if an automatic arrival, it should meet the
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
    cs = cluster_helper(imask & sr.is_S)

    if(match_p_s):
        def filter_record_arrays(from_array: np.ndarray, ref_array: np.ndarray):
            """
            Drops entries from 'from_array' where the source-station cell-pair is
            not found in 'ref_array'.
            @param from_array:
            @param ref_array:
            @return: Trimmed version of 'from_array'
            """
            rset = set()
            for i in np.arange(len(ref_array)):
                rset.add((ref_array['source_block'][i], ref_array['station_block'][i]))
            # end for

            rows = []
            for i in np.arange(len(from_array)):
                if((from_array['source_block'][i], from_array['station_block'][i]) in rset):
                    rows.append(from_array[i])
                # end if
            # end for

            result_array = np.zeros(len(rows), dtype=from_array.dtype)
            for i, row in enumerate(rows):
                result_array[i] = np.array(row, dtype=from_array.dtype)
            # end for

            return result_array
        # end func

        # drop P arrivals if S arrivals are not found in corresponding source-station cell-pair
        cp = filter_record_arrays(cp, cs)

        # drop S arrivals if P arrivals are not found in corresponding source-station cell-pair
        cs = filter_record_arrays(cs, cp)
    # end if

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
@click.option('--outer-grid-min-rays', default=10,
              help='Minimum number of rays that a source-receiver cell-pair, within the outer grid, '
                   'must have to participate in the clustering process.',
              show_default=True)
@click.option('--inner-grid-min-rays', default=3,
              help='Minimum number of rays that a source-receiver cell-pair, with either the source '
                   'or the receiver or both within the inner grid, must have to participate in the '
                   'clustering process.',
              show_default=True)
@click.option('--cluster-only-automatic', default=False, is_flag=True,
              help='Generates summary rays for only automatically picked arrivals, when (1) either the '
                   'source or the receiver or both are within the inner grid -- all preexisting '
                   'arrivals are included in such cases, which helps increase ray-path coverage '
                   'within the inner grid. When both source and receiver are in the outer grid (2), all '
                   'arrivals associated with such source-receiver cell pairs are clustered.'
                   'Note that with this option, for cases in (1), the value of parameter '
                   '--inner-grid-min-rays will apply to automatically picked arrivals only.',
              show_default=True)
@click.option('--match-p-s', default=False, is_flag=True,
              help='After the clustering step, S arrivals not accompanied by corresponding P arrivals '
                   '(matched by cell-ids) are dropped.',
              show_default=True)
def process(input_h5, param_fn, output_file_name, phases, min_slope_ratio,
            p_residual_cutoff, s_residual_cutoff,
            outer_depth_refinement_factor, inner_depth_refinement_factor,
            outer_grid_min_rays, inner_grid_min_rays,
            cluster_only_automatic, match_p_s):
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
                     outer_grid_min_rays=outer_grid_min_rays,
                     inner_grid_min_rays=inner_grid_min_rays,
                     cluster_only_automatic=cluster_only_automatic,
                     match_p_s=match_p_s)

    print('Generating diagnostic plots after clustering..')
    plot_after_cluster(cp, cs, ng, phases, pdf)

    pdf.close()

    print('Saving results..')
    save_clustered(output_file_name, cp, cs)
# end func

if __name__ == "__main__":
    process()
# end if