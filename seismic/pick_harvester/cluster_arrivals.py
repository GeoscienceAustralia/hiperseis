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
from seismic.pick_harvester.cluster_plot_utils import plot_before_cluster
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
            min_slope_ratio=5):

    def cluster_helper(imask):
        source_block = get_source_block_indices(ng, sr, imask)
        station_block = get_station_block_indices(ng, sr, imask)

        cdict = defaultdict(list)
        for i, j, r, eots, elon, elat, edepth_km, slon, slat, selev_km, ott, ecdist \
                in tqdm(zip(source_block, station_block,
                            sr.residual[imask], sr.eorigin_ts[imask],
                            sr.elons[imask], sr.elats[imask], sr.edepths_km[imask],
                            sr.slons[imask], sr.slats[imask], sr.selevs_km[imask],
                            (sr.arrival_ts[imask] - sr.eorigin_ts[imask] - sr.tcorr[imask]),
                            sr.ecdists[imask])):
            cdict[(i, j)].append([i, j, r, eots, elon, elat, edepth_km, slon, slat, selev_km, ott, ecdist])
        # end for

        # convert dict-entries to numpy arrays
        for k in cdict.keys(): cdict[k] = np.array(cdict[k])

        result = []
        for i, (k, item) in enumerate(tqdm(cdict.items())):
            item_residuals = item[:, 2]
            med_idx = np.argwhere(item_residuals == np.quantile(item_residuals, 0.5, interpolation='nearest'))[0][0]

            result.append(tuple(item[med_idx, :]))
        # end for

        fields = {'names': ['source_block', 'station_block', 'residual', 'eorigin_ts', 'elon', 'elat', 'edepth_km',
                            'slon', 'slat', 'selev_km', 'observed_tt', 'ecdist'],
                  'formats': ['i4', 'i4', 'f4', 'f8', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4']}
        result = np.array(result, dtype=fields)
        return result
    # end func

    print('Finding paired P and S arrivals..')
    paired = sr.find_paired()
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
    print('Clustering P hases..')
    cp = cluster_helper(imask & sr.is_P)
    print(len(cp))

    print('Clustering S hases..')
    cs = cluster_helper(imask & paired) # S phases must have accompanying P phases
    print(len(cs))
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('input_h5', type=click.Path(exists=True, dir_okay=False),
                required=True)
@click.argument('param_fn', type=click.Path(exists=True, dir_okay=False),
                required=True)
@click.argument('output_file_name', type=click.Path(exists=False, dir_okay=False),
                required=True)
@click.option('--phases', default='P Pg Pb S Sg Sb', help='A space-separated list of phases (within quotes) to cluster. '
                                                          'Note that phases not in the list are dropped.',
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
def process(input_h5, param_fn, output_file_name, phases, min_slope_ratio,
            p_residual_cutoff, s_residual_cutoff,
            outer_depth_refinement_factor, inner_depth_refinement_factor):
    """
    INPUT_H5: hdf5 input (output of ssst_relocate.py)
    PARAM_FN: Grid parameterization file
    OUTPUT_FILE_NAME: name of output file
    """

    if('CSV' not in output_file_name.upper()): output_file_name = output_file_name + '.csv'
    pdf_output_fn, _ = os.path.splitext(output_file_name)
    pdf_output_fn += '.pdf'

    ng = NestedGrid(param_fn, outer_depth_refinement_factor=outer_depth_refinement_factor,
                    inner_depth_refinement_factor=inner_depth_refinement_factor)
    print('Instantiated clustering grid: {}'.format(ng))

    print('Loading SSST results..')
    sr = SSST_Result(input_h5)

    pdf = PdfPages(pdf_output_fn)

    print('Generating diagnostic plots before clustering..')
    #plot_before_cluster(sr, ng, pdf, min_slope_ratio=min_slope_ratio)
    cluster(ng, sr, phases, p_residual_cutoff=p_residual_cutoff, s_residual_cutoff=s_residual_cutoff)
# end func

if __name__ == "__main__":
    process()
# end if