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

import matplotlib.pyplot as plt
from ordered_set import OrderedSet as set
from matplotlib.backends.backend_pdf import PdfPages
from seismic.pick_harvester.ssst_relocate import h5_to_named_array, get_iters
import cartopy.crs as ccrs
import numpy as np
from pyproj import Geod
from geographiclib.geodesic import Geodesic
from tqdm import tqdm

import os
import click

def plot_results(input_h5, pdf_output_fn, min_slope_ratio=5):
    # get final iteration
    iter = get_iters(input_h5)[-1]

    pdf = PdfPages(pdf_output_fn)

    # distance ranges for which plots are to be generated
    distance_ranges = [[0, 10], [10, 20]] # degrees
    geod = Geod(a=180 / np.pi, f=0)
    geodesic = Geodesic.WGS84

    arrivals = h5_to_named_array(input_h5, '{}/arrivals'.format(iter))
    events = h5_to_named_array(input_h5, '{}/events'.format(iter))
    equality = h5_to_named_array(input_h5, '{}/events/event_quality'.format(iter))
    event_id_to_idx = h5_to_named_array(input_h5, '0/events/event_id_to_idx')
    residual = h5_to_named_array(input_h5, '{}/arrivals/residual'.format(iter))
    is_AUTO_arrival = h5_to_named_array(input_h5, '0/arrivals/is_AUTO_arrival')

    elons = events['lon'][event_id_to_idx[arrivals['event_id']]]
    elats = events['lat'][event_id_to_idx[arrivals['event_id']]]
    slons = arrivals['lon']
    slats = arrivals['lat']

    # define 0.2x0.2 deg grid
    au_sx, au_ex, dx = 110, 155, 0.2  # Australasian region (lon bounds)
    au_sy, au_ey, dy = -45, -10, 0.2  # Australasian region (lat bounds)

    nx = np.int_(np.ceil((au_ex - au_sx) / dx)) + 1
    ny = np.int_(np.ceil((au_ey - au_sy) / dy)) + 1
    res = int(1 / float(dx))  # per deg

    gx, gy = np.meshgrid(np.linspace(au_sx, au_ex, nx), np.linspace(au_sy, au_ey, ny), indexing='ij')

    for i, drange in enumerate(distance_ranges):
        print('Processing distance range {}: '.format(distance_ranges[i]))
        ge = np.zeros(gx.shape)


        local_distance = 20  # deg

        azims, _, ecdists = geod.inv(elons, elats, slons, slats)

        eqfilt = equality[event_id_to_idx[arrivals['event_id']]]

        # spatial clip
        au_filt = (elons > au_sx) & (elons < au_ex) & \
                  (elats > au_sy) & (elats < au_ey) & \
                  (slons > au_sx) & (slons < au_ex) & \
                  (slats > au_sy) & (slats < au_ey) & \
                  ((ecdists > drange[0]) & (ecdists < drange[1]))


        gzp = np.zeros(gx.shape)
        gzs = np.zeros(gx.shape)

        # ===========================================================
        # Process P arrivals
        # ===========================================================
        peset = set()
        psset = set()
        filt = (eqfilt) & (au_filt) & (arrivals['phase'].astype('S1') == b'P') & \
               (~is_AUTO_arrival | (is_AUTO_arrival & (arrivals['quality_measure_slope'] > min_slope_ratio)))

        p_elons = elons[filt]
        p_elats = elats[filt]
        p_slons = slons[filt]
        p_slats = slats[filt]
        for i in tqdm(np.arange(len(p_elons)), desc='P-arrivals: '):
            elon, elat, slon, slat = p_elons[i], p_elats[i], \
                                     p_slons[i], p_slats[i]

            peset.add((elon, elat))
            psset.add((slon, slat))

            g = geodesic.Inverse(elat, elon, slat, slon)
            l = geodesic.Line(g['lat1'], g['lon1'], g['azi1'])
            num = int(np.ceil(g['a12'] * res))

            for j in range(num + 1):
                pos = l.ArcPosition(j * g['a12'] / num)
                cy, cx = pos['lat2'], pos['lon2']

                cxi = np.int_((cx - au_sx) / dx)
                cyi = np.int_((cy - au_sy) / dy)

                gzp[cxi, cyi] += 1
                # end for
        # end for

        print('P-arrivals: Unique #events (%d), Unique #stations (%d)' % (len(peset), len(psset)))

        # ===========================================================
        # Process S arrivals
        # ===========================================================
        seset = set()
        ssset = set()
        filt = eqfilt & au_filt & (arrivals['phase'].astype('S1') == b'S') & \
               (~is_AUTO_arrival | (is_AUTO_arrival & (arrivals['quality_measure_slope'] > min_slope_ratio)))

        s_elons = elons[filt]
        s_elats = elats[filt]
        s_slons = slons[filt]
        s_slats = slats[filt]

        for i in tqdm(np.arange(len(s_elons)), desc='S-arrivals'):
            elon, elat, slon, slat = s_elons[i], s_elats[i], \
                                     s_slons[i], s_slats[i]

            seset.add((elon, elat))
            ssset.add((slon, slat))

            g = geodesic.Inverse(elat, elon, slat, slon)
            l = geodesic.Line(g['lat1'], g['lon1'], g['azi1'])
            num = int(np.ceil(g['a12'] * res))

            for j in range(num + 1):
                pos = l.ArcPosition(j * g['a12'] / num)
                cy, cx = pos['lat2'], pos['lon2']

                cxi = np.int_((cx - au_sx) / dx)
                cyi = np.int_((cy - au_sy) / dy)

                gzs[cxi, cyi] += 1
                # end for
        # end for

        print('S-arrivals: Unique #events (%d), Unique #stations (%d)' % (len(seset), len(ssset)))
        # ===========================================================
        # Plot results
        # ===========================================================
        # Plot p-coverage
        fig = plt.figure()

        fig.set_size_inches(20, 10)
        cax1 = fig.add_axes([0.1, 0.1, 0.4, 0.1])
        cax2 = fig.add_axes([0.5, 0.1, 0.4, 0.1])
        cax1.set_visible(False)
        cax2.set_visible(False)

        llcrnrlat = -45
        urcrnrlat = -10
        llcrnrlon = 110
        urcrnrlon = 155

        ax1 = fig.add_subplot(1, 2, 1, projection=ccrs.PlateCarree())
        ax1.set_extent([llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat], crs=ccrs.PlateCarree())
        ax1.coastlines()

        cbi = ax1.pcolormesh(gx, gy, gzp, vmin=1, vmax=np.max(gzp), cmap=plt.get_cmap('gist_heat_r', 50))

        cbi.cmap.set_under('#99ccff')
        cbar = fig.colorbar(cbi, ax=cax1, format='%d', extend='min',
                            orientation='horizontal')
        cbar.set_label("Hitcount")
        ax1.set_title('P-arrivals')

        # Plot s-coverage
        ax2 = fig.add_subplot(1, 2, 2, projection=ccrs.PlateCarree())
        ax2.set_extent([llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat], crs=ccrs.PlateCarree())
        ax2.coastlines()

        cbi = ax2.pcolormesh(gx, gy, gzs, vmin=1, vmax=np.max(gzs), cmap=plt.get_cmap('gist_heat_r', 50))

        cbi.cmap.set_under('#99ccff')
        cbar = fig.colorbar(cbi, ax=cax2, format='%d', extend='min',
                            orientation='horizontal')
        cbar.set_label("Hitcount")
        ax2.set_title('S-arrivals')

        fig.suptitle('Raypath (epicentral distance {} deg) Coverage in the Australasian Region'.format(drange), size=18)
        pdf.savefig(dpi=300)
        plt.close()
    # end for

    # ===========================================================
    # Plot distributions of residuals for all preexisting arrivals
    # to be exported
    # ===========================================================
    P_CUTOFF = 5
    S_CUTOFF = 10
    pfilt = (eqfilt) & (arrivals['phase'].astype('S1') == b'P') & \
            (~is_AUTO_arrival) & (np.fabs(residual) < P_CUTOFF)
    sfilt = (eqfilt) & (arrivals['phase'].astype('S1') == b'S') & \
            (~is_AUTO_arrival) & (np.fabs(residual) < S_CUTOFF)

    fig, axes = plt.subplots(1, 2)
    fig.set_size_inches(10, 5)

    _ = axes[0].hist(residual[pfilt], bins=20)
    _ = axes[1].hist(residual[sfilt], bins=20)

    axes[0].set_xlabel('Residual [s]'); axes[0].set_ylabel('Frequency')
    axes[1].set_xlabel('Residual [s]'); axes[1].set_ylabel('Frequency')
    axes[0].set_title('P-arrivals'); axes[1].set_title('S-arrivals')

    axes[0].text(0.7, 0.7, 'N: {}'.format(np.sum(pfilt)), transform=axes[0].transAxes)
    axes[0].text(0.7, 0.65, 'std: {:0.3f}'.format(np.std(residual[pfilt])), transform=axes[0].transAxes)
    axes[1].text(0.7, 0.7, 'N: {}'.format(np.sum(sfilt)), transform=axes[1].transAxes)
    axes[1].text(0.7, 0.65, 'std: {:0.3f}'.format(np.std(residual[sfilt])), transform=axes[1].transAxes)

    fig.suptitle('All Preexisting Arrivals', fontsize=18)
    pdf.savefig(dpi=300)
    plt.close()

    # ===========================================================
    # Plot distributions of residuals for all automatic
    # arrivals to be exported
    # ===========================================================
    pfilt = (eqfilt) & (arrivals['phase'].astype('S1') == b'P') & \
            (is_AUTO_arrival & (arrivals['quality_measure_slope'] > min_slope_ratio)) & \
            (np.fabs(residual) < P_CUTOFF)
    sfilt = (eqfilt) & (arrivals['phase'].astype('S1') == b'S') & \
            (is_AUTO_arrival & (arrivals['quality_measure_slope'] > min_slope_ratio)) & \
            (np.fabs(residual) < S_CUTOFF)

    fig, axes = plt.subplots(1, 2)
    fig.set_size_inches(10, 5)

    _ = axes[0].hist(residual[pfilt], bins=20)
    _ = axes[1].hist(residual[sfilt], bins=20)

    axes[0].set_xlabel('Residual [s]'); axes[0].set_ylabel('Frequency')
    axes[1].set_xlabel('Residual [s]'); axes[1].set_ylabel('Frequency')
    axes[0].set_title('P-arrivals'); axes[1].set_title('S-arrivals')

    axes[0].text(0.7, 0.7, 'N: {}'.format(np.sum(pfilt)), transform=axes[0].transAxes)
    axes[0].text(0.7, 0.65, 'std: {:0.3f}'.format(np.std(residual[pfilt])), transform=axes[0].transAxes)
    axes[1].text(0.7, 0.7, 'N: {}'.format(np.sum(sfilt)), transform=axes[1].transAxes)
    axes[1].text(0.7, 0.65, 'std: {:0.3f}'.format(np.std(residual[sfilt])), transform=axes[1].transAxes)

    fig.suptitle('All Automatic Arrivals (slope-ratio > {})'.format(min_slope_ratio), fontsize=18)
    pdf.savefig(dpi=300)
    plt.close()

    pdf.close()
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('input_h5', type=click.Path(exists=True, dir_okay=False),
                required=True)
@click.argument('output_file_name', type=click.Path(exists=False, dir_okay=False),
                required=True)
@click.option('--min-slope-ratio', default=5, help='Automatic arrivals with quality_measure_slope less than this '
                                                   'value are discarded.',
              show_default=True)
def process(input_h5, output_file_name, min_slope_ratio):
    """
    INPUT_H5: hdf5 input (output of ssst_relocate.py)
    OUTPUT_FILE_NAME: name of output file
    """

    if('CSV' not in output_file_name.upper()): output_file_name = output_file_name + '.csv'
    pdf_output_file_name, _ = os.path.splitext(output_file_name)
    pdf_output_file_name += '.pdf'

    plot_results(input_h5, pdf_output_file_name, min_slope_ratio)
# end func

if __name__ == "__main__":
    process()
# end if