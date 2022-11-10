"""
Description:

References:

CreationDate:   15/03/22
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     15/03/22   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import os
from collections import defaultdict
import tqdm.auto as tqdm

import seismic.receiver_fn.rf_util as rf_util
from seismic.receiver_fn.rf_plot_utils import pdf_merge
from seismic.stream_io import get_obspyh5_index

from seismic.network_event_dataset import NetworkEventDataset
from seismic.swp_station_orientations import analyze_station_orientations as swp_station_orientations, load_grv
from seismic.rf_station_orientations import analyze_station_orientations as rf_station_orientations
import logging
import click
import uuid, json
import cartopy.crs as ccrs
import numpy as np
from scipy.stats import circmean as cmean
import matplotlib
from mpi4py import MPI

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from shutil import rmtree

logging.basicConfig()

paper_size_A4 = (8.27, 11.69)  # inches

def get_station_coords(ned):
    #assuming there is only one stations in ned
    result = defaultdict(list)
    for sta, evid, stream in ned:
        result['.'.join([ned.network, sta])] = [stream[0].stats.station_longitude, stream[0].stats.station_latitude]
        break
    # end for

    return result
# end func

def plot_summary(station_coords_dict, corrections_dict, output_fn):
    minLon = 1e32
    maxLon = -1e32
    minLat = 1e32
    maxLat = -1e32

    coords = defaultdict()
    lons = []
    for k, v in station_coords_dict.items():
        lon, lat = v[0], v[1]

        minLon = min(lon, minLon)
        maxLon = max(lon, maxLon)
        minLat = min(lat, minLat)
        maxLat = max(lat, maxLat)

        coords[k] = [lon, lat]
        lons.append(lon)
    # end for

    minLon -= 0.5
    maxLon += 0.5
    minLat -= 0.5
    maxLat += 0.5

    lons = np.array(lons)
    clon = cmean(lons, high=180, low=-180)

    crs = ccrs.PlateCarree(central_longitude=clon)

    fig = plt.figure(figsize=(paper_size_A4[0], paper_size_A4[1]))
    ax1 = plt.subplot(2, 1, 1, projection=crs)
    ax2 = plt.subplot(2, 1, 2, projection=crs)

    for ax, ttext in zip([ax1, ax2], ['Azimuth corrections (RF)', 'Azimuth corrections (Surface-wave Polarization)']):
        # draw coastlines.
        ax.coastlines('50m')
        #ax.set_extent([minLon, maxLon, minLat, maxLat], crs)
        gl = ax.gridlines(draw_labels=True,
                          linewidth=1, color='gray',
                          alpha=0.5, linestyle='--')
        ax.set_title(ttext, fontsize=20, pad=30)
    # end for

    # find cmap range
    diffs = []
    rf_corrs = corrections_dict['rf']
    swp_corrs = corrections_dict['swp']
    for nsl in coords.keys():
        lon, lat = coords[nsl]

        corr_rf = None
        corr_swp = None
        try:
            corr_rf = rf_corrs[nsl]['azimuth_correction']
            corr_swp = swp_corrs[nsl]['azimuth_correction']

            while (corr_rf > 180): corr_rf -= 360
            while (corr_rf < -180): corr_rf += 360

            while (corr_swp > 180): corr_swp -= 360
            while (corr_swp < -180): corr_swp += 360

        except:
            continue
        # end try
        diffs.append(np.fabs(corr_rf - corr_swp))
    # end for

    if(len(diffs)):
        cmap = matplotlib.cm.jet
        norm = matplotlib.colors.Normalize(0, np.max(diffs))

        # plot stations
        for nsl in coords.keys():
            net, sta, loc = nsl.split('.')
            lon, lat = coords[nsl]

            corr_rf = None
            corr_swp = None
            try:
                corr_swp = swp_corrs[nsl]['azimuth_correction']
                corr_rf = rf_corrs[nsl]['azimuth_correction']

                while (corr_rf > 180): corr_rf -= 360
                while (corr_rf < -180): corr_rf += 360

                while (corr_swp > 180): corr_swp -= 360
                while (corr_swp < -180): corr_swp += 360

            except:
                continue
            # end try

            for ax, corr in zip([ax1, ax2], [corr_rf, corr_swp]):
                # print (np.fabs(corr_rf-corr_swp))
                color = cmap(norm(np.fabs(corr_rf - corr_swp)))
                px, py = lon-clon, lat
                pxl, pyl = lon-clon + 0.02, lat - 0.1
                ax.scatter(px, py, 2, transform=crs, marker='o', c='g', edgecolor='none', zorder=10)
                ax.annotate(sta, xy=(pxl, pyl), fontsize=3)

                ux = np.cos(np.radians(corr))
                uy = np.sin(np.radians(corr))

                # print(netsta, corr, ux, uy)

                ax.quiver(px, py, -uy, ux, transform=crs, scale_units='inches',
                          color=color, scale=5, width=0.002, pivot='middle')
            # end for
        # end for
        cbax = fig.add_axes([0.35, 0.5, 0.3, 0.01])

        cb = matplotlib.colorbar.ColorbarBase(cbax, cmap=cmap, norm=norm, orientation='horizontal')
        cbax.set_title("Angular difference between methods [°]")
    # end for
    plt.tight_layout(h_pad=5)
    plt.savefig(output_fn, dpi=300)
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('src-h5-event-file', type=click.Path(exists=True, dir_okay=False),
                required=True)
@click.argument('network', type=str, required=True)
@click.option('--output-basename', type=click.Path(dir_okay=False),
              help='Output file basename to store results in JSON format and plots in pdf format')
@click.option('--station-list', default='*', help='A space-separated list of stations (within quotes) to process.', type=str,
              show_default=True)
def main(src_h5_event_file, network, output_basename, station_list):
    """
    Run station orientation checks.

    Example usage::

    """
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()
    proc_hdfkeys = None
    rf_h5_root = 'waveforms/P'
    sw_h5_root = 'waveforms/SW'

    tempdir = None
    if(rank == 0):
        # retrieve all available hdf_keys
        proc_hdfkeys = set(get_obspyh5_index(src_h5_event_file, seeds_only=True, root=rf_h5_root))
        proc_hdfkeys = proc_hdfkeys.union(set(get_obspyh5_index(src_h5_event_file, seeds_only=True, root=sw_h5_root)))
        proc_hdfkeys = list(proc_hdfkeys)

        # trim stations to be processed based on the user-provided network- and station-list
        proc_hdfkeys = rf_util.trim_hdf_keys(proc_hdfkeys, network, station_list)

        # split work-load over all procs
        proc_hdfkeys = rf_util.split_list(proc_hdfkeys, nproc)
        tempdir = os.path.join(os.path.dirname(output_basename), str(uuid.uuid4()))
        os.makedirs(tempdir, exist_ok=True)
    # end if

    tempdir = comm.bcast(tempdir, root=0)
    proc_hdfkeys = comm.bcast(proc_hdfkeys, root=0)

    local_results_rf = defaultdict(dict)
    local_results_swp = defaultdict(dict)
    pbar = tqdm.tqdm(total=len(proc_hdfkeys[rank]))
    grv_dict = load_grv()
    pdf_names = []
    local_coords = defaultdict(dict)
    for nsl in proc_hdfkeys[rank]:
        pbar.set_description("Rank {}: {}".format(rank, nsl))
        net, sta, loc = nsl.split('.')

        # note that ned contains a single station
        ned_rf = NetworkEventDataset(src_h5_event_file, network=net, station=sta, location=loc, root=rf_h5_root)
        ned_swp = NetworkEventDataset(src_h5_event_file, network=net, station=sta, location=loc, root=sw_h5_root)

        curr_output_file = os.path.join(tempdir, '{}.pdf'.format(nsl))

        results_rf = None
        results_swp = None
        with PdfPages(curr_output_file) as pdf:
            fig, (ax1, ax2) = plt.subplots(2, 1)
            fig.set_size_inches(paper_size_A4[1], paper_size_A4[0]) #landscape

            fig.suptitle(nsl, fontsize=16)
            ax1.set_title('Receiver Function')
            ax2.set_title('Surface-wave Polarization')

            results_rf = rf_station_orientations(ned_rf, ax=ax1)
            results_swp = swp_station_orientations(ned_swp, grv_dict, ax=ax2)

            plt.tight_layout()
            pdf.savefig(dpi=300, orientation='portrait')
            plt.close()
        # end with

        local_results_rf.update(results_rf)
        local_results_swp.update(results_swp)
        pdf_names.append(curr_output_file)

        local_coords.update(get_station_coords(ned_rf))

        pbar.update()
    # end for

    comm.barrier()

    global_results_rf = comm.gather(local_results_rf, root=0)
    global_results_swp = comm.gather(local_results_swp, root=0)
    pdf_names = comm.gather(pdf_names, root=0)
    global_coords = comm.gather(local_coords, root=0)

    if (rank == 0):
        flat_results_dict = defaultdict(lambda: defaultdict(dict))
        flat_coords_dict = defaultdict(dict)
        for d in global_results_rf: flat_results_dict['rf'].update(d)
        for d in global_results_swp: flat_results_dict['swp'].update(d)

        for d in global_coords: flat_coords_dict.update(d)

        json_fn = output_basename + '.json'
        with open(json_fn, 'w') as f:
            json.dump(flat_results_dict, f, indent=4)
        #end with

        # flatten list and merge pdfs
        pdf_fn = output_basename + '.pdf'
        pdf_names = [item for items in pdf_names for item in items]

        # plot summary
        logger.info('Plotting summary..')
        summary_ofn = os.path.join(tempdir, 'summary.pdf')
        plot_summary(flat_coords_dict, flat_results_dict, summary_ofn)
        pdf_names = [summary_ofn] + pdf_names

        pdf_merge(pdf_names, pdf_fn)

        rmtree(tempdir)

        logger.info("Finishing...")
        logger.info("bulk_station_orientations SUCCESS!")
    # end if
# end func

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
# end if
