import matplotlib.pyplot as plt
import numpy
from ordered_set import OrderedSet as set
from matplotlib.backends.backend_pdf import PdfPages
import cartopy.crs as ccrs
import numpy as np
from pyproj import Geod
from geographiclib.geodesic import Geodesic
from tqdm import tqdm
from cycler import cycler
from seismic.pick_harvester.ssst_utils import SSST_Result
from seismic.pick_harvester.cluster_utils import NestedGrid
from collections import defaultdict
P_CUTOFF = 5
S_CUTOFF = 10

def cover_page(pdf:PdfPages, text):
    page = plt.figure(figsize=(11.69,8.27))
    page.clf()
    page.text(0.5, 0.5, text, transform=page.transFigure, size=24, ha="center")
    pdf.savefig()
# end func

def plot_before_cluster(sr:SSST_Result, ng:NestedGrid, pdf:PdfPages, min_slope_ratio=5):
    cover_page(pdf, "Before Clustering")

    equality = sr.equality
    residual = sr.residual
    is_AUTO_arrival = sr.is_AUTO_arrival
    ecdists = sr.ecdists
    is_P = sr.is_P
    is_S = sr.is_S
    slope_ratio = sr.slope_ratio

    if(0):
        geod = Geod(a=180 / np.pi, f=0)

        elons = sr.elons
        elats = sr.elats
        slons = sr.slons
        slats = sr.slats

        # define 2D grid
        sx, ex, dx = -180, 180, ng.ig.dx
        sy, ey, dy = 90, -90, ng.ig.dy

        nx = np.int_(np.ceil((ex - sx) / dx)) + 1
        ny = np.int_(np.ceil((sy - ey) / dy)) + 1
        res = 1 / float(dx)  # per deg

        gx, gy = np.meshgrid(np.linspace(sx, ex, nx), np.linspace(sy, ey, ny), indexing='ij')

        gzp = np.zeros(gx.shape)
        gzs = np.zeros(gx.shape)

        # ===========================================================
        # Process P arrivals
        # ===========================================================
        imask = (equality) & (is_P) & (np.fabs(residual) < P_CUTOFF) & \
                (~is_AUTO_arrival | (is_AUTO_arrival & (slope_ratio > min_slope_ratio)))

        p_elons = elons[imask]
        p_elats = elats[imask]
        p_slons = slons[imask]
        p_slats = slats[imask]
        p_ecdists = ecdists[imask]
        for i in tqdm(np.arange(0, len(p_elons), 1), desc='P-arrivals: '):
            elon, elat, slon, slat = p_elons[i], p_elats[i], \
                                     p_slons[i], p_slats[i]

            npts = int(np.ceil(p_ecdists[i] * res))
            xy = np.vstack([ np.array([elon, elat]),
                             np.array(geod.npts(elon, elat, slon, slat, npts)),
                             np.array([slon, slat]) ])
            cxi = np.int_((xy[:, 0] - sx) / dx)
            cyi = np.int_((sy - xy[:, 1]) / dy)
            gzp[cxi, cyi] += 1
        # end for

        # ===========================================================
        # Process S arrivals
        # ===========================================================
        imask = (equality) & (is_S) & (np.fabs(residual) < S_CUTOFF) & \
                (~is_AUTO_arrival | (is_AUTO_arrival & (slope_ratio > min_slope_ratio)))

        s_elons = elons[imask]
        s_elats = elats[imask]
        s_slons = slons[imask]
        s_slats = slats[imask]
        s_ecdists = ecdists[imask]
        for i in tqdm(np.arange(0, len(s_elons), 1), desc='S-arrivals'):
            elon, elat, slon, slat = s_elons[i], s_elats[i], \
                                     s_slons[i], s_slats[i]

            npts = int(np.ceil(s_ecdists[i] * res))
            xy = np.vstack([ np.array([elon, elat]),
                             np.array(geod.npts(elon, elat, slon, slat, npts)),
                             np.array([slon, slat]) ])

            cxi = np.int_((xy[:, 0] - sx) / dx)
            cyi = np.int_((sy - xy[:, 1]) / dy)
            gzs[cxi, cyi] += 1
        # end for

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

        llcrnrlat = -54
        urcrnrlat = 0
        llcrnrlon = 100
        urcrnrlon = 180
        vmin = 20

        ax1 = fig.add_subplot(1, 2, 1, projection=ccrs.PlateCarree())
        ax1.set_extent([llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat], crs=ccrs.PlateCarree())
        ax1.coastlines('50m')

        cbi = ax1.pcolormesh(gx, gy, gzp, vmin=vmin, vmax=np.max(gzp), cmap=plt.get_cmap('gist_heat_r', 50))

        cbi.cmap.set_under('green')
        cbar = fig.colorbar(cbi, ax=cax1, format='%d', extend='min',
                            orientation='horizontal')
        cbar.set_label("Hitcount")
        ticks = list(cbar.get_ticks())
        cbar.set_ticks([vmin] + ticks)
        cbar.ax.tick_params(labelsize=6)
        ax1.set_title('P-arrivals')

        # Plot s-coverage
        ax2 = fig.add_subplot(1, 2, 2, projection=ccrs.PlateCarree())
        ax2.set_extent([llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat], crs=ccrs.PlateCarree())
        ax2.coastlines('50m')

        cbi = ax2.pcolormesh(gx, gy, gzs, vmin=vmin, vmax=np.max(gzs), cmap=plt.get_cmap('gist_heat_r', 50))

        cbi.cmap.set_under('green')
        cbar = fig.colorbar(cbi, ax=cax2, format='%d', extend='min',
                            orientation='horizontal')
        cbar.set_label("Hitcount")
        ticks = list(cbar.get_ticks())
        cbar.set_ticks([vmin] + ticks)
        cbar.ax.tick_params(labelsize=6)
        ax2.set_title('S-arrivals')

        fig.suptitle('Surface-projected Raypath Coverage in the Australasian Region', size=18)
        pdf.savefig(dpi=300)
        plt.close()
    # end if

    # ===========================================================
    # Gather indices of different phases
    # ===========================================================
    phase_indices_dict = {}
    for p in [b'P', b'Pg', b'Pb', b'Pn', b'S', b'Sg', b'Sb', b'Sn']:
        phase_indices_dict[p] = sr.phase == p
    # end for

    t = sr.arrival_ts - sr.eorigin_ts
    custom_cycler = (cycler(color=['r', 'g', 'b', 'm']))

    # ===========================================================
    # Plot distributions of residuals for all preexisting arrivals
    # to be exported
    # ===========================================================
    p_imask = (equality) & (is_P) & \
              (~is_AUTO_arrival) & (np.fabs(residual) < P_CUTOFF)
    s_imask = (equality) & (is_S) & \
              (~is_AUTO_arrival) & (np.fabs(residual) < S_CUTOFF)

    fig, axes = plt.subplots(1, 4)
    fig.set_size_inches(20, 5)

    _ = axes[0].hist(residual[p_imask], bins=20)
    _ = axes[1].hist(residual[s_imask], bins=20)

    axes[0].set_xlabel('Residual [s]'); axes[0].set_ylabel('Frequency')
    axes[1].set_xlabel('Residual [s]'); axes[1].set_ylabel('Frequency')
    axes[0].set_title('P-arrivals'); axes[1].set_title('S-arrivals')

    axes[0].text(0.7, 0.7, 'N: {}'.format(np.sum(p_imask)), transform=axes[0].transAxes)
    axes[0].text(0.7, 0.65, 'std: {:0.3f}'.format(np.std(residual[p_imask])), transform=axes[0].transAxes)
    axes[1].text(0.7, 0.7, 'N: {}'.format(np.sum(s_imask)), transform=axes[1].transAxes)
    axes[1].text(0.7, 0.65, 'std: {:0.3f}'.format(np.std(residual[s_imask])), transform=axes[1].transAxes)

    # distance-time scatter plots
    axes[2].set_prop_cycle(custom_cycler)
    axes[3].set_prop_cycle(custom_cycler)

    for k in phase_indices_dict.keys():
        x1 = y1 = x2 = y2 = None
        if (b'P' in k):
            x1 = ecdists[phase_indices_dict[k] & p_imask]
            y1 = t[phase_indices_dict[k] & p_imask]
            axes[2].scatter(x1, y1, s=0.1, label=k.decode(), rasterized=True)
        elif (b'S' in k):
            x2 = ecdists[phase_indices_dict[k] & s_imask]
            y2 = t[phase_indices_dict[k] & s_imask]
            axes[3].scatter(x2, y2, s=0.1, label=k.decode(), rasterized=True)
        # end if
    # end for
    lg0 = axes[2].legend()
    lg1 = axes[3].legend()

    for handle in lg0.legendHandles: handle.set_sizes([5.0])
    for handle in lg1.legendHandles: handle.set_sizes([5.0])
    axes[2].set_xlabel('Distance [°]')
    axes[2].set_ylabel('Time [s]')
    axes[3].set_xlabel('Distance [°]')
    axes[3].set_ylabel('Time [s]')

    fig.suptitle('All Preexisting Arrivals', fontsize=18)
    pdf.savefig(dpi=300)
    plt.close()

    # ===========================================================
    # Plot distributions of residuals for all automatic
    # arrivals to be exported
    # ===========================================================
    p_imask = (equality) & (is_P) & \
              (is_AUTO_arrival & (slope_ratio > min_slope_ratio)) & \
              (np.fabs(residual) < P_CUTOFF)
    s_imask = (equality) & (is_S) & \
              (is_AUTO_arrival & (slope_ratio > min_slope_ratio)) & \
              (np.fabs(residual) < S_CUTOFF)

    fig, axes = plt.subplots(1, 4)
    fig.set_size_inches(20, 5)

    _ = axes[0].hist(residual[p_imask], bins=20)
    _ = axes[1].hist(residual[s_imask], bins=20)

    axes[0].set_xlabel('Residual [s]'); axes[0].set_ylabel('Frequency')
    axes[1].set_xlabel('Residual [s]'); axes[1].set_ylabel('Frequency')
    axes[0].set_title('P-arrivals'); axes[1].set_title('S-arrivals')

    axes[0].text(0.7, 0.7, 'N: {}'.format(np.sum(p_imask)), transform=axes[0].transAxes)
    axes[0].text(0.7, 0.65, 'std: {:0.3f}'.format(np.std(residual[p_imask])), transform=axes[0].transAxes)
    axes[1].text(0.7, 0.7, 'N: {}'.format(np.sum(s_imask)), transform=axes[1].transAxes)
    axes[1].text(0.7, 0.65, 'std: {:0.3f}'.format(np.std(residual[s_imask])), transform=axes[1].transAxes)

    # distance-time scatter plots
    axes[2].set_prop_cycle(custom_cycler)
    axes[3].set_prop_cycle(custom_cycler)

    for k in phase_indices_dict.keys():
        x1 = y1 = x2 = y2 = None
        if (b'P' in k):
            x1 = ecdists[phase_indices_dict[k] & p_imask]
            y1 = t[phase_indices_dict[k] & p_imask]
            axes[2].scatter(x1, y1, s=0.1, label=k.decode(), rasterized=True)
        elif (b'S' in k):
            x2 = ecdists[phase_indices_dict[k] & s_imask]
            y2 = t[phase_indices_dict[k] & s_imask]
            axes[3].scatter(x2, y2, s=0.1, label=k.decode(), rasterized=True)
        # end if
    # end for
    lg0 = axes[2].legend()
    lg1 = axes[3].legend()

    for handle in lg0.legendHandles: handle.set_sizes([5.0])
    for handle in lg1.legendHandles: handle.set_sizes([5.0])
    axes[2].set_xlabel('Distance [°]')
    axes[2].set_ylabel('Time [s]')
    axes[3].set_xlabel('Distance [°]')
    axes[3].set_ylabel('Time [s]')

    fig.suptitle('All Automatic Arrivals (slope-ratio > {})'.format(min_slope_ratio), fontsize=18)
    pdf.savefig(dpi=300)
    plt.close()
# end func

def plot_after_cluster(p_clustered:np.ndarray, s_clustered:numpy.ndarray, ng: NestedGrid,
                       phases:str, pdf: PdfPages):
    cover_page(pdf, "After Clustering")

    if(0):
        geod = Geod(a=180 / np.pi, f=0)

        # define 2D grid
        sx, ex, dx = -180, 180, ng.ig.dx
        sy, ey, dy = 90, -90, ng.ig.dy

        nx = np.int_(np.ceil((ex - sx) / dx)) + 1
        ny = np.int_(np.ceil((sy - ey) / dy)) + 1
        res = 1 / float(dx)  # per deg

        gx, gy = np.meshgrid(np.linspace(sx, ex, nx), np.linspace(sy, ey, ny), indexing='ij')

        gzp = np.zeros(gx.shape)
        gzs = np.zeros(gx.shape)

        # ===========================================================
        # Process clustered P arrivals
        # ===========================================================
        p_ig_orig_dest = 0
        p_ig_outside = 0
        p_ig_intersecting = 0
        for i in tqdm(np.arange(0, len(p_clustered), 1), desc='P-arrivals'):
            elon, elat, slon, slat, ecdist, edepth_km = p_clustered['elon'][i], p_clustered['elat'][i], \
                                                        p_clustered['slon'][i], p_clustered['slat'][i], \
                                                        p_clustered['ecdist'][i], p_clustered['edepth_km'][i]

            npts = int(np.ceil(ecdist * res))
            xy = np.vstack([np.array([elon, elat]),
                            np.array(geod.npts(elon, elat, slon, slat, npts)),
                            np.array([slon, slat])])
            cxi = np.int_((xy[:, 0] - sx) / dx)
            cyi = np.int_((sy - xy[:, 1]) / dy)
            gzp[cxi, cyi] += 1

            # gather statistics on ray-geometry
            if(elon < 0): elon += 360
            if(slon < 0): slon += 360
            contains_orig = ng.ig._is_within_grid(elon, 90 - elat, edepth_km)
            contains_dest = ng.ig._is_within_grid(slon, 90 - slat, 0)

            if(contains_orig or contains_dest):
                p_ig_orig_dest += 1
                continue
            # end if

            intersects = False
            for cxy in xy:
                lon = cxy[0]
                colat = 90 - cxy[1]
                if(lon < 0): lon += 360

                if(ng.ig._is_within_grid(lon, colat, 0)):
                    intersects = True
                    break
                # end if
            # end for

            if(intersects):
                p_ig_intersecting += 1
            else:
                p_ig_outside += 1
            # end if
        # end for

        # ===========================================================
        # Process clustered S arrivals
        # ===========================================================
        s_ig_orig_dest = 0
        s_ig_outside = 0
        s_ig_intersecting = 0
        for i in tqdm(np.arange(0, len(s_clustered), 1), desc='S-arrivals'):
            elon, elat, slon, slat, ecdist, edepth_km = s_clustered['elon'][i], s_clustered['elat'][i], \
                                                        s_clustered['slon'][i], s_clustered['slat'][i], \
                                                        s_clustered['ecdist'][i], s_clustered['edepth_km'][i]

            npts = int(np.ceil(ecdist * res))
            xy = np.vstack([np.array([elon, elat]),
                            np.array(geod.npts(elon, elat, slon, slat, npts)),
                            np.array([slon, slat])])
            cxi = np.int_((xy[:, 0] - sx) / dx)
            cyi = np.int_((sy - xy[:, 1]) / dy)
            gzs[cxi, cyi] += 1

            # gather statistics on ray-geometry
            if(elon < 0): elon += 360
            if(slon < 0): slon += 360
            contains_orig = ng.ig._is_within_grid(elon, 90 - elat, edepth_km)
            contains_dest = ng.ig._is_within_grid(slon, 90 - slat, 0)

            if(contains_orig or contains_dest):
                s_ig_orig_dest += 1
                continue
            # end if

            intersects = False
            for cxy in xy:
                lon = cxy[0]
                colat = 90 - cxy[1]
                if(lon < 0): lon += 360

                if(ng.ig._is_within_grid(lon, colat, 0)):
                    intersects = True
                    break
                # end if
            # end for

            if(intersects):
                s_ig_intersecting += 1
            else:
                s_ig_outside += 1
            # end if
        # end for

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

        llcrnrlat = -54
        urcrnrlat = 0
        llcrnrlon = 100
        urcrnrlon = 180
        vmin = 1

        ax1 = fig.add_subplot(1, 2, 1, projection=ccrs.PlateCarree())
        ax1.set_extent([llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat], crs=ccrs.PlateCarree())
        ax1.coastlines('50m')

        cbi = ax1.pcolormesh(gx, gy, gzp, vmin=vmin, vmax=np.max(gzp), cmap=plt.get_cmap('gist_heat_r', 50))

        cbi.cmap.set_under('green')
        cbar = fig.colorbar(cbi, ax=cax1, format='%d', extend='min',
                            orientation='horizontal')
        cbar.set_label("Hitcount")
        ticks = list(cbar.get_ticks())
        cbar.set_ticks([vmin] + ticks)
        cbar.ax.tick_params(labelsize=6)
        ax1.set_title('P-arrivals')

        ax1.text(0, -0.05, '1. Rays with event origin or station within inner grid: {}'.format(p_ig_orig_dest),
                 fontsize=10, transform=ax1.transAxes)
        ax1.text(0, -0.1, '2. Rays with both event origin and station outside inner grid: {}'.format(p_ig_outside),
                 fontsize=10, transform=ax1.transAxes)
        ax1.text(0, -0.15, '3. Rays in 2. intersecting inner grid: {}'.format(p_ig_intersecting),
                 fontsize=10, transform=ax1.transAxes)

        # Plot s-coverage
        ax2 = fig.add_subplot(1, 2, 2, projection=ccrs.PlateCarree())
        ax2.set_extent([llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat], crs=ccrs.PlateCarree())
        ax2.coastlines('50m')

        cbi = ax2.pcolormesh(gx, gy, gzs, vmin=vmin, vmax=np.max(gzs), cmap=plt.get_cmap('gist_heat_r', 50))

        cbi.cmap.set_under('green')
        cbar = fig.colorbar(cbi, ax=cax2, format='%d', extend='min',
                            orientation='horizontal')
        cbar.set_label("Hitcount")
        ticks = list(cbar.get_ticks())
        cbar.set_ticks([vmin] + ticks)
        cbar.ax.tick_params(labelsize=6)
        ax2.set_title('S-arrivals')

        ax2.text(0, -0.05, '1. Rays with event origin or station within inner grid: {}'.format(s_ig_orig_dest),
                 fontsize=10, transform=ax2.transAxes)
        ax2.text(0, -0.1, '2. Rays with both event origin and station outside inner grid: {}'.format(s_ig_outside),
                 fontsize=10, transform=ax2.transAxes)
        ax2.text(0, -0.15, '3. Rays in 2. intersecting inner grid: {}'.format(s_ig_intersecting),
                 fontsize=10, transform=ax2.transAxes)

        fig.suptitle('Surface-projected Clustered Raypath Coverage in the Australasian Region', size=18)
        pdf.savefig(dpi=300)
    # end if

    # ===========================================================
    # Plot distributions of residuals for clustered arrivals
    # ===========================================================
    fig, axes = plt.subplots(1, 2)
    fig.set_size_inches(10, 5)

    _ = axes[0].hist(p_clustered['residual'], bins=20)
    _ = axes[1].hist(s_clustered['residual'], bins=20)

    axes[0].set_xlabel('Residual [s]'); axes[0].set_ylabel('Frequency')
    axes[1].set_xlabel('Residual [s]'); axes[1].set_ylabel('Frequency')
    axes[0].set_title('P-arrivals'); axes[1].set_title('S-arrivals')

    axes[0].text(0.7, 0.7, 'N: {}'.format(len(p_clustered)), transform=axes[0].transAxes)
    axes[0].text(0.7, 0.65, 'std: {:0.3f}'.format(np.std(p_clustered['residual'])), transform=axes[0].transAxes)
    axes[1].text(0.7, 0.7, 'N: {}'.format(len(s_clustered)), transform=axes[1].transAxes)
    axes[1].text(0.7, 0.65, 'std: {:0.3f}'.format(np.std(s_clustered['residual'])), transform=axes[1].transAxes)

    p_phases = [phase for phase in phases.split(' ') if 'P' in phase]
    s_phases = [phase for phase in phases.split(' ') if 'S' in phase]

    phc = {}
    for p in p_phases: phc[p] = np.sum(p.encode() == p_clustered['phase'].astype('S10'))
    axes[0].set_title('Phase counts P: {}\n'.format(phc), fontdict={'fontsize': 8},
                      pad=30)

    phc = {}
    for p in s_phases: phc[p] = np.sum(p.encode() == s_clustered['phase'].astype('S10'))
    axes[1].set_title('Phase counts S: {}\n'.format(phc), fontdict={'fontsize': 8},
                      pad=30)

    fig.suptitle('Distributions of residuals for clustered arrivals', fontsize=18)
    plt.tight_layout()
    pdf.savefig(dpi=300)

    # ===========================================================
    # Plot distribution of clustered P arrivals
    # ===========================================================
    fig, axes = plt.subplots(1, 2)
    fig.set_size_inches(10, 5)

    count_dict_regional = defaultdict(int)
    count_dict_global = defaultdict(int)
    for i, j in zip(p_clustered['source_block'], p_clustered['station_block']):
        if (ng.is_inner_block(i) or ng.is_inner_block(j)):
            count_dict_regional[(i, j)] += 1
        else:
            count_dict_global[(i, j)] += 1
        # end if
    # end for

    counts_regional = []
    counts_global = []
    for k, c in count_dict_regional.items():
        counts_regional.append(c)
    # end for
    for k, c in count_dict_global.items():
        counts_global.append(c)
    # end for

    counts_regional = np.array(counts_regional)
    counts_global = np.array(counts_global)

    bins = np.histogram_bin_edges(counts_regional)
    logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
    axes[0].hist(counts_regional, bins=logbins)
    axes[0].set_xscale('log')
    axes[0].set_xlim(1, np.max(logbins))
    axes[0].set_title('Regional Grid')
    axes[0].set_xlabel('Number of rays in cluster')
    axes[0].set_ylabel('Frequency')
    axes[0].text(0.7, 0.7, 'N: {}\nmin: {}\nmax: {}\nmean: {:.3f}\nstd: {:.3f}\n'.format(
        len(counts_regional), np.min(counts_regional),
        np.max(counts_regional), np.mean(counts_regional),
        np.std(counts_regional)), transform=axes[0].transAxes)

    bins = np.histogram_bin_edges(counts_global)
    logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
    axes[1].hist(counts_global, bins=logbins)
    axes[1].set_xscale('log')
    axes[1].set_xlim(1, np.max(logbins))
    axes[1].set_title('Global Grid')
    axes[1].set_xlabel('Number of rays in cluster')
    axes[1].text(0.7, 0.7, 'N: {}\nmin: {}\nmax: {}\nmean: {:.3f}\nstd: {:.3f}\n'.format(
        len(counts_global), np.min(counts_global),
        np.max(counts_global), np.mean(counts_global),
        np.std(counts_global)), transform=axes[1].transAxes)

    fig.suptitle('Clustered P-arrivals')
    pdf.savefig(dpi=300)

    # ===========================================================
    # Plot distribution of clustered S arrivals
    # ===========================================================
    fig, axes = plt.subplots(1, 2)
    fig.set_size_inches(10, 5)

    count_dict_regional = defaultdict(int)
    count_dict_global = defaultdict(int)
    for i, j in zip(s_clustered['source_block'], s_clustered['station_block']):
        if (ng.is_inner_block(i) or ng.is_inner_block(j)):
            count_dict_regional[(i, j)] += 1
        else:
            count_dict_global[(i, j)] += 1
        # end if
    # end for

    counts_regional = []
    counts_global = []
    for k, c in count_dict_regional.items():
        counts_regional.append(c)
    # end for
    for k, c in count_dict_global.items():
        counts_global.append(c)
    # end for

    counts_regional = np.array(counts_regional)
    counts_global = np.array(counts_global)

    bins = np.histogram_bin_edges(counts_regional)
    logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
    axes[0].hist(counts_regional, bins=logbins)
    axes[0].set_xscale('log')
    axes[0].set_xlim(1, np.max(logbins))
    axes[0].set_title('Regional Grid')
    axes[0].set_xlabel('Number of rays in cluster')
    axes[0].set_ylabel('Frequency')
    axes[0].text(0.7, 0.7, 'N: {}\nmin: {}\nmax: {}\nmean: {:.3f}\nstd: {:.3f}\n'.format(
        len(counts_regional), np.min(counts_regional),
        np.max(counts_regional), np.mean(counts_regional),
        np.std(counts_regional)), transform=axes[0].transAxes)

    bins = np.histogram_bin_edges(counts_global)
    logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
    axes[1].hist(counts_global, bins=logbins)
    axes[1].set_xscale('log')
    axes[1].set_xlim(1, np.max(logbins))
    axes[1].set_title('Global Grid')
    axes[1].set_xlabel('Number of rays in cluster')
    axes[1].text(0.7, 0.7, 'N: {}\nmin: {}\nmax: {}\nmean: {:.3f}\nstd: {:.3f}\n'.format(
        len(counts_global), np.min(counts_global),
        np.max(counts_global), np.mean(counts_global),
        np.std(counts_global)), transform=axes[1].transAxes)

    fig.suptitle('Clustered S-arrivals')
    pdf.savefig(dpi=300)

    plt.close()
# end func
