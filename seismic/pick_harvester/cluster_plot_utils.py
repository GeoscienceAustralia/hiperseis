import matplotlib.pyplot as plt
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
    geod = Geod(a=180 / np.pi, f=0)

    equality = sr.equality
    residual = sr.residual
    is_AUTO_arrival = sr.is_AUTO_arrival

    elons = sr.elons
    elats = sr.elats
    slons = sr.slons
    slats = sr.slats
    ecdists = sr.ecdists
    is_P = sr.is_P
    is_S = sr.is_S
    slope_ratio = sr.slope_ratio

    # define 2D grid
    sx, ex, dx = -180, 180, ng.ig.dx
    sy, ey, dy = 90, -90, ng.ig.dy

    nx = np.int_(np.ceil((ex - sx) / dx)) + 1
    ny = np.int_(np.ceil((sy - ey) / dy)) + 1
    res = int(1 / float(dx))  # per deg

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
    for i in tqdm(np.arange(0, len(p_elons), 100), desc='P-arrivals: '):
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
    for i in tqdm(np.arange(0, len(s_elons), 10), desc='S-arrivals'):
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

    pdf.close()
# end func
