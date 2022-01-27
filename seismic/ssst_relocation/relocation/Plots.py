# -*- coding: utf-8 -*-
"""
Created on Thu Dec 23 14:36:56 2021

@author: U37509
"""

import argparse, os, time
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER

def compute_geodesics(lon0, colat0, lon1, colat1):
    # cot(colat) = u*cos(lon - v)
    #            = u*cos(lon)*cos(v) + u*sin(lon)*sin(v)
    #            = a*cos(lon) + b*sin(lon)
    # a = u*cos(v), b = u*sin(v)
    # u = sqrt(a^2 + b^2)
    # v = arctan(b/a)
    
    # a*cos(0) + b*sin(0) = 1/tan(colat0)
    # a*cos(L) + b*sin(L) = 1/tan(colat1)
    
    # a*cos(lon0) + b*sin(lon0) = 1/tan(colat0)
    # a*cos(lon1) + b*sin(lon1) = 1/tan(colat1)
    # [[cos(lon0), sin(lon0)], [cos(lon1), sin(lon1)]] * [a, b] = \
    #     [1/tan(colat0), 1/tan(colat1)]
    # M*[a, b] = [c[0], c[1]]
    # [a, b] = [[cos(lon0), sin(lon0)], 
    #           [cos(lon1), sin(lon1)]]^{-1} * [1/tan(colat0), 1/tan(colat1)]
    #        = inv(M)*c
    
    # M = [[cos(lon0), sin(lon0)], 
    #      [cos(lon1), sin(lon1)]]
    # det(M) = cos(lon0)*sin(lon1) - cos(lon1)*sin(lon0)
    # N = M^{-1}
    #   = [[sin(lon1), -sin(lon0)], 
    #      [-cos(lon1), cos(lon0)]]/(cos(lon0)*sin(lon1) - cos(lon1)*sin(lon0))
    
    degrad = np.pi/180
    lon0 = lon0*degrad
    colat0 = colat0*degrad
    lon1 = lon1*degrad
    colat1 = colat1*degrad
    
    a = (np.sin(lon1)/np.tan(colat0) - np.sin(lon0)/np.tan(colat1))/ \
            (np.cos(lon0)*np.sin(lon1) - np.cos(lon1)*np.sin(lon0))
    b = (np.cos(lon0)/np.tan(colat1) - np.cos(lon1)/np.tan(colat0))/ \
            (np.cos(lon0)*np.sin(lon1) - np.cos(lon1)*np.sin(lon0))
            
    u = np.sqrt(a**2 + b**2)
    v = np.arctan(b/a)
       
    val1 = (np.pi - np.arctan(1/(u*np.cos(lon0 - v)))) % np.pi
    val2 = (np.pi - np.arctan(1/(u*np.cos(lon0 - v - np.pi)))) % np.pi
    
    ind = np.abs(np.pi - colat0 - val1) > np.abs(np.pi - colat0 - val2)
    v[ind] = v[ind] + np.pi
    
    return u, v/degrad
#end func
    
def compute_rays(picks):
    elon = picks['elon']
    ecolat = 90.0 - picks['elat']
    slon = picks['slon']
    scolat = 90.0 - picks['slat']
    
    u, v = compute_geodesics(elon, ecolat, slon, scolat)
    
    x0 = np.min([elon, slon], axis=0)
    x1 = np.max([elon, slon], axis=0)
    
    ind = (x1 - x0) > 180.0
    start = np.zeros_like(x0)
    end = np.zeros_like(x0)
    start[~ind] = x0[~ind]
    end[~ind] = x1[~ind]
    start[ind] = x1[ind]
    end[ind] = x0[ind]
    
    return u, v, start, end
#end func
    
def filter_rays(rays, lonmin, lonmax, latmin, latmax):
    degrad = np.pi/180.0
    u, v, start, end = rays
    
    lons = (np.linspace(0.0, (end - start) % 360.0).T + \
            np.expand_dims(start, 1)) % 360.0
    lons[lons > 180.0] = lons[lons > 180.0] - 360.0
    lons = lons*degrad
    v_temp = np.expand_dims(v, 1)*degrad
    
    # colat = arccot(u*cos(lon - v))
    #       = arctan(1/(u*cos(lon - v)))
    colats = np.arctan(1/(np.expand_dims(u, 1)*np.cos(lons - v_temp)))
    lats = (np.pi - colats) % np.pi - np.pi/2
    
    lons = lons/degrad
    lats = lats/degrad
    
    if lonmax == 180.0:
        ind1 = (lons >= lonmin)
    else:
        ind1 = np.logical_and(lons >= lonmin, 
                              (lons + 180.0) % 360.0 <= \
                              (lonmax + 180.0) % 360.0)
    #end if
    ind2 = np.logical_and(lats >= latmin, lats <= latmax)
    ind = np.logical_and(ind1, ind2)
    rows = np.argwhere(np.any(ind, axis=1)).T[0]
    
    return u[rows], v[rows], start[rows], end[rows]
#end func

def filter_hypocentres(hypocentres, lonmin=-180, lonmax=180, latmin=-90, 
                       latmax=90, depthmin=0, depthmax=6400000, timemin=0, 
                       timemax=1e10):
    filt = list()
    filt.append(hypocentres['elon'] >= lonmin)
    filt.append(hypocentres['elon'] <= lonmax)
    filt.append(hypocentres['elat'] >= latmin)
    filt.append(hypocentres['elat'] <= latmax)
    filt.append(hypocentres['edepth'] >= depthmin)
    filt.append(hypocentres['edepth'] <= depthmax)
    filt.append(hypocentres['origin_time'] >= timemin)
    filt.append(hypocentres['origin_time'] <= timemax)
    filt = np.all(filt, axis=0)
    return hypocentres[filt]
#end func

def hypocentres_from_picks(picks):
    lon = picks['elon']
    lat = picks['elat']
    depth = picks['edepth']
    ot = picks['origin_time']
    ind = np.unique(picks['event_id'], return_index=True)[1]
    temp = [tuple(row) for row in np.array([lon[ind], lat[ind], depth[ind], 
                                            ot[ind]]).T]
    hypocentres = np.array(temp, dtype=[('elon', float), ('elat', float), 
                                        ('edepth', float), 
                                        ('origin_time', float)])
    return hypocentres
#end func

def epicentre_scatterplot(hypocentres, title_prefix, output_path, show=False):
    plt.figure(figsize=(10, 10))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.stock_img()
    ax.coastlines()
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, 
                      color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_left = True
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    plt.scatter(hypocentres['elon'], hypocentres['elat'], s=0.5, c='b')
    plt.xlim(np.floor(np.min(hypocentres['elon'])), 
             np.ceil(np.max(hypocentres['elon'])))
    plt.ylim(np.floor(np.min(hypocentres['elat'])),
             np.ceil(np.max(hypocentres['elat'])))
    plt.title(str(title_prefix + ' epicentres'))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.savefig(os.path.join(output_path, 
                             str(title_prefix + '.epicentres.png')))
    if show:
        plt.show()
    #end if
#end func
    
def residual_histogram(picks, title_prefix, output_path, tcor_available=False,
                       show=False):
    ind = ~np.isin(picks['phase'], ['Px', 'Sx', 'X'])
    plt.figure()
    if tcor_available:
        plt.hist(picks['residual'][ind] - picks['tcor'][ind], bins=100, 
                 range=(-10, 10))
    else:
        plt.hist(picks['residual'][ind], bins=100, range=(-10, 10))
    #end if
    plt.title(str(title_prefix + ' residual histogram'))
    plt.xlabel('Residual (seconds)')
    plt.ylabel('Frequency')
    plt.savefig(os.path.join(output_path, 
                             str(title_prefix + '.residuals.png')))
    if show:
        plt.show()
    #end if
#end func
    
def raypath_plot(rays, lonmin, lonmax, latmin, latmax, ray_interval, bin_width, 
                 bin_min, title_prefix, output_path, show=False):
    degrad = np.pi/180.0
    u, v, start, end = rays
    u = u[::ray_interval]
    v = v[::ray_interval]
    start = start[::ray_interval]
    end = end[::ray_interval]
    
    if lonmax - lonmin == 360.0:
        lons = lonmin + np.arange(0.0, (bin_width + lonmax - lonmin), 
                                  bin_width)
    else:
        lons = lonmin + np.arange(0.0, (bin_width + lonmax - lonmin) % 360.0, 
                                  bin_width)
    #end if
    lat_bins = np.arange(latmin, latmax + bin_width, bin_width)
    
    lons = lons*degrad
    v_temp = v*degrad
    
    # colat = arccot(u*cos(lon - v))
    #       = arctan(1/(u*cos(lon - v)))
    colats = np.arctan(1/(u*np.cos(np.expand_dims(lons, 1) - v_temp)))
    lats = (np.pi - colats) % np.pi - np.pi/2
    
    lons = lons/degrad
    lats = lats/degrad
    
    end = start + (end - start) % 360.0
    
    ind1 = start > np.expand_dims(lons, 1)
    ind2 = end < np.expand_dims(lons, 1)
    ind = np.logical_or(ind1, ind2)
    lats[ind] = np.nan
    
    hist = np.apply_along_axis(lambda a: np.histogram(a, bins=lat_bins)[0], 
                               1, lats).astype(float)
    hist[hist < bin_min] = np.nan
    
    plt.figure(figsize=(10, 10))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.stock_img()
    ax.coastlines()
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, 
                      color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_left = True
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    plt.plot(lons, lats, c='m', lw=0.01)
    if lonmax - lonmin == 360.0:
        plt.xlim([lonmin, lonmax])
    else:
        plt.xlim([lonmin, lonmin + (lonmax - lonmin) % 360.0])
    #end if
    plt.ylim([latmin, latmax])
    plt.title(str(title_prefix + ' rays'))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.savefig(os.path.join(output_path, str(title_prefix + '.rays.png')))
    if show:
        plt.show()
    #end if
    
    plt.figure(figsize=(10, 10))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.stock_img()
    ax.coastlines()
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, 
                      color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_left = True
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    im = plt.pcolormesh(lons, lat_bins[:-1], hist.T, vmin=1, 
                        vmax=np.nanmax(hist), cmap='hot_r')
    cbar = plt.colorbar(im, orientation='horizontal', fraction=0.046, pad=0.04,
                        extend='min')
    cbar.set_label('Hitcount')
    if lonmax - lonmin == 360.0:
        plt.xlim([lonmin, lonmax])
    else:
        plt.xlim([lonmin, lonmin + (lonmax - lonmin) % 360.0])
    #end if
    plt.ylim([latmin, latmax])
    plt.title(str(title_prefix + ' rays heatplot'))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.savefig(os.path.join(output_path, 
                             str(title_prefix + '.rays_heatplot.png')))
    if show:
        plt.show()
    #end if
#end func
    
def read_input_txt(filename, tcor_available):
    with open(filename, 'r') as file:
        rows = file.read().splitlines()
    #end with
    rows = [row.split() for row in rows[1:]]
    for row in rows:
        row[0] = row[0].split('/')[-1]
    #end for
    
    dtype=[('event_id', int), ('origin_time', float), ('mag', float), 
           ('elon', float), ('elat', float), ('edepth', float), ('net', 'U2'), 
           ('stat', 'U5'), ('cha', 'U3'), ('arrival_time', float), 
           ('phase', 'U8'), ('slon', float), ('slat', float), ('az', float), 
           ('baz', float), ('dist', float), ('residual', float), ('a', float), 
           ('b', float), ('c', float), ('d', float), ('e', float), 
           ('f', float)]
    if tcor_available:
        dtype.append(('tcor', float))
    #end if
    
    arr = np.array([tuple(row) for row in rows[1:]], dtype=dtype)
    return arr
#end func    

def process():
    parser = argparse.ArgumentParser(description='Plots')
    
    parser.add_argument("--files", type=str, nargs='+', required=True)
    parser.add_argument("--output_path", type=str, default='.')
    parser.add_argument("--show", type=bool, default=False)
    parser.add_argument("--lonmin", type=float, default=-180.0)
    parser.add_argument("--lonmax", type=float, default=180.0)
    parser.add_argument("--latmin", type=float, default=-90.0)
    parser.add_argument("--latmax", type=float, default=90.0)
    parser.add_argument("--depthmin", type=float, default=0.0)
    parser.add_argument("--depthmax", type=float, default=6400000.0)
    parser.add_argument("--timemin", type=float, default=0.0)
    parser.add_argument("--timemax", type=float, default=1e10)
    parser.add_argument("--plot_epicentres", type=bool, default=False)
    parser.add_argument("--plot_residuals", type=bool, default=False)
    parser.add_argument("--plot_raypaths", type=bool, default=False)
    parser.add_argument("--tcor_available", type=bool, default=False)
    parser.add_argument("--ray_interval", type=int, default=1)
    parser.add_argument("--ray_heatplot_bin_width", type=float, default=0.1)
    parser.add_argument("--ray_heatplot_bin_min", type=int, default=1)
    
    """
    import sys
    sys.argv = ['/g/data/ha3/la8536/SSST/relocation/Plots.py',
                "--files", 
                '/g/data/ha3/la8536/SSST/output_events/ensemble.p_copy.txt',
                "--output_path", '/g/data/ha3/la8536/SSST/output_events/',
                "--lonmin", '110', "--lonmax", '155', "--latmin", '-45', 
                "--latmax", '-10', "--plot_residuals", 'True', 
                "--plot_epicentres", 'True', "--show", 'True']
    """
    
    args = parser.parse_args()
    files = args.files
    output_path = args.output_path
    show = args.show
    plot_epicentres = args.plot_epicentres
    plot_residuals = args.plot_residuals
    plot_raypaths = args.plot_raypaths
    tcor_available = args.tcor_available
    ray_interval = args.ray_interval
    ray_heatplot_bin_width = args.ray_heatplot_bin_width
    ray_heatplot_bin_min = args.ray_heatplot_bin_min
    lonmin = args.lonmin
    lonmax = args.lonmax
    latmin = args.latmin
    latmax = args.latmax
    depthmin = args.depthmin
    depthmax = args.depthmax
    timemin = args.timemin
    timemax = args.timemax
    
    for file in files:
        t0 = time.time()
        print('Loading file, time = 0.0')
        #picks = np.load(file)
        picks = read_input_txt(file, tcor_available)
        title_prefix = os.path.basename(file)[:-4]
        
        if plot_epicentres:
            print('Finding hypocentres, time =', time.time() - t0)
            hypocentres = hypocentres_from_picks(picks)
            print('Filtering hypocentres, time =', time.time() - t0)
            hypocentres = filter_hypocentres(hypocentres, lonmin=lonmin, 
                                             lonmax=lonmax, latmin=latmin, 
                                             latmax=latmax, depthmin=depthmin, 
                                             depthmax=depthmax, 
                                             timemin=timemin, timemax=timemax)
            print('Producing epicentre scatterplot, time =', time.time() - t0)
            epicentre_scatterplot(hypocentres, title_prefix, output_path, 
                                  show=show)
        #end if
        
        if plot_raypaths:
            print('Computing', len(picks), 'geodesics, time =', 
                  time.time() - t0)
            rays = compute_rays(picks)
            print('Filtering geodesics, time =', time.time() - t0)
            rays = filter_rays(rays, lonmin, lonmax, latmin, latmax)
            print('Producing geodesic plot, time =', time.time() - t0)
            raypath_plot(rays, lonmin, lonmax, latmin, latmax, ray_interval,
                         ray_heatplot_bin_width, ray_heatplot_bin_min,
                         title_prefix, output_path, show=show)
        #end if
        
        if plot_residuals:
            print('Producing residual histogram, time =', time.time() - t0)
            residual_histogram(picks, title_prefix, output_path, 
                               tcor_available=tcor_available, show=show)
        #end if
        
        print('Done, time =', time.time() - t0)
    #end for
#end func

if __name__ == '__main__':
    process()
#end if
