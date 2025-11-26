#!/usr/bin/env python

import matplotlib as mpl
mpl.use('agg')
from matplotlib.colors import LightSource
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib import cm, colors

from netCDF4 import Dataset as NetCDFFile
import numpy as np
import click
import os
import h5py
from obspyh5 import iterh5
from seismic.receiver_fn.rf_ccp_util import Gravity
from adjustText import adjust_text

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

# below is work around of NCI python Basemap configuration problem, use line above for normal setup
# from mympl_toolkits.basemap import Basemap

def gmtColormap(fileName):
    """ Borrowed from AndrewStraw, scipy-cookbook
      this subroutine converts GMT cpt file to matplotlib palette """

    import colorsys

    try:
      f = open(fileName)
    except:
      print("file ",fileName, "not found")
      return None

    lines = f.readlines()
    f.close()

    x = []
    r = []
    g = []
    b = []
    colorModel = "RGB"
    for l in lines:
      ls = l.split()
      if l[0] == "#":
         if ls[-1] == "HSV":
             colorModel = "HSV"
             continue
         else:
             continue
      if ls[0] == "B" or ls[0] == "F" or ls[0] == "N":
         pass
      else:
          x.append(float(ls[0]))
          r.append(float(ls[1]))
          g.append(float(ls[2]))
          b.append(float(ls[3]))
          xtemp = float(ls[4])
          rtemp = float(ls[5])
          gtemp = float(ls[6])
          btemp = float(ls[7])

    x.append(xtemp)
    r.append(rtemp)
    g.append(gtemp)
    b.append(btemp)

    nTable = len(r)
    x = np.array( x )
    r = np.array( r )
    g = np.array( g )
    b = np.array( b )

    if colorModel == "HSV":
     for i in range(r.shape[0]):
         rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
         r[i] = rr ; g[i] = gg ; b[i] = bb
    if colorModel == "HSV":
     for i in range(r.shape[0]):
         rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
         r[i] = rr ; g[i] = gg ; b[i] = bb
    if colorModel == "RGB":
      r = r/255.
      g = g/255.
      b = b/255.


    xNorm = (x - x[0])/(x[-1] - x[0])

    red = []
    blue = []
    green = []
    for i in range(len(x)):
      red.append([xNorm[i],r[i],r[i]])
      green.append([xNorm[i],g[i],g[i]])
      blue.append([xNorm[i],b[i],b[i]])
    colorDict = {"red":red, "green":green, "blue":blue}

    return (x, colors.LinearSegmentedColormap('my_colormap',colorDict,255))
# end func

def rf_get_coords(rf_fn):
    coords = []
    names = []

    print('Reading RF station coordinates from {}..'.format(rf_fn))
    hf = h5py.File(rf_fn, 'r')

    for k1 in hf.keys():
        for k2 in hf[k1].keys():
            for k3 in hf[k1][k2].keys():
                group = '{}/{}/{}'.format(k1,k2,k3)
                for t in iterh5(rf_fn, group=group, headonly=True, mode='r'):
                    name = '.'.join(k2.split('.')[:-1])
                    if(name not in names):
                        names.append(name)
                        coords.append([t.stats.station_longitude, t.stats.station_latitude])
                    # end if
                    break
                # end func
                break
            # end for
        # end for
    # end for
    coords = np.array(coords)
    return names, coords
# end func

def plot_topo(coords, topo_grid, cpt_file, show_colorbar=False, alpha=1):

    fig=plt.figure(figsize=(11.69,8.27))
    plt.tick_params(labelsize=8)
    
    ax = fig.add_subplot(111)
    lon_min=min(coords[:,0])-1.
    lon_max=max(coords[:,0])+1.

    lat_1=min(coords[:,1])-1.
    lat_2=min(coords[:,1])
    lat_min=min(coords[:,1])-1.
    lat_max=max(coords[:,1])+1.
    
    lat_0=(lat_max+lat_min)/2.
    lon_0=(lon_max+lon_min)/2.

    m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,\
                lon_0=lon_0, lat_0=lat_0, \
                llcrnrlon=lon_min,llcrnrlat=lat_min, \
                urcrnrlon=lon_max,urcrnrlat=lat_max,\
                rsphere=6371200.,resolution='i')
    
    try:
        m.drawcoastlines()
    except:
        pass
    # end try

    #m.drawstates()
    #m.drawcountries()
    m.drawparallels(np.arange(-90.,90.,1.), labels=[1,0,0,0],fontsize=8, dashes=[2, 2], color='0.5', linewidth=0.75)
    m.drawmeridians(np.arange(0.,360.,1.), labels=[0,0,0,1], fontsize=8, dashes=[2, 2], color='0.5', linewidth=0.75)

    
#   ne below can draw scale but must be adjusted
#   m.drawmapscale(lon_min, lat_max, lon_max, lat_min, 400, fontsize = 16, barstyle='fancy', zorder=100)

    # Loading topography
    print('Reading topography grid {}..'.format(topo_grid))
    nc = NetCDFFile(topo_grid)
    
    zscale =20. #gray
    zscale =50. #colour
    data = nc.variables['elevation'][:] / zscale
    lons = nc.variables['lon'][:]
    lats = nc.variables['lat'][:]

    # transform to metres
    nx = int((m.xmax-m.xmin)/500.)+1
    ny = int((m.ymax-m.ymin)/500.)+1

    topodat = m.transform_scalar(data,lons,lats,nx,ny)

    zvals,cmap = gmtColormap(cpt_file)

    # make shading

    ls = LightSource(azdeg = 180, altdeg = 45)
    norm = colors.Normalize(vmin=-8000/zscale, vmax=5000/zscale)#myb
    rgb = ls.shade(topodat, cmap=cmap, norm=norm)
    im = m.imshow(rgb, alpha=alpha)
    
    if(show_colorbar): cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap))

    return fig, m
# end func

def plot_grav(coords, grav_grid, cpt_file, resolution=1, show_colorbar=False, alpha=1):
    DEG2KM = 111.
    dlonlat = resolution/DEG2KM

    # Loading gravity grid
    gravity = None
    try:
        print('Reading gravity grid {}..'.format(grav_grid))
        gravity = Gravity(grav_grid)
    except Exception as e:
        print(str(e))
        assert 0, 'Failed to load gravity grid. Aborting..'
    # end try

    # initialize plot
    fig=plt.figure(figsize=(11.69,8.27))
    plt.tick_params(labelsize=8)

    ax = fig.add_subplot(111)
    lon_min=min(coords[:,0])-1
    lon_max=max(coords[:,0])+1

    lat_1=min(coords[:,1])-1
    lat_2=min(coords[:,1])
    lat_min=min(coords[:,1])-1
    lat_max=max(coords[:,1])+1

    lat_0=(lat_max+lat_min)/2.
    lon_0=(lon_max+lon_min)/2.

    m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,\
                lon_0=lon_0, lat_0=lat_0, \
                llcrnrlon=lon_min,llcrnrlat=lat_min, \
                urcrnrlon=lon_max,urcrnrlat=lat_max,\
                rsphere=6371200.,resolution='i')

    try:
        m.drawcoastlines()
    except:
        pass
    # end try

    #m.drawstates()
    #m.drawcountries()
    m.drawparallels(np.arange(-90.,90.,1.), labels=[1,0,0,0],fontsize=8, dashes=[2, 2], color='0.5', linewidth=0.75)
    m.drawmeridians(np.arange(0.,360.,1.), labels=[0,0,0,1], fontsize=8, dashes=[2, 2], color='0.5', linewidth=0.75)

    nx = int((lon_max - lon_min)/dlonlat + 1)
    ny = int((lat_max - lat_min)/dlonlat + 1)
    lons = np.linspace(lon_min, lon_max, nx)
    lats = np.linspace(lat_min, lat_max, ny)

    glons, glats = np.meshgrid(lons, lats)
    vals = np.zeros(glons.shape)

    for i in np.arange(glons.shape[0]):
        for j in np.arange(glons.shape[1]):
            vals[i,j] = gravity.query(glons[1,j], glats[i,j])
        # end for
    # end for

    zvals,cmap = gmtColormap(cpt_file)
    cbinfo = m.pcolormesh(glons, glats, vals, latlon=True, cmap=cmap,
                          shading='auto', rasterized=True, alpha=alpha)

    if(show_colorbar): cbar = fig.colorbar(cbinfo)

    return fig, m
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('plot-type', required=True,
                type=click.Choice(['topo', 'grav'], case_sensitive=False))
@click.argument('rf-file', required=True,
                type=click.Path('r'))
@click.argument('grid', required=True,
                type=click.Path('r'))
@click.argument('cpt-file', required=True,
                type=click.Path('r'))
@click.argument('output-path', required=True,
                type=click.Path(exists=True))
@click.option('--alpha', default=1., show_default=True,
              help='Alpha value for topography/gravity plot')
@click.option('--marker-size', default=2, show_default=True,
              help='Marker size for stations')
@click.option('--label-font-size', default=5, show_default=True,
              help='Font size for station names')
@click.option('--adjust-labels', is_flag=True, default=False, show_default=True,
              help='Adjusts text labels. Note that this option can slow down performance significantly')
@click.option('--add-arrows', is_flag=True, default=False, show_default=True,
              help='Add arrows between station labels and markers. Useful for scenarios in which '
                   'station labels need to be moved far apart to reduce clutter')
@click.option('--show-colorbar', is_flag=True, default=False, show_default=True,
              help='Shows colorbar')
def process(plot_type, rf_file, grid, cpt_file, output_path, alpha,
            marker_size, label_font_size, adjust_labels, add_arrows, show_colorbar):

    """
    PLOT_TYPE : Plot-type; can be either topo (topography) or grav (gravity) \n
    RF_FILE : Receiver function data in HDF5 format \n
    GRID: Topography grid in NetCDF format or Gravity grid in ERS format \n
    CPT_FILE: GMT color palette to be used \n
    OUTPUT_PATH: Output folder \n

    example usage:
    python rf_plot_map.py topo OA_Yr2_event_waveforms.h5 au_gebco.nc mby_topo-bath.cpt tmp/

    """
    names, coords = rf_get_coords(rf_file)

    #names = names[:20]
    #coords = coords[:20]

    # initialization of map
    fig = m = None
    if(plot_type == 'topo'):
        fig, m = plot_topo(coords, grid, cpt_file, show_colorbar, alpha)
    elif(plot_type == 'grav'):
        fig, m = plot_grav(coords, grid, cpt_file, show_colorbar, alpha)
    # end if

    ax = fig.axes[0]
    lon, lat = m(coords[:,0], coords[:, 1])
    markers = ax.plot(lon, lat, 'y^', markeredgecolor='k', markeredgewidth=0.2, \
                       markersize=marker_size)
    labels = []
    for name, x, y in zip(names, lon, lat):
        name = name.split('.')[1]
        labels.append(plt.text(x, y, name, fontdict={'fontsize':label_font_size}))
    # end for
    
    if(adjust_labels):
        adjust_text(labels, 
                    lon, lat, ax=ax, lim=50,
                    add_object=markers)
    # end if

    if(add_arrows):
        for t, xi, yi in zip(labels, lon, lat):
            ax.annotate('', xy=(xi, yi), xytext=t.get_position(),
                        arrowprops=dict(
                            arrowstyle="-",
                            color='red',
                            lw=0.2,
                            mutation_scale=2,
                            shrinkA=2,
                            shrinkB=2))
    # end if

    net = names[0].split('.')[0]

    fname = os.path.join(output_path, net + '-{}-map.pdf'.format(plot_type))
    plt.savefig(fname)
    plt.close()
# end func

#-------------Main---------------------------------

if __name__=='__main__':
    """ It is an example of how to plot nice maps """

    process()
# end if