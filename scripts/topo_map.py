"""
Alexei's example code, together with data in /g/data/ha3/SHARED_DATA
$ ls
cpt  GEBCO  Raster_Mosaic  SRTM_data  topo_map.py  topo_map.py~

Requires additional modules for imports
Date: 2018-05-23
"""
from os import path, walk, system
from matplotlib.mlab import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource
from numpy import arange, mean, percentile, array, unique, where, argsort, floor, ceil
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from misc_tools import remove_last_cmap_colour

plt.rcParams['pdf.fonttype'] = 42

fig = plt.figure(figsize=(18,10))
plt.tick_params(labelsize=16)
ax = fig.add_subplot(111)

m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='h',area_thresh=2000.)

# draw coastlines, state and country boundaries, edge of map.

m.drawcoastlines()
m.drawstates()
m.drawcountries()
m.drawparallels(arange(-90.,90.,1.), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,1.), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)

m.drawmapscale(144, -34.8, 146., -38.5, 400, fontsize = 16, barstyle='fancy', zorder=100)

nc = NetCDFFile('/nas/active/ops/aap_seismic/Passive/SHARED_DATA/GEBCO/au_gebco.nc')

zscale =20. #gray
zscale =50. #colour
data = nc.variables['elevation'][:] / zscale
lons = nc.variables['lon'][:]
lats = nc.variables['lat'][:]

# transform to metres
nx = int((m.xmax-m.xmin)/500.)+1
ny = int((m.ymax-m.ymin)/500.)+1

topodat = m.transform_scalar(data,lons,lats,nx,ny)

cptfile = '/nas/active/ops/aap_seismic/Passive/SHARED_DATA/cpt/mby_topo-bath.cpt'

cmap, zvals = cpt2colormap(cptfile, 256)
cmap = remove_last_cmap_colour(cmap)
        
# make shading

ls = LightSource(azdeg = 180, altdeg = 45)
norm = mpl.colors.Normalize(vmin=-8000/zscale, vmax=5000/zscale)#myb
rgb = ls.shade(topodat, cmap=cmap, norm=norm)
im = m.imshow(rgb)

plt.show()

