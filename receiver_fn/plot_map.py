from netCDF4 import Dataset as NetCDFFile
from matplotlib.colors import LightSource
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
#from mpl_toolkits.basemap import Basemap
from mympl_toolkits.basemap import Basemap

def gmtColormap(fileName):
      """ Borrowed from AndrewStraw, scipy-cookbook
          this subroutine converts GMT cpt file to matplotlib palette """

      import colorsys

      try:
          f = open(fileName)
      except:
          print "file ",fileName, "not found"
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
      return (x,mpl.colors.LinearSegmentedColormap('my_colormap',colorDict,255))




def plot_map(data):

    fig=plt.figure(figsize=(11.69,8.27))
    plt.tick_params(labelsize=12)
    
    ax = fig.add_subplot(111)
    data=np.array(data)
    lon_min=min(data[:,0])-1.
    lon_max=max(data[:,0])+1.

    lat_1=min(data[:,1])-1.
    lat_2=min(data[:,1])
    lat_min=min(data[:,1])-1.
    lat_max=max(data[:,1])+1.
    
    lat_0=(lat_max+lat_min)/2.
    lon_0=(lon_max+lon_min)/2.

    m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,\
                lon_0=lon_0, lat_0=lat_0, \
                llcrnrlon=lon_min,llcrnrlat=lat_min, \
                urcrnrlon=lon_max,urcrnrlat=lat_max,\
                rsphere=6371200.,resolution='h',area_thresh=2000.)
    
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    m.drawparallels(np.arange(-90.,90.,1.), labels=[1,0,0,0],fontsize=14, dashes=[2, 2], color='0.5', linewidth=0.75)
    m.drawmeridians(np.arange(0.,360.,1.), labels=[0,0,0,1], fontsize=14, dashes=[2, 2], color='0.5', linewidth=0.75)

    
    
#   m.drawmapscale(lon_min, lat_max, lon_max, lat_min, 400, fontsize = 16, barstyle='fancy', zorder=100)
    nc = NetCDFFile('/g/data/ha3/Passive/SHARED_DATA/GEBCO/au_gebco.nc')
    
    zscale =20. #gray
    zscale =50. #colour
    data = nc.variables['elevation'][:] / zscale
    lons = nc.variables['lon'][:]
    lats = nc.variables['lat'][:]
    
    # transform to metres
    nx = int((m.xmax-m.xmin)/500.)+1
    ny = int((m.ymax-m.ymin)/500.)+1
    
    topodat = m.transform_scalar(data,lons,lats,nx,ny)
    
    cptfile = '/g/data/ha3/Passive/SHARED_DATA/cpt/mby_topo-bath.cpt'
    
    zvals,cmap = gmtColormap(cptfile)
    
    # make shading
    
    ls = LightSource(azdeg = 180, altdeg = 45)
    norm = mpl.colors.Normalize(vmin=-8000/zscale, vmax=5000/zscale)#myb
    rgb = ls.shade(topodat, cmap=cmap, norm=norm)
    im = m.imshow(rgb)
    return m


""" It is an example of how to plot nice maps """

data=[[144,-34.8],[146,-38.5]]
data=np.array(data)
m=plot_map(data)
lon,lat=m(data[:,0],data[:,1])
plt.plot(lon,lat,'yo',markeredgecolor='k', markeredgewidth=0.5, \
          markersize=10.0,  alpha=0.7)
plt.show()
