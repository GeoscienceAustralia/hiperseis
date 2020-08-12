#!/usr/bin/env python
"""
Use cartopy to plot point dataset onto map.
"""

import os

import click
import numpy as np
import matplotlib.pyplot as plt
import cartopy as cp
import numpy as np
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable

COLORMAP = "./moho_kennett.cpt"

def plot_spatial_map(point_dataset, projection_code=None, gradient=False, title=None, 
                     feature_label=None, bounds=None, scale=None):
    """
    Make spatial plot of point dataset with filled contours overlaid on map.

    :param point_dataset: Name of point dataset file. Should be in format produced by
        script `pointsets2grid.py`
    :param projection_code: EPSG projection code, e.g. 3577 for Australia
    :param title: Title string for top of plot
    :param feature_label: Label for the color bar of the plotted feature
    :param bounds: 4 element tuple of (L, R, B, T) for limiting plot extent
    :param scale: 2 element tuple of (vmin, vmax) for limiting colormap scale
    :return:
    """
    with open(point_dataset, 'r') as f:
        nx = int(f.readline())
        ny = int(f.readline())
        ds = np.loadtxt(f, delimiter=',')
    # end with
    # For each column in ds, reshape((ny, nx))
    x = ds[:, 0].reshape((ny, nx))
    y = ds[:, 1].reshape((ny, nx))
    z = ds[:, 2].reshape((ny, nx))

    # Source data is defined in lon/lat coordinates.
    data_crs = cp.crs.PlateCarree()

    if projection_code is not None:
        map_projection = cp.crs.epsg(projection_code)
    else:
        map_projection = data_crs

    resolution = '50m'
    if os.path.exists(COLORMAP):
        _, cmap = _gmt_colormap(COLORMAP)
    else:
        cmap = 'magma_r'

    # Figure out bounds in map coordinates
    if bounds:
        l, b, r, t = bounds
        xy_min = map_projection.transform_point(l, b, data_crs)
        xy_max = map_projection.transform_point(r, t, data_crs)
    else:
        xy_map = map_projection.transform_points(data_crs, x.flatten(), y.flatten())[:, :2]
        xy_min = xy_map.min(axis=0)
        xy_max = xy_map.max(axis=0)
        span = xy_max - xy_min
        xy_min -= 0.2*span
        xy_max += 0.2*span

    
    fig = plt.figure(figsize=(16,9))
    ax = plt.axes(projection=map_projection)
    ax.set_xlim(xy_min[0], xy_max[0])
    ax.set_ylim(xy_min[1], xy_max[1])
    _gridliner = ax.gridlines(draw_labels=True, linestyle=':')
    ax.add_feature(cp.feature.OCEAN.with_scale(resolution))
    ax.add_feature(cp.feature.LAND.with_scale(resolution))
    ax.add_feature(cp.feature.STATES.with_scale(resolution), linewidth=0.7,
                   edgecolor="#808080a0", antialiased=True)

    if scale:
        norm = Normalize(scale[0], scale[1])
        sm = ScalarMappable(norm, cmap)
    else:
        norm = None
        sm = None

    if gradient:
        u, v = np.gradient(z)
        ax.quiver(x, y, v, u, transform=data_crs, zorder=3, angles='xy')
    else:
        _contr = ax.contourf(x, y, z, levels=20, transform=data_crs, cmap=cmap, norm=norm,
                             alpha=0.8, zorder=2)
        sm = sm if sm is not None else _contr
        _cb = plt.colorbar(sm, pad=0.1, shrink=0.6, fraction=0.05)
        _cb.ax.invert_yaxis()
        if feature_label is not None:
            _cb.set_label(feature_label)
    if title is not None:
        plt.title(title)
    return fig
# end func


@click.command()
@click.option('--gradient', is_flag=True)
@click.option('--projection-code', type=int, required=False,
              help='EPSG projection code to reproject map to, e.g. 3577 for Australia')
@click.option('--bounds', nargs=4, type=float, required=False,
              help='Bounding box for limiting map extents, of format xmin, ymin, xmax, ymax')
@click.option('--scale', nargs=2, type=float, required=False,
              help='vmin, vmax for colormap scaling')
@click.argument('point-dataset', type=click.Path(dir_okay=False, exists=True),
                required=True)
@click.argument('output-file', type=click.Path(dir_okay=False, exists=False),
                required=False)
def main(point_dataset, gradient, projection_code=None, bounds=None, scale=None,
         output_file=None):
    _f = plot_spatial_map(point_dataset, 
                          projection_code=projection_code,
                          gradient=gradient,
                          title='Moho depth from blended data',
                          feature_label='Moho depth (km)',
                          bounds=bounds,
                          scale=scale)
    if output_file is not None:
        _, ext = os.path.splitext(output_file)
        assert ext and ext.lower() in ['.png', '.pdf'], 'Provide output file extension to specify output format!'
        plt.savefig(output_file, dpi=300)
        print('Saved plot in file', output_file)
    # end if
    plt.show()
    plt.close()
# end func

def _gmt_colormap(filename):
      """ Borrowed from AndrewStraw, scipy-cookbook
          this subroutine converts GMT cpt file to matplotlib palette """

      import colorsys

      try:
          f = open(filename)
      except:
          print("file ",filename, "not found")
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
      return (x,LinearSegmentedColormap('my_colormap',colorDict,255))

if __name__ == '__main__':
    main()
# end if



