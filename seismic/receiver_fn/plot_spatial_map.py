#!/usr/bin/env python
"""
Use cartopy to plot moho grid and gradient onto a map.
"""
import json
import os

import click
import numpy as np
import matplotlib.pyplot as plt
import cartopy as cp
import numpy as np
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable


COLORMAP = "./moho_kennett.cpt"


def plot_spatial_map(grid_data, gradient_data, projection_code=None, 
                     title=None, feature_label=None, bounds=None, scale=None):
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
    with open(grid_data, 'r') as f:
        nx = int(f.readline())
        ny = int(f.readline())
        grid_ds = np.loadtxt(f, delimiter=',')
    with open(gradient_data, 'r') as f:
        grad_ds = np.loadtxt(f, delimiter=',', skiprows=2)

    # For each column in data, reshape((ny, nx))
    x = grid_ds[:, 0].reshape((ny, nx))
    y = grid_ds[:, 1].reshape((ny, nx))
    z = grid_ds[:, 2].reshape((ny, nx))
    u = grad_ds[:, 2].reshape((ny, nx))
    v = grad_ds[:, 3].reshape((ny, nx))

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
        xy_min -= 0.1*span
        xy_max += 0.1*span

    fig = plt.figure(figsize=(14,12))

    cont_ax = plt.subplot(2, 1, 1, projection=map_projection)
    grad_ax = plt.subplot(2, 1, 2, projection=map_projection)

    for ax in [cont_ax, grad_ax]:
        ax.set_xlim(xy_min[0], xy_max[0])
        ax.set_ylim(xy_min[1], xy_max[1])
        ax.gridlines(draw_labels=True, linestyle=':')
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

    contr = cont_ax.contourf(x, y, z, levels=20, transform=data_crs, cmap=cmap, norm=norm,
                              alpha=0.8, zorder=2)
    sm = sm if sm is not None else contr
    cont_div = make_axes_locatable(cont_ax)
    cax = cont_div.append_axes('bottom', size='5%', pad='7%', axes_class=plt.Axes)
    cb = plt.colorbar(sm, cax=cax, orientation="horizontal")
    if feature_label is not None:
        cb.set_label(feature_label)
    
    grad_ax.quiver(x, y, u, v, transform=data_crs, zorder=2, angles='xy', units='xy')
    # Hack: shrink gradient ax so it's the same size as contour ax (contour ax shrinks due to 
    # appending colorbar ax)
    grad_div = make_axes_locatable(grad_ax)
    grad_div.append_axes('bottom', size='5%', pad='7%', axes_class=plt.Axes).set_visible(False)

    if title is not None:
        cont_ax.set_title(title)

    return fig


def from_config(config_file):
    """
    Create plots from config as part of moho workflow.
    """
    print("Plotting Moho grid and gradient map")
    with open(config_file, mode='r') as f:
        job_config = json.load(f)

    plotting = job_config['plotting']
    cp = plotting.get('cartopy_parameters', {})
    scale = cp.get('scale')
    fmt = cp.get('format', 'png')
    show = cp.get('show', False)
    title = cp.get('title', 'Moho depth from blended data')
    cb_label = cp.get('cb_label', 'Moho depth (km)')
    outdir = job_config.get('output_dir', os.getcwd())
    grid_data = os.path.join(outdir, 'moho_grid.csv')
    gradient_data = os.path.join(outdir, 'moho_gradient.csv')
    fig = plot_spatial_map(grid_data, gradient_data, scale=scale, 
                           title=title, feature_label=cb_label)
    if show:
        print("Showing plot, close display window to continue")
        plt.show()
    outfile = os.path.join(outdir, f'moho_plot.{fmt}')
    fig.savefig(outfile, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Complete! Plot saved to '{outfile}'")
                

@click.command()
@click.option('--projection-code', type=int, required=False,
              help='EPSG projection code to reproject map to, e.g. 3577 for Australia')
@click.option('--bounds', nargs=4, type=float, required=False,
              help='Bounding box for limiting map extents, of format xmin, ymin, xmax, ymax')
@click.option('--scale', nargs=2, type=float, required=False,
              help='vmin, vmax for colormap scaling')
@click.argument('grid-data', type=click.Path(dir_okay=False, exists=True),
                required=True)
@click.argument('grad-data', type=click.Path(dir_okay=False, exists=True),
                required=True)
@click.argument('output-file', type=click.Path(dir_okay=False, exists=False),
                required=False)
def main(grid_data, grad_data, projection_code=None, bounds=None, scale=None,
         output_file=None):
    plot_spatial_map(grid_data, grad_data,
                     projection_code=projection_code,
                     title='Moho depth from blended data',
                     feature_label='Moho depth (km)',
                     bounds=bounds,
                     scale=scale)

    if output_file is not None:
        _, ext = os.path.splitext(output_file)
        assert ext and ext.lower() in ['.png', '.pdf'], \
            'Provide output file extension to specify output format!'
        plt.savefig(output_file, dpi=300)
        print('Saved plot in file', output_file)

    plt.show()
    plt.close()


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
