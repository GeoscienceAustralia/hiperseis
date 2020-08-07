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

from seismic.receiver_fn.legacy.plot_map import gmtColormap

def plot_spatial_map(point_dataset, projection_code, gradient=False, title=None, 
                     feature_label=None, cpt_colormap=None, bounds=None):
    """
    Make spatial plot of point dataset with filled contours overlaid on map.

    :param point_dataset: Name of point dataset file. Should be in format produced by
        script `pointsets2grid.py`
    :param projection_code: EPSG projection code, e.g. 3577 for Australia
    :param title: Title string for top of plot
    :param feature_label: Label for the color bar of the plotted feature
    :param cpt_colormap: A CPT file to use as the matplotlib colormap
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

    map_projection = cp.crs.epsg(projection_code)
    resolution = '50m'
    if cpt_colormap is not None:
        _, cmap = gmtColormap(cpt_colormap)
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
    if gradient:
        u, v = np.gradient(z)
        ax.quiver(x, y, v, u, transform=data_crs, zorder=3, angles='xy')
    else:
        _contr = ax.contourf(x, y, z, levels=20, transform=data_crs, cmap=cmap,
                             alpha=0.8, zorder=2)
        _cb = plt.colorbar(_contr, pad=0.1, shrink=0.6, fraction=0.05)
        _cb.ax.invert_yaxis()
        if feature_label is not None:
            _cb.set_label(feature_label)
    if title is not None:
        plt.title(title)
    return fig
# end func


@click.command()
@click.option('--projection-code', type=int, required=True,
              help='EPSG projection code, e.g. 3577 for Australia')
@click.option('--gradient', is_flag=True)
@click.option('--cpt-colormap', type=click.Path(dir_okay=False, exists=True))
@click.option('--bounds', nargs=4, type=float, required=False,
              help='Bounding box for limiting map extents, of format xmin, ymin, xmax, ymax')
@click.argument('point-dataset', type=click.Path(dir_okay=False, exists=True),
                required=True)
@click.argument('output-file', type=click.Path(dir_okay=False, exists=False),
                required=False)
def main(point_dataset, projection_code, gradient, cpt_colormap=None, bounds=None, 
         output_file=None):
    _f = plot_spatial_map(point_dataset, projection_code, 
                          gradient=gradient,
                          title='Moho depth from blended data',
                          feature_label='Moho depth (km)',
                          cpt_colormap=cpt_colormap,
                          bounds=bounds)
    if output_file is not None:
        _, ext = os.path.splitext(output_file)
        assert ext and ext.lower() in ['.png', '.pdf'], 'Provide output file extension to specify output format!'
        plt.savefig(output_file, dpi=300)
        print('Saved plot in file', output_file)
    # end if
    plt.show()
    plt.close()
# end func


if __name__ == '__main__':
    main()
# end if
