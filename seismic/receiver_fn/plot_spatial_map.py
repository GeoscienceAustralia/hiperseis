#!/usr/bin/env python
"""
Use cartopy to plot point dataset onto map.
"""

import os

import click
import numpy as np
import matplotlib.pyplot as plt
import cartopy as cp


def plot_spatial_map(point_dataset, projection_code, title=None, feature_label=None):
    """
    Make spatial plot of point dataset with filled contours overlaid on map.

    :param point_dataset: Name of point dataset file. Should be in format produced by
        script `pointsets2grid.py`
    :param projection_code: EPSG projection code, e.g. 3577 for Australia
    :param title: Title string for top of plot
    :param feature_label: Label for the color bar of the plotted feature
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
    cmap = 'magma_r'

    # Figure out bounds in map coordinates
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
    # # Compute map extents in lon/lat for defining range of background image
    # p0 = data_crs.transform_point(xy_min[0], xy_min[1], map_projection)
    # p1 = data_crs.transform_point(xy_min[0], xy_max[1], map_projection)
    # p2 = data_crs.transform_point(xy_max[0], xy_min[1], map_projection)
    # p3 = data_crs.transform_point(xy_max[0], xy_max[1], map_projection)
    # map_bounds_lonlat = np.array([p0, p1, p2, p3])
    # map_bounds_lonlat_min = map_bounds_lonlat.min(axis=0)
    # map_bounds_lonlat_max = map_bounds_lonlat.max(axis=0)
    # ax.background_img(resolution='med',
    #                   extent=[map_bounds_lonlat_min[0], map_bounds_lonlat_max[0],
    #                           map_bounds_lonlat_min[1], map_bounds_lonlat_max[1]])
    _gridliner = ax.gridlines(draw_labels=True, linestyle=':')
    ax.add_feature(cp.feature.OCEAN.with_scale(resolution))
    ax.add_feature(cp.feature.LAND.with_scale(resolution))
    ax.add_feature(cp.feature.STATES.with_scale(resolution), linewidth=0.7,
                   edgecolor="#808080a0", antialiased=True)
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
@click.argument('point-dataset', type=click.Path(dir_okay=False, exists=True),
                required=True)
@click.argument('output-file', type=click.Path(dir_okay=False, exists=False),
                required=False)
def main(point_dataset, projection_code, output_file=None):
    _f = plot_spatial_map(point_dataset, projection_code,
                          title='Moho depth from blended data',
                          feature_label='Moho depth (km)')
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
