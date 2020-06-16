#!/usr/bin/env python
"""

"""

import click
import numpy as np
import matplotlib.pyplot as plt
import cartopy as cp


def plot_spatial_map(point_dataset, projection_code, title=None, feature_label=None):
    """

    :param point_dataset:
    :param projection_code:
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

    _f = plt.figure(figsize=(16,9))
    ax = plt.axes(projection=map_projection)
    ax.set_xlim(xy_min[0], xy_max[0])
    ax.set_ylim(xy_min[1], xy_max[1])
    # ax.stock_img()
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
    plt.show()
    plt.close()
# end func


@click.command()
@click.option('--projection-code', type=int, required=True)
@click.argument('point-dataset', type=click.Path(dir_okay=False, exists=True),
                required=True)
def main(point_dataset, projection_code):
    plot_spatial_map(point_dataset, projection_code,
                     title='Moho depth from blended data',
                     feature_label='Moho depth (km)')
# end func


if __name__ == '__main__':
    main()
# end if
