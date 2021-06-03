"""
Generates geotiffs of the interpolated moho grid and gradient,
and shapefiles of the station/method locations.
"""
import os

import numpy as np
import rasterio
import shapefile

CRS = rasterio.crs.CRS.from_epsg(4326)

def _profile(data, nx, ny, bands=1, bounds=None):
    """
    Profile for writing depth and gradient. Dtype and band count needs
    to be set depending on data.
    """
    if bounds is not None:
        l, b, r, t = bounds
    else:
        l, b = np.min(data[:, 0]), np.min(data[:, 1])
        r, t = np.max(data[:, 0]), np.max(data[:, 1])

    with rasterio.Env():
        profile = rasterio.profiles.DefaultGTiffProfile()
        transform = rasterio.transform.from_bounds(l, b, r, t, nx, ny)
        profile.update(crs=CRS, transform=transform, width=nx, height=ny,
                       count=bands, dtype=data.dtype)

    return profile


def from_params(params):
    if not os.path.exists(params.gis_dir):
        os.mkdir(params.gis_dir)

    write_depth_grid(params.grid_data, params.bounds, params.gis_grid_file)
    write_gradient_grid(params.grad_data, params.bounds, params.gis_grad_file)
    write_sample_locations(params.method_datasets, params.gis_loc_file)


def write_depth_grid(grid_data, bounds, outfile):
    """
    Writes the interpolated depth grid as a geotiff.
    """
    print("Writing depth grid geotiff")

    with open(grid_data, 'r') as fr:
        nx = int(fr.readline())
        ny = int(fr.readline())
        grid_ds = np.loadtxt(fr, delimiter=',')

    gtiff_profile = _profile(grid_ds, nx, ny, bands=1, bounds=bounds)

    with rasterio.Env():
        # GDAL origin is top-left, so we need to flip the data so first element is top-left cell
        data = np.flipud(grid_ds[:, 2].reshape((ny, nx)))
        gtiff_profile.update(count=1, dtype=data.dtype)
        with rasterio.open(outfile, 'w', **gtiff_profile) as dst:
            dst.write(data, 1)

    print(f"Complete! File saved to '{outfile}'")


def write_gradient_grid(grad_data, bounds, outfile):
    """
    Writes the gradient grid as a two band raster, first band is U
    components and second band is V components.
    """
    print("Writing gradient grid geotiff")

    with open(grad_data, 'r') as fr:
        nx = int(fr.readline())
        ny = int(fr.readline())
        grad_ds = np.loadtxt(fr, delimiter=',')

    gtiff_profile = _profile(grad_ds, nx, ny, bands=2, bounds=bounds)

    with rasterio.Env():
        # GDAL origin is top-left, so we need to flip the data so first element is top-left cell
        u_data = np.flipud(grad_ds[:, 2].reshape((ny, nx)))
        v_data = np.flipud(grad_ds[:, 3].reshape((ny, nx)))
        gtiff_profile.update(count=2, dtype=u_data.dtype)
        with rasterio.open(outfile, 'w', **gtiff_profile) as dst:
            dst.write(u_data, 1)
            dst.write(v_data, 2)

    print(f"Complete! File saved to '{outfile}'")


def write_sample_locations(methods, outfile):
    print("Writing location shapefile")

    for data in methods:
        w = shapefile.Writer(outfile.format(data.name), shapeType=1)
        w.field('WEIGHT', 'N', decimal=2)
        w.field('DEPTH', 'N', decimal=2)
        w.field('STA', 'C')
        if data.net is None:
            net = np.chararray(data.val.shape, 3)
            net.fill('N/A')
        else:
            net = data.net
        if data.sta is None:
            sta = np.chararray(data.val.shape, 3)
            sta.fill('N/A')
        else:
            sta = data.sta
        formatted_data = np.array((data.lon, data.lat, data.total_weight, data.val, net, sta)).T
        for d in formatted_data:
            w.point(float(d[0]), float(d[1]))
            w.record(WEIGHT=d[2], DEPTH=d[3], STA='.'.join((str(d[4]), str(d[5]))))
        w.close()
        # Write .prj file
        with open(f'{outfile.format(data.name)}.prj', 'w') as prj:
            epsg = 'GEOGCS["WGS 84",'
            epsg += 'DATUM["WGS_1984",'
            epsg += 'SPHEROID["WGS 84",6378137,298.257223563]]'
            epsg += ',PRIMEM["Greenwich",0],'
            epsg += 'UNIT["degree",0.0174532925199433]]'
            prj.write(epsg)

    print(f"Complete! Location shapefiles written to '{outfile.format('method_name')}'")
