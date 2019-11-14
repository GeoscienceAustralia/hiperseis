#!/usr/bin/env python
"""
Generate common conversion point (CCP) plot as per *C.Sippl, "Moho geometry along a north-south passive seismic
transect through Central Australia", Technophysics 676 (2016), pp.56-69, 
DOI https://doi.org/10.1016/j.tecto.2016.03.031*

This code adapted from Christian Sippl's original code.

Workflow:
    extract_event_traces.py --> generate_rf.py --> rf_quality_filter.py --> plot_ccp.py (this script)

Example usage:
    python seismic/receiver_fn/plot_ccp.py --start-latlon -19.5 133.0 --end-latlon -19.5 140.0 --width 120 \
        --channels T --stacked-scale 0.3 --title "Network OA CCP T-stacking (profile BS24-CF24)" \
        /software/hiperseis/seismic/receiver_fn/DATA/OA-ZRT-cleaned.h5 \
        /software/hiperseis/seismic/receiver_fn/DATA/OA-ZRT-T_CCP_stack_BS24-CF24_2km_spacing.png
"""

# pylint: disable=too-many-locals,too-many-arguments,invalid-name

import os
import numpy as np
# from future.utils import iteritems

import click
from obspy.taup import TauPyModel
from obspy.taup.taup_create import TauPCreate
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from mpl_toolkits import basemap
import rf
from tqdm import tqdm

from seismic.receiver_fn.rf_util import KM_PER_DEG
from seismic.ASDFdatabase import FederatedASDFDataSet


def plot_ccp(matrx, length, max_depth, spacing, ofile=None, vlims=None, metadata=None, title=None):
    """Plot results of CCP stacking procedure.

    :param matrx: [description]
    :type matrx: numpy.array
    :param length: Length (km) of the transect line
    :type length: float
    :param max_depth: Maximum depth (km) of the slice to plot
    :type max_depth: float
    :param spacing: Grid spacing (km)
    :type spacing: float
    :param ofile: Output file name to save to plot to, defaults to None
    :type ofile: str, optional
    :param vlims: Min and max values for the color scale, defaults to None
    :type vlims: tuple(float, float), optional
    :param metadata: Metadata generated by ccp_compute_station_params(), defaults to None
    :type metadata: dict, optional
    :param title: Title text to add to the plot, defaults to None
    :type title: str, optional
    """
    tickstep_x = 50
    tickstep_y = 25

    plt.figure(figsize=(16, 9))
    interpolation = 'hanning'
    extent = (0, length, 0, max_depth)
    assert not np.any(np.isnan(matrx))
    if vlims is not None:
        im = plt.imshow(matrx, cmap='seismic', aspect='auto', vmin=vlims[0], vmax=vlims[1], extent=extent,
                        interpolation=interpolation, origin='lower')
    else:
        im = plt.imshow(matrx, cmap='seismic', aspect='auto', extent=extent, interpolation=interpolation,
                        origin='lower')
    # end if
    cb = plt.colorbar(im)
    cb.set_label('Stacked amplitude (arb. units)')

    if title is not None:
        plt.title(title, fontsize=16, y=1.02)

    plt.xlim(0, length)
    plt.ylim(max_depth*1.0001, 0)

    plt.xlabel('Distance (km)', fontsize=12)
    plt.ylabel('Depth (km)', fontsize=12)

    plt.xticks(np.arange(0.0, length*1.0001, tickstep_x), fontsize=12)
    plt.yticks(np.arange(0, max_depth*1.0001, tickstep_y), fontsize=12)
    plt.tick_params(right=True, labelright=True, axis='y', labelsize=12)
    plt.gca().yaxis.set_minor_locator(MultipleLocator(1))
    plt.gca().yaxis.set_tick_params(which='minor', right=True)
    plt.grid(color='#80808080', linestyle=':', axis='y')

    if metadata is not None:
        stn_labels = []
        # Extract station label data into a list first so we can sort by distance along transect
        for stn, meta in metadata.items():
            if meta is None:
                continue
            stn_labels.append((stn, meta))
        # end for
        stn_labels.sort(key=lambda md: md[1]['sta_offset'])
        i = 0
        for stn, meta in stn_labels:
            x = meta['sta_offset']
            y = 1.0 + (i%2)*2.5
            plt.text(x, y, "{} ({})".format(stn, meta['event_count']), horizontalalignment='center',
                     verticalalignment='top', fontsize=9, backgroundcolor='#ffffffa0')
            i += 1
        # end for
    # end if

    if ofile:
        plt.savefig(ofile, dpi=300)
    else:
        plt.show()
    plt.close()
# end func


def setup_ccp_profile(length, spacing, maxdep):
    """Construct the grid for a CCP stacking profile

    :param length: Length (km) of the profile
    :type length: float
    :param spacing: Grid spacing (km)
    :type spacing: float
    :param maxdep: Maximum depth (km) of the slice to plot
    :type maxdep: float
    :return: Zeroed matrix and mesh coordinates
    :rtype: numpy.array, numpy.array, numpy.array
    """
    #calculate number of cells in x and y direction 
    n_y = int(round(maxdep / spacing, 0))
    n_x = int(round(length / spacing, 0))

    #get center values
    depstep = np.arange(0 + round(spacing / 2.0, 1), maxdep, spacing)
    lenstep = np.arange(0 + round(spacing / 2.0, 1), length, spacing)

    #create matrix
    mtrx = np.zeros([n_x, n_y])

    return mtrx, depstep, lenstep
# end func


def get_amplitude(trace, time):
    """
    retrieve amplitude value
    """
    rf_offset = trace.stats.onset - trace.stats.starttime
    indx = (time + rf_offset) * trace.stats.sampling_rate
    amp = trace.data[int(indx)]/trace.stats.amp_max
    return amp
# end func


def add_ccp_trace(trace, inc_p, matrx, matrx_entry, vmod, depstep, lenstep, sta_offset, az):
    """
    project amplitudes from all RFs onto the profile...2D rot:
    """
    # start at zero: inc_p given, inc_s needs to be calculated  
    h = 0
    c = 0
    d = 0
    tpz = 0
    tsz = 0
    rpz = 0
    rsz = 0

    for j in range(1, len(depstep)):
        c = d
        if j == 1:
            h = h_tot = 1
        else:
            h = depstep[j] - depstep[j - 1]
            h_tot += h
        # check in velocity model
        for f in range(len(vmod[0])):
            if vmod[0][f] < depstep[j]:
                d = f

        # derive P incidence from previous P incidence, then current S from current P
        inc_p = np.arcsin((np.sin(inc_p * np.pi/180.) * vmod[1][d]) / vmod[1][c]) * 180 / np.pi
        inc_s = np.arcsin((np.sin(inc_p * np.pi/180.) * vmod[2][d]) / vmod[1][d]) * 180 / np.pi

        # horizontal distances (attention: still needs azimuth normalization)
        rpz += h*np.tan(inc_p/180.*np.pi)
        rsz += h*np.tan(inc_s/180.*np.pi)

        rd = rpz - rsz

        tpz += h/np.cos(inc_p/180.*np.pi) * (1/vmod[1][d])
        tsz += h/np.cos(inc_s/180.*np.pi) * (1/vmod[2][d])

        td = np.sin(inc_p/180.*np.pi)/vmod[1][d] * rd
        # check if velocity jump, if yes get new angles
        tps = tsz + td - tpz

        amp = get_amplitude(trace, tps)

        # project, put into correct bin in matrix
        xsz = rsz * np.cos(az * np.pi / 180.)  # relative to station
        indx_x, indx_y = matrx_lookup(xsz, sta_offset, h_tot, depstep, lenstep)

        matrx[indx_x, indx_y] += amp
        matrx_entry[indx_x, indx_y] += 1

    return matrx,matrx_entry
# end func


def matrx_lookup(xsz, sta_offset, h, depstep, lenstep):
    """
    return index values for amplitude contrbution in profile matrix
    """
    distance_offset = sta_offset - xsz # because zero is in the north

    diff_x = 9999.0
    diff_y = 9999.0
    indx_x = 0
    indx_y = 0

    #find horizontal position
    for j in range(len(lenstep)):
        if abs(lenstep[j] - distance_offset) < diff_x:
            diff_x = abs(lenstep[j] - distance_offset)
            indx_x = j

    for k in range(len(depstep)):
        if abs(depstep[k] - h) < diff_y:
            diff_y = abs(depstep[k] - h)
            indx_y = k

    return indx_x, indx_y
# end func


def bounding_box(startpoint, endpoint):
    """Compute a bounding box from start and end points.

    :param startpoint: Coordinates of starting point
    :type startpoint: pair of float
    :param endpoint: Coordinates of end point
    :type endpoint: pair of float
    :return: Bounding box (left, bottom, right, top)
    :rtype: tuple(float, float, float, float)
    """
    ybig = max(startpoint[0], endpoint[0])
    ysmall = min(startpoint[0], endpoint[0])
    xbig = max(startpoint[1], endpoint[1])
    xsmall = min(startpoint[1], endpoint[1])
    return (xsmall, ysmall, xbig, ybig)
# end func


def equirectangular_projection(x0, y0, x1, y1):
    """Perform equirectangular projection of a pair of latitude, longitude coordinates to cartesian coordinates.

    This length calculation uses the forward equirectangular projection
    (https://en.wikipedia.org/wiki/Equirectangular_projection).
    (See also https://www.movable-type.co.uk/scripts/latlong.html)

    :param x0: Point 0 longitude (deg)
    :type x0: float
    :param y0: Point 0 latitude (deg)
    :type y0: float
    :param x1: Point 1 longitude (deg)
    :type x1: float
    :param y1: Point 1 latitude (deg)
    :type y1: float
    :return: Lengths of sides of rectangle and the diagonal. The diagonal is the distance between points 0 and 1.
    :rtype: float, float, float
    """
    profile_x_length = (x1 - x0) * KM_PER_DEG * np.cos((y1 + y0) / 2. * np.pi / 180.)
    profile_y_length = (y1 - y0) * KM_PER_DEG
    # This is an approximate great circle arc length
    length = np.sqrt(profile_x_length**2 + profile_y_length**2)

    return profile_x_length, profile_y_length, length
# end func


def bearing(p1, p2):
    """Compute bearing (forward azimuth) in degrees from p1 to p2.

    Math reference: https://www.movable-type.co.uk/scripts/latlong.html

    :param p1: (latitude, longitude) in degrees
    :type p1: tuple(float, float)
    :param p2: (latitude, longitude) in degrees
    :type p2: tuple(float, float)
    """
    p1_rad = (p1[0]*np.pi/180.0, p1[1]*np.pi/180.0)
    p2_rad = (p2[0]*np.pi/180.0, p2[1]*np.pi/180.0)
    delta_long = p2_rad[1] - p1_rad[1]
    y = np.sin(delta_long)*np.cos(p2_rad[0])
    x = np.cos(p1_rad[0]) * np.sin(p2_rad[0]) - np.sin(p1_rad[0]) * np.cos(p2_rad[0]) * np.cos(delta_long)
    b = np.arctan2(y, x)
    b = b * 180.0 / np.pi
    b = (b + 360.0) % 360.0
    return b
# end func


def angular_distance(p1, p2):
    """Compute the angular distance (in degrees) between two points p1 and p2 using Haversine formula.

    Math reference: https://www.movable-type.co.uk/scripts/latlong.html

    :param p1: (latitude, longitude) in degrees
    :type p1: tuple(float, float)
    :param p2: (latitude, longitude) in degrees
    :type p2: tuple(float, float)
    """
    p1_rad = (p1[0]*np.pi/180.0, p1[1]*np.pi/180.0)
    p2_rad = (p2[0]*np.pi/180.0, p2[1]*np.pi/180.0)
    delta_lat = p2_rad[0] - p1_rad[0]
    delta_long = p2_rad[1] - p1_rad[1]
    sin_dlat = np.sin(0.5*delta_lat)
    sin_dlong = np.sin(0.5*delta_long)
    a = sin_dlat*sin_dlat + np.cos(p1_rad[0]) * np.cos(p2_rad[0]) * sin_dlong * sin_dlong
    a = np.clip(a, 0.0, 1.0)
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))
    c = c * 180.0 / np.pi
    return c
# end func


def cross_along_track_distance(p1, p2, p3):
    """Compute the cross-track distance and the along-track distance of p3 in relation to great circle from p1 to p2.

    Math reference: https://www.movable-type.co.uk/scripts/latlong.html

    :param p1: (latitude, longitude) in degrees
    :type p1: tuple(float, float)
    :param p2: (latitude, longitude) in degrees
    :type p2: tuple(float, float)
    :param p3: (latitude, longitude) in degrees
    :type p3: tuple(float, float)
    :return: Cross-track distance (orthogonal to arc from p1 to p2) and along-track distance (along arc from
        p1 to p2) of the location p3.
    :rtype: tuple(float, float)
    """
    delta_13 = angular_distance(p1, p3) * np.pi / 180.0
    theta_13 = bearing(p1, p3) * np.pi / 180.0
    theta_12 = bearing(p1, p2) * np.pi / 180.0
    relative_bearing = theta_13 - theta_12
    while relative_bearing > np.pi:
        relative_bearing -= 2*np.pi
    while relative_bearing < -np.pi:
        relative_bearing += 2*np.pi
    ct_angle = np.arcsin(np.sin(delta_13) * np.sin(relative_bearing))
    at_angle = np.arccos(np.cos(delta_13) / np.cos(ct_angle))
    # If bearing from p1 to p3 is more than 90 degrees from bearing to p2, then the along-track angle
    # is negative.
    if np.abs(relative_bearing) > np.pi/2.0:
        at_angle = -at_angle
    ct_angle = ct_angle * 180.0 / np.pi
    at_angle = at_angle * 180.0 / np.pi
    return (ct_angle * KM_PER_DEG, at_angle * KM_PER_DEG)
# end func


def ccp_compute_station_params(rf_stream, startpoint, endpoint, width, bm=None):
    """Determines which stations are between startpoint and endpoint great circle profile line,
       and within *width* distance of that profile line. Generates a dictionary of distance
       along profile line (*sta_offset*) and orthogonal distance from profile line (*dist*).
       For stations outside the region of interest, dictionary has the value None.

    :param rf_stream: RFStreams whose stations will be processed.
    :type rf_stream: rf.rfstream.RFStream
    :param startpoint: (latitude, longitude) in degrees of the profile line start point.
    :type startpoint: tuple(float, float)
    :param endpoint: (latitude, longitude) in degrees of the profile line end point.
    :type endpoint: tuple(float, float)
    :param width: Width of ROI in km to either side of the profile line.
    :type width: float
    :param bm: Optional basemap upon which to plot stations coded according to whether they are within the ROI.
    :type bm: mpl_toolkits.basemap.Basemap
    :return: Dictionary keyed by station code containing distances relative to profile line, or
        None if station is outside the ROI.
    :rtype: dict
    """
    stn_params = {}
    length = angular_distance(startpoint, endpoint) * KM_PER_DEG
    print("Transect length = {} km".format(length))
    pbar = tqdm(total=len(rf_stream), ascii=True)
    for tr in rf_stream:
        pbar.update()
        stat_code = tr.stats.station
        pbar.set_description("Station {}".format(stat_code))
        if stat_code not in stn_params:
            sta_lat = tr.stats.station_latitude
            sta_lon = tr.stats.station_longitude
            (dist, sta_offset) = cross_along_track_distance(startpoint, endpoint, (sta_lat, sta_lon))
            within_dist_tolerance = (abs(dist) <= width)
            within_profile_length = (sta_offset >= 0 and sta_offset <= length)

            if within_dist_tolerance and within_profile_length:
                pbar.write("Station " + stat_code + " included: (offset, dist) = ({}, {})".format(sta_offset, dist))
                if bm is not None:
                    x, y = bm(sta_lon, sta_lat)
                    bm.plot(x, y, 'ro')
                    plt.text(x, y, stat_code, fontsize=6, color='#202020a0')
                stn_params[stat_code] = {'dist': dist, 'sta_offset': sta_offset}
            else:
                if bm is not None:
                    x, y = bm(sta_lon, sta_lat)
                    bm.plot(x, y, '^', color='#808080')
                stn_params[stat_code] = None
            # end if
        # end if
    # end for
    pbar.close()

    return stn_params
# end func


def ccp_generate(rf_stream, startpoint, endpoint, width, spacing, max_depth, channels=None, v_background='ak135',
                 station_map_file=None):
    """Main function for processing RF collection and plotting common conversion point (CCP) stack of receiver
       functions (RFs) along a specific line between startpoint and endpoint.

    :param rf_stream: Sequence of receiver function traces
    :type rf_stream: rf.rfstream.RFStream
    :param startpoint: Starting point for the transect line in latitude, longitude degrees
    :type startpoint: tuple(float, float)
    :param endpoint: End point for the transect line in latitude, longitude degrees
    :type endpoint: tuple(float, float)
    :param width: Width (km) of region about the transect line from which to project RFs to the slice.
    :type width: float
    :param spacing: Size (km) of each discrete sample cell in the spatial slice beneath the transect.
    :type spacing: float
    :param max_depth: Maximum depth (km) of the slice beneath the transect
    :type max_depth: float
    :param channels: Filter for channels, defaults to None in which case only 'R' channel is used.
        Allowed values are in ['Z', 'R', 'T'].
    :type channels: list(str), optional
    :param v_background: Assumed background 1D velocity model, defaults to 'ak135'
    :type v_background: str, optional
    :param station_map_file: File name to save station map to, defaults to None
    :type station_map_file: str, optional
    :return: Normalized stack matrix, normalization factors, profile length and station metadata
    :rtype: numpy.array, numpy.array, float, dict
    """
    if channels is None:
        channels = ['R']

    xsmall, ysmall, xbig, ybig = bounding_box(startpoint, endpoint)
    _, _, length = equirectangular_projection(xsmall, ysmall, xbig, ybig)

    #set velocity model (other options can be added)
    # AK135
    vp = [5.8, 6.5, 8.04, 8.045, 8.05, 8.175, 8.3, 8.4825, 8.665,
          8.8475, 9.03, 9.36, 9.528, 9.696, 9.864, 10.032, 10.2]
    vs = [3.46, 3.85, 4.48, 4.49, 4.5, 4.509, 4.518, 4.609, 4.696,
          4.783, 4.87, 5.08, 5.186, 5.292, 5.398, 5.504, 5.61]
    if v_background.lower() == 'ak135':
        # NOTE: Conventional default AK135
        b_lay = [0., 20., 35., 77.5, 120., 165., 210., 260.,
                 310., 360., 410., 460., 510., 560., 610., 660.]
        model = TauPyModel(model=v_background.lower())
    elif v_background.lower() == 'ak135_60':
        # NOTE: AK135 with Moho at 60 km
        b_lay = [0., 20., 60., 77.5, 120., 165., 210., 260.,
                 310., 360., 410., 460., 510., 560., 610., 660.]
        model_dir = os.path.join(os.path.split(__file__)[0], 'models')
        model_file = os.path.join(model_dir, 'ak135_60.tvel')
        built_file = os.path.join(model_dir, 'ak135_60.npz')
        print('Using background crustal model at {}'.format(model_file))
        if not os.path.isfile(built_file):
            model_factory = TauPCreate(model_file, built_file)
            model_factory.load_velocity_model()
            model_factory.run()
        # end if
        model = TauPyModel(model=built_file)
    else:
        assert False, 'NYI'
    #end if
    vmod = (b_lay, vp, vs)

    #first set up a matrix
    profile_mesh, depstep, lenstep = setup_ccp_profile(length, spacing, max_depth)

    # snapshot of initial matrix
    mesh_entries = profile_mesh.copy()

    #create map plot file
    expand_width = 1.0 + width/KM_PER_DEG
    m = basemap.Basemap(projection='merc', urcrnrlat=ybig + expand_width, urcrnrlon=xbig + expand_width,
                        llcrnrlon=xsmall - expand_width, llcrnrlat=ysmall - expand_width, resolution='i')
    m.drawcoastlines()
    x1, y1 = m(startpoint[1], startpoint[0])
    x2, y2 = m(endpoint[1], endpoint[0])
    m.plot([x1, x2], [y1, y2], 'r--')
    parallels = np.arange(np.floor(ysmall - 2), np.ceil(ybig + 2))
    m.drawparallels(parallels, color="#a0a0a0", labels=[1, 1, 0, 0], fontsize=8)
    meridians = np.arange(np.floor(xsmall - 2), np.ceil(xbig + 2))
    m.drawmeridians(meridians, rotation=45, color="#a0a0a0", labels=[0, 0, 1, 1], fontsize=8)

    # Precompute the station parameters for a given code, as this is the same for every trace.
    print("Computing included stations...")
    stn_params = ccp_compute_station_params(rf_stream, startpoint, endpoint, width, m)

    # Processing/extraction of rf_stream data
    print("Projecting included stations to slice...")
    pbar = tqdm(total=len(rf_stream), ascii=True)
    for tr in rf_stream:
        pbar.update()
        stat_code = tr.stats.station
        pbar.set_description("{} event {}".format(stat_code, tr.stats.event_id))

        if tr.stats.channel[-1] in channels and stn_params[stat_code] is not None:
            # Updating of matrix
            sta_offset = stn_params[stat_code]['sta_offset']
            try:
                back_azi = tr.stats.back_azimuth
                arr_p = model.get_travel_times(tr.stats.event_depth, tr.stats.distance, phase_list=['P'])
                inc_p = arr_p[0].incident_angle
                profile_mesh, mesh_entries = add_ccp_trace(tr, inc_p, profile_mesh, mesh_entries, vmod,
                                                           depstep, lenstep, sta_offset, back_azi)
                if 'event_count' in stn_params[stat_code]:
                    stn_params[stat_code]['event_count'] += 1
                else:
                    stn_params[stat_code]['event_count'] = 1
            except IndexError as err:
                print(err)
                continue
        # end if
    # end for
    pbar.close()

    # Plot parameters on the map
    ax = plt.gca()
    plt.text(0.01, 0.98, 'Start = ({:3.3f},{:3.3f}) deg'.format(*startpoint), verticalalignment='top',
             transform=ax.transAxes, fontsize=6, backgroundcolor='#ffffffa0')
    plt.text(0.01, 0.93, 'End = ({:3.3f},{:3.3f}) deg'.format(*endpoint), verticalalignment='top',
             transform=ax.transAxes, fontsize=6, backgroundcolor='#ffffffa0')
    plt.text(0.01, 0.88, 'Width = {:3.1f} km'.format(width), verticalalignment='top',
             transform=ax.transAxes, fontsize=6, backgroundcolor='#ffffffa0')

    # Show or save basemap plot
    if station_map_file is None:
        plt.show()
    else:
        plt.savefig(station_map_file, dpi=300)
    plt.close()

    if not np.all(mesh_entries[:] == 0):
        #normalize, then plot
        matrx_norm = (profile_mesh / mesh_entries).transpose()
        matrx_norm[np.isnan(matrx_norm)] = 0
        return matrx_norm, mesh_entries.transpose(), length, stn_params
    else:
        return None, None, 0, stn_params
    # end if
# end func


def run(rf_stream, output_file, start_latlon, end_latlon, width, spacing, max_depth, channels,
        background_model='ak135', stacked_scale=None, title=None, plot_density=False):
    """Run CCP generation on a given dataset of RFs.

    :param rf_stream: Set of RFs to use for CCP plot
    :type rf_stream: rf.RFStream
    :param output_file: Name of output file (png format)
    :type output_file: str or Path
    :param start_latlon: Starting (latitude, longitude) coordinates of line transect
    :type start_latlon: tuple(float, float)
    :param end_latlon: End (latitude, longitude) coordinates of line transect
    :type end_latlon: tuple(float, float)
    :param width: Width of transect (km)
    :type width: float
    :param spacing: Discretization size (km) for RF ray sampling
    :type spacing: float
    :param max_depth: Maximum depth of slice below the transect line (km)
    :type max_depth: float
    :param channels: String of comma-separated component IDs to source for the RF amplitude
    :type channels: str, comma separated
    :param background_model: 1D background model to assume
    :type background_model: str
    :param stacked_scale: Max value to represent on color scale of CCP plot
    :type stacked_scale: float
    :param title: Title to place at top of CCP plot
    :type title: str
    :return: None
    """

    channels = channels.split(',')

    output_file_base, ext = os.path.splitext(output_file)
    if ext != ".png":
        output_file += ".png"
    matrix_norm, sample_density, length, stn_params = \
        ccp_generate(rf_stream, start_latlon, end_latlon, width=width, spacing=spacing, max_depth=max_depth,
                     channels=channels, v_background=background_model, station_map_file=output_file_base + '_MAP.png')

    if matrix_norm is not None:
        # Range of stacked amplitude for imshow to get best contrast
        if stacked_scale is not None:
            vlims = (-stacked_scale, stacked_scale)
        else:
            vlims = None
        # endif
        plot_ccp(matrix_norm, length, max_depth, spacing, ofile=output_file, vlims=vlims, metadata=stn_params,
                 title=title)
        if plot_density and sample_density is not None:
            sample_density_file = output_file_base + '_SAMPLE_DENSITY.png'
            # Use median of number of events per station to set the scale range.
            sc = sorted([s['event_count'] for s in stn_params.values() if s is not None])
            median_samples = sc[len(sc)//2]
            plot_ccp(sample_density, length, max_depth, spacing, ofile=sample_density_file, vlims=(0, median_samples),
                     metadata=stn_params, title=title + ' [sample density]' if title else None)
        # end if
    # end if
# end func


def run_batch(transect_file, rf_waveform_file, fed_db_file, stack_scale=0.4, width=30.0, spacing=2.0,
              max_depth=200.0, channels='R'):
    """Run CCP generation in batch mode along a series of transects.

    :param transect_file: File containing specification of network and station locations of ends of transects
    :type transect_file: str or Path
    :param rf_waveform_file: HDF5 file of QA'd receiver functions for the network matching the transect file
    :type rf_waveform_file: str or Path
    :param fed_db_file: Name of file with which to initialize FederatedASDFDataBase
    :type fed_db_file: str or Path
    :param stack_scale: Max value to represent on color scale of CCP plot
    :type stack_scale: float
    :param width: Width of transect (km)
    :type width: float
    :param spacing: Discretization size (km) for RF ray sampling
    :type spacing: float
    :param max_depth: Maximum depth of slice below the transect line (km)
    :type max_depth: float
    :param channels: String of comma-separated component IDs to source for the RF amplitude
    :type channels: str, comma separated
    :return: None
    """

    print("Reading HDF5 file...")
    rf_stream = rf.read_rf(rf_waveform_file, 'H5')

    db = FederatedASDFDataSet.FederatedASDFDataSet(fed_db_file)
    sta_coords = db.unique_coordinates

    with open(transect_file, 'r') as f:
        net = f.readline().strip()
        for transect in f.readlines():
            if not transect.strip():
                continue
            sta_start, sta_end = transect.split(',')
            sta_start = sta_start.strip()
            sta_end = sta_end.strip()
            start = '.'.join([net, sta_start])
            end = '.'.join([net, sta_end])
            start = np.array(sta_coords[start])
            end = np.array(sta_coords[end])
            # Offset ends slightly to make sure we don't lose end stations due to truncation error.
            # Note: for simplicity this treats lat/lon like cartesian coords, but this is approximate
            # and will break down near poles, for long transects, or if transect crosses the antimeridian.
            dirn = (end - start)
            dirn = dirn/np.linalg.norm(dirn)
            start -= 25*dirn/KM_PER_DEG
            end += 25*dirn/KM_PER_DEG
            start_latlon = (start[1], start[0])
            end_latlon = (end[1], end[0])

            outfile = '{}-ZRT-R_CCP_stack_{}-{}_{}km_spacing.png'.format(net, sta_start, sta_end, spacing)
            title = 'Network {} CCP R-stacking (profile {}-{})'.format(net, sta_start, sta_end)

            run(rf_stream, outfile, start_latlon, end_latlon, width, spacing, max_depth, channels,
                stacked_scale=stack_scale, title=title)

    # end for
    # end with
# end func

# ---------------- MAIN ----------------
@click.command()
@click.option('--start-latlon', nargs=2, type=float, required=True,
              help='Start coordinates of the profile line as latitude longitude (degrees),'
                   ' using space as separator, e.g. -22 133')
@click.option('--end-latlon', nargs=2, type=float, required=True,
              help='End coordinates of the profile line as latitude longitude (degrees),'
                   ' using space as separator, e.g. -19 140')
@click.option('--width', type=float, default=100.0, show_default=True,
              help='Width of the strip around the transect line (km)')
@click.option('--spacing', type=float, default=2.0, show_default=True,
              help='Spacing of cells in the 2D square mesh covering the slice below the transect line (km)')
@click.option('--max-depth', type=float, default=200.0, show_default=True,
              help='Maximum depth of slice below the transect line (km)')
@click.option('--stacked-scale', type=float,
              help='Max value to represent on color scale of CCP plot. Adjust for optimal contrast.')
@click.option('--channels', type=str, default='R', show_default=True,
              help='Comma separated list of channels to use for stacking, e.g. R,T')
@click.option('--background-model', type=click.Choice(['ak135', 'ak135_60']), default='ak135',
              show_default=True, help='1D background model to assume.')
@click.option('--title', type=str, help='Title text applied to the plots.')
@click.argument('rf-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output-file', type=click.Path(exists=False, dir_okay=False), required=True)
def main(rf_file, output_file, start_latlon, end_latlon, width, spacing, max_depth, channels,
         background_model, stacked_scale=None, title=None):
    # rf_file is the clean H5 file of ZRT receiver functions, generated by rf_quality_filter.py
    print("Reading HDF5 file...")
    stream = rf.read_rf(rf_file, 'H5')
    run(stream, output_file, start_latlon, end_latlon, width, spacing, max_depth, channels, background_model,
        stacked_scale, title)
# end func main


# ---------------- MAIN ----------------
if __name__ == "__main__":
    # run_batch('AQ_CCP_transects.txt', 'AQT_rfs_20151128T042000-20191108T000317_ZRT_it_rev1_qual.h5',
    #           '/g/data1a/ha3/Passive/SHARED_DATA/Index/asdf_files.txt', stack_scale=0.4, width=40.0,
    #           spacing=2.0, max_depth=100.0)
    # run_batch('7X_CCP_transects.txt', '7X_rfs_20090616T034200-20110401T231849_ZRT_it_rev2_qual.h5',
    #           '/g/data1a/ha3/Passive/SHARED_DATA/Index/asdf_files.txt', stack_scale=0.3, width=40.0,
    #           spacing=2.0, max_depth=100.0)
    main()  # pylint: disable=no-value-for-parameter
# end if
