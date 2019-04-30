#!/usr/bin/env python
"""
Generate common conversion point (CCP) plot as per *C.Sippl, "Moho geometry along a north-south passive seismic
transect through Central Australia", Technophysics 676 (2016), pp.56-69, 
DOI https://doi.org/10.1016/j.tecto.2016.03.031*

This code adapted from Christian Sippl's original code.

Workflow:
    prepare_rf_data.py --> generate_rf.py --> rf_smart_bin.py --> plot_ccp.py (this script)
"""

# Remove this after initial development
# pylint: skip-file

import os
import sys
import numpy as np

import click
import obspy
from obspy.taup import TauPyModel
import matplotlib.pyplot as plt
from mpl_toolkits import basemap
import rf

from rf_util import KM_PER_DEG
from tqdm import tqdm

PY2 = sys.version_info[0] < 3
if PY2:
    import cPickle as pkl
else:
    import pickle as pkl

USE_PICKLE = False

def plot_ccp(matrx, length, max_depth, spacing, ofile=None, vlims=None):
    """
    plot results of CCP stacking procedure
    """
    tickstep_x = 50
    tickstep_y = 25

    plt.figure(figsize=(16,9))
    if vlims is not None:
        plt.imshow(matrx, aspect='equal', cmap='jet', vmin=vlims[0], vmax=vlims[1])
    else:
        plt.imshow(matrx, aspect='equal', cmap='jet')

    plt.ylim([int(max_depth/spacing), 0])

    plt.xlabel('distance [km]')
    plt.ylabel('depth [km]') 

    plt.xticks(range(0, int(length/spacing), int(tickstep_x/spacing)), range(0, int(length), tickstep_x))
    plt.yticks(range(0, int(max_depth/spacing), int(tickstep_y/spacing)), range(0, int(max_depth), tickstep_y))

    if ofile:
        plt.savefig(ofile, dpi=300)
    else:
        plt.show()
    plt.close()


def setup_ccp_profile(length, spacing, maxdep):
    """
    construct the grid for a CCP stacking profile
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


def get_amplitude(trace, time, rf_offset=5.0):
    """
    retrieve amplitude value
    """
    indx = (time + rf_offset) * trace.stats.sampling_rate
    amp = trace.data[int(indx)]/trace.stats.amax
    return amp


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


def matrx_lookup(xsz, sta_offset, h, depstep, lenstep):
    """
    return index values for amplitude contrbution in profile matrix
    """
    distance_offset = sta_offset - xsz # because zero is in the north

    diff_x = 999.0
    diff_y = 999.0
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


def bounding_box(startpoint, endpoint):
    ybig = max(startpoint[0], endpoint[0])
    ysmall = min(startpoint[0], endpoint[0])
    xbig = max(startpoint[1], endpoint[1])
    xsmall = min(startpoint[1], endpoint[1])
    return (xsmall, ysmall, xbig, ybig)


def equirectangular_projection(x0, y0, x1, y1):
    # This length calculation appears to use the forward equirectangular projection (https://en.wikipedia.org/wiki/Equirectangular_projection)
    # (See also https://www.movable-type.co.uk/scripts/latlong.html)
    profile_x_length = (x1 - x0) * KM_PER_DEG * np.cos((y1 + y0) / 2. * np.pi / 180.)
    profile_y_length = (y1 - y0) * KM_PER_DEG
    # This is an approximate great circle arc length
    length = np.sqrt(profile_x_length**2 + profile_y_length**2)

    return profile_x_length, profile_y_length, length


def ccp_compute_station_params(rf_stream, startpoint, endpoint, width, bm):

    xsmall, ysmall, xbig, ybig = bounding_box(startpoint, endpoint)
    profile_x_length, profile_y_length, length = equirectangular_projection(xsmall, ysmall, xbig, ybig)

    # TODO: Reverse engineer the undocumented coordinate system here. Maybe zero azimuth is cartographic north? Clockwise like a compass bearing?
    try:
        # Why not using atan2 here? Then we wouldn't need to "find quadrant" block below.
        profile_azimuth = (np.arctan(profile_y_length / (float(profile_x_length))) * 180.0) / np.pi
    except ZeroDivisionError:
        profile_azimuth = 0.0
    #find quadrant (additive term to angle...)
    if startpoint[0] >= endpoint[0] and startpoint[1] < endpoint[1]:
        add = 90.
    elif startpoint[0] < endpoint[0] and startpoint[1] <= endpoint[1]:
        add = 0.
    elif startpoint[0] > endpoint[0] and startpoint[1] >= endpoint[1]:
        add = 180.
    elif startpoint[0] <= endpoint[0] and startpoint[1] > endpoint[1]:
        add = 270.
    profile_azimuth += add
    # profile_azimuth = np.arctan2(profile_y_length, profile_x_length) * 180.0 / np.pi

    stn_params = {}
    angle_norm = profile_azimuth % 90
    xstart = 0
    ystart = 0
    xend, yend, _ = equirectangular_projection(startpoint[1], startpoint[0], endpoint[1], endpoint[0])
    pbar = tqdm(total=len(rf_stream), ascii=True)
    for tr in rf_stream:
        pbar.update()
        stat_code = tr.stats.station
        pbar.set_description("Station {}".format(stat_code))
        if stat_code not in stn_params:
            # TODO: Replace obfuscated trigonometric calculations here with vector geometry and projection calculations.
            #       (eliminate the trig and improve code transparency)
            sta_lat = tr.stats.station_latitude
            sta_lon = tr.stats.station_longitude
            xstat, ystat, _ = equirectangular_projection(startpoint[1], startpoint[0], sta_lon, sta_lat)
            # Assumption to make sense of undocumented code: "dist" here is the orthogonal Euclidean distance of the station from the profile line.
            dist = ((yend - ystart) * xstat - (xend - xstart) * ystat + xend*ystart - yend*xstart) / np.sqrt((yend - ystart)**2 + (xend - xstart)**2)
            within_dist_tolerance = (abs(dist) <= width)
            within_profile_length = False
            if within_dist_tolerance:
                if abs(profile_azimuth - 90.) > 30. and abs(profile_azimuth - 270.) > 30.:
                    #calculate position on profile (length)
                    sta_offset_n = ((startpoint[0] - sta_lat) * KM_PER_DEG) / np.sin(angle_norm * np.pi/180.)
                    n_correction = dist / np.tan(angle_norm*np.pi/180.)
                    sta_offset = sta_offset_n + n_correction
                else:
                    sta_offset_e = -((startpoint[1] - sta_lon) * KM_PER_DEG * np.cos(startpoint[0] * np.pi / 180.)) / np.sin((90. - angle_norm) * np.pi / 180.)
                    e_correction = dist / np.tan((90 - angle_norm)*np.pi/180.)
                    sta_offset = sta_offset_e - e_correction
                # end if
                within_profile_length = (sta_offset > 0 and sta_offset < length)
            # end if
            
            if within_dist_tolerance and within_profile_length:
                pbar.write("Station " + stat_code + " included: (offset, dist) = ({}, {})".format(sta_offset, dist))
                #add station to CCP stack
                x, y = bm(sta_lon, sta_lat)
                bm.plot(x, y, 'ro')
                stn_params[stat_code] = {'dist': dist, 'sta_offset': sta_offset}
            else:
                x, y = bm(sta_lon, sta_lat)
                bm.plot(x, y, 'k^')
                stn_params[stat_code] = None
            # end if
        # end if
    # end for
    pbar.close()

    return stn_params


def ccp_generate(rf_stream, startpoint, endpoint, width, spacing, max_depth, v_background='ak135'):
    # rf_stream should be of type rf.rfstream.RFStream

    xsmall, ysmall, xbig, ybig = bounding_box(startpoint, endpoint)
    _, _, length = equirectangular_projection(xsmall, ysmall, xbig, ybig)

    #set velocity model (other options can be added) 
    if v_background == 'ak135':
        b_lay = [0., 20., 35., 77.5, 120., 165., 210., 260.,
                 310., 360., 410., 460., 510., 560., 610., 660.]
        vp = [5.8, 6.5, 8.04, 8.045, 8.05, 8.175, 8.3, 8.4825, 8.665,
              8.8475, 9.03, 9.36, 9.528, 9.696, 9.864, 10.032, 10.2]
        vs = [3.46, 3.85, 4.48, 4.49, 4.5, 4.509, 4.518, 4.609, 4.696,
              4.783, 4.87, 5.08, 5.186, 5.292, 5.398, 5.504, 5.61]
        vmod = (b_lay, vp, vs)

    #first set up a matrix
    profile_mesh, depstep, lenstep = setup_ccp_profile(length, spacing, max_depth)

    # snapshot of initial matrix
    mesh_entries = profile_mesh.copy()

    #create map plot file
    m = basemap.Basemap(projection='merc', urcrnrlat=ybig + 1.0, urcrnrlon=xbig + 1.0,
                        llcrnrlon=xsmall - 1.0, llcrnrlat=ysmall - 1.0, resolution='i')
    m.drawcoastlines()
    x1, y1 = m(startpoint[1], startpoint[0])
    x2, y2 = m(endpoint[1], endpoint[0])
    m.plot([x1, x2], [y1, y2], 'r--')

    # Precompute the station parameters for a given code, as this is the same for every trace.
    print("Computing included stations...")
    stn_params = ccp_compute_station_params(rf_stream, startpoint, endpoint, width, m)

    # Processing/extraction of rf_stream data
    print("Projecting included stations to slice...")
    model = TauPyModel(model=v_background)
    pbar = tqdm(total=len(rf_stream), ascii=True)
    # TARGET_STNS = ['BS24', 'BS25', 'BS26', 'BS27', 'BS28']
    for tr in rf_stream:
        pbar.update()
        stat_code = tr.stats.station
        # if stat_code not in TARGET_STNS:
        #     continue
        pbar.set_description("{} event {}".format(stat_code, tr.stats.event_id))

        if stn_params[stat_code] is not None:
            # Updating of matrix
            sta_offset = stn_params[stat_code]['sta_offset']
            try:
                back_azi = tr.stats.back_azimuth
                arr_p = model.get_travel_times(tr.stats.event_depth, tr.stats.distance, phase_list=['P'])
                inc_p = arr_p[0].incident_angle
                profile_mesh, mesh_entries = add_ccp_trace(tr, inc_p, profile_mesh, mesh_entries, vmod,
                                                           depstep, lenstep, sta_offset, back_azi)
            except IndexError as err:
                print(err)
                pass
        # end if
    # end for
    pbar.close()

    # Show basemap plot
    plt.show()
    plt.close()

    if not np.all(mesh_entries[:] == 0):
        #normalize, then plot
        matrx_norm = (profile_mesh / mesh_entries).transpose()
        return matrx_norm, length
    else:
        return None, 0


# ---------------- MAIN ----------------
if __name__ == "__main__":

    # rf_file is the clean H5 file of ZRT receiver functions, generated by rf_smart_bin.py
    rf_file = 'seismic/receiver_fn/DATA/OA-ZRT-R-cleaned.h5'

    # start_latlon = (-21.5, 133.0)
    # end_latlon = (-18.5, 139.0)
    start_latlon = (-22.0, 133.0)
    end_latlon = (-19.0, 133.0)

    width = 50.0
    spacing = 1.0
    max_depth = 200.0
    vmin, vmax = (-0.10, 0.10)

    rf_file_base, _ = os.path.splitext(rf_file)
    pkl_file = rf_file_base + '.pkl'
    if os.path.exists(pkl_file) and USE_PICKLE:
        with open(pkl_file, 'rb') as f:
            matrix_norm, length = pkl.load(f)
    else:
        print("Reading HDF5 file...")
        stream = rf.read_rf(rf_file, 'H5')
        matrix_norm, length = ccp_generate(stream, start_latlon, end_latlon, width=width, spacing=spacing, max_depth=max_depth)
        with open(pkl_file, 'wb') as f:
            pkl.dump((matrix_norm, length), f, pkl.HIGHEST_PROTOCOL)

    if matrix_norm is not None:
        outfile = rf_file_base + '.png'
        plot_ccp(matrix_norm, length, max_depth, spacing, ofile=outfile, vlims=(vmin, vmax))
