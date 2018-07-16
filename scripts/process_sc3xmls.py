# infolder: Githubz/passive-seismic (master)
# conda install obspy
# conda install mpi4py
# export PYTHONPATH=.:$PYTHONPATH
# export ELLIPCORR=/Softlab/Githubz/passive-seismic/ellip-corr
# python2 scripts/process_sc3xmls.py /Softlab/Data/PST/sc3xmls/events_xmls_sc3ml/

import os
import sys
import logging
import glob
from obspy import read_events
from obspy.geodetics import locations2degrees, gps2dist_azimuth
from seismic import mpiops
import ellipcorr

from seismic.cluster.cluster import Grid
from inventory.parse_inventory import read_all_stations


log= logging.getLogger()

def process_sc3xml_files(path2dir, output_csv):
    """
    Process a directory of sc3 xml files to get events attributes
    Ref:/g/data/ha3/Passive/perm_stations/extracted_events/analyst-reviewed

    :param path2dir:
    :param output_csv:
    :return:
    """
    grid = Grid()

    stations = read_all_stations() # mpiops.run_once(read_all_stations)

    xml_file_pattern=os.path.join(path2dir, '*.xml')
    evt_files = glob.glob(xml_file_pattern)

    with open(output_csv, 'w') as csv_out:
        for evt_file in evt_files:
            print ("processing", evt_file)
            evts = read_events(evt_file)
            for evt in evts:
                print(str(evt))

                p_arr, s_arr, missing_stations, arr_stations=process_event(evt, stations, grid, 'P S', 1000)

                print (p_arr)
                print(s_arr)
                print(arr_stations)

                prefor = evt.preferred_origin() if evt.preferred_origin() else evt.origins[0]
                magPresent = True
                if evt.preferred_magnitude() or len(evt.magnitudes) > 0:
                    prefmag = evt.preferred_magnitude() or evt.magnitudes[0]
                else:
                    magPresent = False
                csv_out.write(str(prefor.latitude) + ',' + str(prefor.longitude) + ',' + str(prefor.depth/1000) + ',' + (str(prefmag.mag) if magPresent else '') + ',' + str(prefor.time) + '\n')

    return 0

def process_event(event, stations, grid, wave_type, counter):
    """
    :param event: obspy.core.event.Event class instance
    :param stations: dict
        stations dict
    :param grid: Grid class instance
    :param wave_type: str
        Wave type pair to generate inversion inputs. See `gather` function.
    :param counter: int
        event counter in this process
    """
    p_type, s_type = wave_type.split()

    # use preferred origin timestamp as the event number
    # if preferred origin is not populated, use the first origin timestamp
    origin = event.preferred_origin() or event.origins[0]

    # TODO: longer term uncomment the following line when fortran TT code is
    # updated to use longer integers
    # ev_number = int(origin.time.timestamp)

    # TODO: delete this definition of ev_number when fortran code can use
    # longer integers

    # the following definition ensures a 8 digit event number that is also
    # unique
    assert counter < 100000, 'Counter must be less than 100000'
    ev_number = int(str(counter) + '{0:0=3d}'.format(mpiops.rank))

    p_arrivals = []
    s_arrivals = []
    missing_stations = []
    arrival_staions = []

    # other event parameters we need
    ev_latitude = origin.latitude
    ev_longitude = origin.longitude
    ev_depth = origin.depth

    if ev_latitude is None or ev_longitude is None or ev_depth is None:
        return p_arrivals, s_arrivals, missing_stations, arrival_staions

    event_block = grid.find_block_number(ev_latitude, ev_longitude, z=ev_depth)

    for arr in origin.arrivals:
        sta_code = arr.pick_id.get_referred_object(
        ).waveform_id.station_code

        # ignore arrivals not in stations dict, workaround for now for
        # ENGDAHL/ISC events
        # TODO: remove this condition once all ISC/ENGDAHL stations are
        # available
        # Actually it does not hurt retaining this if condition. In case,
        # a station comes in which is not in the dict, the data prep will
        # still work
        # Note some stations are still missing even after taking into account
        #  of all seiscomp3 stations, ISC and ENGDAHL stations
        if sta_code not in stations:
            log.warning('Station {} not found in inventory'.format(sta_code))
            missing_stations.append(str(sta_code))
            continue
        sta = stations[sta_code]

        degrees_to_source = locations2degrees(ev_latitude, ev_longitude,
                                              sta.latitude,
                                              sta.longitude)

        # ignore stations more than 90 degrees from source
        if degrees_to_source > 90.0:
            # log.info('Ignored this station arrival as distance from source '
            #          'is {} degrees'.format(degrees_to_source))
            continue

        # TODO: use station.elevation information
        station_block = grid.find_block_number(sta.latitude, sta.longitude, z=0.0)

        if arr.phase in wave_type.split():

            ellipticity_corr = ellipcorr.ellipticity_corr(
                phase=arr.phase,
                edist=degrees_to_source,
                edepth=ev_depth / 1000.0,
                # TODO: check co-latitude definition
                # no `ecolat` bounds check in fortran ellipcorr subroutine
                # no `origin.latitude` bounds check in obspy
                ecolat=90 - ev_latitude,  # conversion to co-latitude
                azim=gps2dist_azimuth(ev_latitude, ev_longitude,
                                      sta.latitude, sta.longitude)[1]
            )
            t_list = [event_block, station_block, arr.time_residual,
                      ev_number, ev_longitude, ev_latitude, ev_depth,
                      sta.longitude, sta.latitude,
                      (arr.pick_id.get_referred_object().time.timestamp -
                       origin.time.timestamp) + ellipticity_corr,
                      degrees_to_source,
                      sta_code]
            arrival_staions.append(sta_code)
            p_arrivals.append(t_list + [1]) if arr.phase == p_type else \
                s_arrivals.append(t_list + [2])
        else:  # ignore the other phases
            pass
    return p_arrivals, s_arrivals, missing_stations, arrival_staions

# ======================================================
# python2 scripts/process_sc3xmls.py /Softlab/Data/PST/sc3xmls/events_xmls_sc3ml/
if __name__ == "__main__":
    if len(sys.argv)>1:
        in_dir=sys.argv[1]
    else:
        in_dir="../tests/mocks/events/analyst_event_samples/"

    process_sc3xml_files(in_dir, "myoutput.csv")
