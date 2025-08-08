"""
Description:
    Wrapper Class for providing fast access to data contained within a set of ASDF files

References:

CreationDate:   12/12/18
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     12/12/18   RH
    LastUpdate:     2020-04-10 Fei Zhang  clean up + added example run for the script
"""

from collections import defaultdict

# from mpi4py import MPI
import numpy as np
from scipy.spatial import cKDTree

from seismic.ASDFdatabase._FederatedASDFDataSetImpl import _FederatedASDFDataSetImpl
from seismic.misc import rtp2xyz, setup_logger
from obspy.core import UTCDateTime
import click

class FederatedASDFDataSet():
    def __init__(self, asdf_source, force_reindex=False, logger=None,
                 single_item_read_limit_in_mb=1024,
                 single_threaded_access=True):
        """
        Initializer for FederatedASDFDataSet.

        :param asdf_source: Path to a text file containing a list of ASDF files. \
               Entries can be commented out with '#'
        :param force_reindex: Force reindex even if a preexisting db file is found
        :param logger: logger instance
        :param single_item_read_limit_in_mb: buffer size for Obspy reads
        :param single_threaded_access: By default, data are read via unthreaded MPI-processes.
               This can be relaxed for threaded GUI applications, though data access will still
               remain single-threaded.
        """
        self.logger = logger
        self.asdf_source = asdf_source
        self._earth_radius = 6371  # km

        # Instantiate implementation class
        self.fds = _FederatedASDFDataSetImpl(asdf_source, force_reindex=force_reindex, logger=logger,
                                             single_item_read_limit_in_mb=single_item_read_limit_in_mb,
                                             single_threaded_access=single_threaded_access)

        # Populate coordinates
        rtps_dict = defaultdict()
        for ds_dict in self.fds.asdf_station_coordinates:
            for key in list(ds_dict.keys()):

                lon, lat, _ = ds_dict[key]
                rtps_dict[key] = [self._earth_radius,
                                  np.radians(90 - lat),
                                  np.radians(lon)]
            # end for
        # end for

        rtps_list = []
        for k in list(rtps_dict.keys()):
            rtps_list.append(rtps_dict[k])
        # end for
        rtps = np.array(rtps_list)
        xyzs = rtp2xyz(rtps[:, 0], rtps[:, 1], rtps[:, 2])

        self._tree = cKDTree(xyzs)
        self._key_list = np.array(list(rtps_dict.keys()))
    # end func

    @property
    def unique_coordinates(self):
        """

        :return: dictionary containing [lon, lat] coordinates indexed by 'net.sta'
        """
        return self.fds._unique_coordinates
    # end func

    def corrections_enabled(self):
        """

        @return: whether GPS clock-corrections have been enabled by setting
                 the environment variable GPS_CLOCK_CORRECTION=1
        """
        return self.fds.corrections_enabled
    # end func

    def get_closest_stations(self, lon, lat, nn=1):
        """

        :param lon: longitude (degree)
        :param lat: latitude (degrees)
        :param nn: number of closest stations to fetch
        :return: A tuple containing a list of closest 'network.station' names and a list of distances
                 (in ascending order) in kms
        """
        assert nn > 0, 'nn must be > 0'

        xyz = rtp2xyz(np.array([self._earth_radius]),
                      np.array([np.radians(90 - lat)]),
                      np.array([np.radians(lon)]))
        d, l = self._tree.query(xyz, nn)

        if isinstance(l, int):
            l = [l]

        if (len(d.shape) == 1):
            d = np.expand_dims(d, axis=0)

        l = l[l < len(self.unique_coordinates)]

        if isinstance(l, int):
            l = [l]

        return (list(self._key_list[l]), d[0, :len(l)])

    # end func

    def get_recording_timespan(self, network, station=None, location=None, channel=None):
        """
        :param network: network code
        :param station: station code
        :param location: location code (optional)
        :param channel: channel code (optional)
        :return: tuple containing min and max times as UTCDateTime objects. If no matching records are found
                 min is set to 2100-01-01T00:00:00.000000Z and max is set to 1900-01-01T00:00:00.000000Z
        """

        return self.fds.get_recording_timespan(network, station=station, location=location, channel=channel)
    # end func

    def get_all_recording_timespans(self):
        """
        Get a structured numpy array with named columns
        'net', 'sta', 'loc', 'cha', 'min_st', 'max_et'
        representing contents of the database
        @return:
        """

        results = self.fds.get_all_recording_timespans()
        return results
    # end if

    def get_stations(self, starttime, endtime, network=None, station=None, location=None, channel=None):
        """
        :param starttime: start time string in UTCDateTime format; can also be an instance of obspy.UTCDateTime
        :param endtime: end time string in UTCDateTime format; can also be an instance of obspy.UTCDateTime
        :param network: network code (optional)
        :param station: station code (optional)
        :param location: location code (optional)
        :param channel: channel code (optional)

        :return: a list containing [net, sta, loc, cha, lon, lat, elev_m] in each row
        """
        results = self.fds.get_stations(starttime, endtime, network, station, location, channel)
        return results

    # end func

    def get_waveform_count(self, network, station, location, channel, starttime,
                           endtime):
        """
        Count the number of traces within the given parameters of network, station, etc..
        and date range. This is a fast method of determing whether any trace data exists
        in a given time period, if you don't actually need the waveform data itself.

        :param network: network code
        :param station: station code
        :param location: location code
        :param channel: channel code
        :param starttime: start time string in UTCDateTime format; can also be an instance of obspy.UTCDateTime
        :param endtime: end time string in UTCDateTime format; can also be an instance of obspy.UTCDateTime
        :return: The number of streams containing waveform data over the time-range provided
        """
        return self.fds.get_waveform_count(network, station, location, channel,
                                           starttime, endtime)

    # end func

    def get_waveforms(self, network, station, location, channel, starttime,
                      endtime, trace_count_threshold=200, nearest_sample=True):
        """
        :param network: network code
        :param station: station code
        :param location: location code
        :param channel: channel code
        :param starttime: start time string in UTCDateTime format; can also be an instance of obspy.UTCDateTime
        :param endtime: end time string in UTCDateTime format; can also be an instance of obspy.UTCDateTime
        :param trace_count_threshold: returns an empty Stream if the number of traces within the time-range provided
                                      exceeds the threshold (default 200). This is particularly useful for filtering
                                      out data from bad stations, e.g. those from the AU.Schools network
        :param nearest_sample: fetches nearest sample if True
        :return: an obspy.Stream containing waveform data over the time-rage provided
        """
        s = self.fds.get_waveforms(network, station, location, channel, starttime,
                                   endtime, trace_count_threshold, nearest_sample)
        return s

    # end func

    def get_location_codes(self, network, station, starttime=None, endtime=None):
        """
        :param network: network code
        :param station: station code
        :param starttime: start time string in UTCDateTime format; can also be an instance of obspy.UTCDateTime
        :param endtime: end time string in UTCDateTime format; can also be an instance of obspy.UTCDateTime

        :return: a list containing unique location codes within the timeframe specified
        """

        return self.fds.get_location_codes(network, station, starttime=starttime, endtime=endtime)

    # end func

    def stations_iterator(self, network_list=[], station_list=[]):
        """
        This function provides an iterator over the entire data volume contained in all the ASDF files listed in the
        text file during instantiation. When FederatedASDFDataSet is instantiated in an MPI-parallel environment,
        meta-data for the entire data volume are equally partitioned over all processors -- in such instances, this
        function provides an iterator over the data allocated to a given processor. This functionality underpins
        parallel operations, e.g. picking arrivals.

        :param network_list: a list of networks to process
        :param station_list: a list of stations to process

        :return: tuples containing [net, sta, start_time, end_time]; start- and end-times are instances of obspy.UTCDateTime
        """
        for item in self.fds.stations_iterator(network_list=network_list, station_list=station_list):
            yield item
        # end for
    # end func

    def get_inventory(self, network=None, station=None):
        """
        This function returns the combined (for all underlying ASDF files) xml inventory when both 'network' and 'station'
        are set to None, otherwise a subset is returned. Some processing workflows (e.g. RF) require an Obspy inventory
        to iterate over data -- this function is intended to cater for those requirements, while the more comprehensive
        'get_stations' function should be used for fetching matching stations that have waveform data within a given
        time interval

        :param network: network code
        :param station: station code
        """

        inv = self.fds.get_inventory(network=network, station=station)
        return inv
    # end func

    def find_gaps(self, network=None, station=None, location=None,
                  channel=None, start_date_ts=None, end_date_ts=None,
                  min_gap_length=86400):
        """
        This function returns gaps in data as a numpy array with columns: net, sta, loc, cha, start_timestamp,
        end_timestamp.
        @param network: network code
        @param station: station code
        @param location: location code
        @param channel: channel code
        @param start_date_ts: start timestamp
        @param end_date_ts: end timestamp
        @param min_gap_length: minimum length of gap in seconds; smaller gaps in data are ignored
        @return:
        """
        return self.fds.find_gaps(network, station, location, channel, start_date_ts, end_date_ts, min_gap_length)
    # end func

    def get_recording_duration(self, network=None, station=None, location=None, channel=None,
                                     starttime=None, endtime=None, cumulative=False):
        """
        Fetches total recording duration in seconds. Note that 'duration_seconds' in the output exclude data-gaps

        @param network:
        @param station:
        @param location:
        @param channel:
        @param starttime:
        @param endtime:
        @param cumulative: returns cumulative recording times, otherwise blocks of start- and end-times
        @return: Numpy record array with columns, if cumulative=False:
                 net, sta, loc, cha, block_st, block_et
                 , otherwise:
                 net, sta, loc, cha, lon, lat, min_st, max_et, duration_seconds
        """

        rows = self.fds.get_recording_duration(network=network, station=station, location=location, channel=channel,
                                               starttime=starttime, endtime=endtime, cumulative=cumulative)
        return rows
    # end func
# end class


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('asdf-source', required=True,
                type=click.Path(exists=True))
@click.option('--force-reindex', default=False, is_flag=True,
              help='Force reindex, even if a preexisting database is found')
@click.option('--generate-summary', default=False, is_flag=True,
              help='Generate coverage and data availability summaries')
def process(asdf_source, force_reindex, generate_summary):
    """
    ASDF_SOURCE: Text file containing a list of paths to ASDF files
    """

    ofn = 'FederatedASDFDataSet.Indexer.log'
    logger = setup_logger('', ofn)
    ds = FederatedASDFDataSet(asdf_source, force_reindex=force_reindex, logger=logger)

    if(generate_summary):
        if(ds.fds.rank == 0):
            ts = UTCDateTime().strftime("%Y-%m-%dT%H.%M.%S")
            logger.info('Generating coverage summary..')
            ofn = 'FederatedASDFDataSet.Summary.{}.txt'.format(ts)

            with open(ofn, 'w') as fh:
                fh.write('# net, sta, loc, cha, lon, lat, min_starttime, max_endtime, duration_months\n')

                rows = ds.get_recording_duration(cumulative=True)
                for row in rows:
                    net, sta, loc, cha, min_st, max_et, duration_seconds = row
                    duration_months = duration_seconds / (86400 * 30)

                    lon, lat = ds.unique_coordinates['{}.{}'.format(net, sta)]
                    line = '{},{},{},{},{:3.4f},{:3.4f},{},{},{:5.3f}\n'.\
                           format(net, sta, loc, cha, lon, lat,
                                  UTCDateTime(min_st).strftime('%Y-%m-%dT%H:%M:%S'),
                                  UTCDateTime(max_et).strftime('%Y-%m-%dT%H:%M:%S'),
                                  duration_months)
                    
                    if(duration_seconds > (max_et - min_st)): 
                        logger.warn('Potential overlapping data found: {}'.format(line.strip()))
                    # end if
                    
                    fh.write(line)
                # end for
            # end with
        # end if
    # end if
    logger.info('Done..')
# end func

if __name__ == "__main__":
    process()
# end func
