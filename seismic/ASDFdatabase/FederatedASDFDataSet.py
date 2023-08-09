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
from seismic.misc import rtp2xyz

class FederatedASDFDataSet():
    def __init__(self, asdf_source, logger=None,
                 single_item_read_limit_in_mb=1024,
                 single_threaded_access=True):
        """
        Initializer for FederatedASDFDataSet.

        :param asdf_source: Path to a text file containing a list of ASDF files. \
               Entries can be commented out with '#'
        :param logger: logger instance
        :param single_item_read_limit_in_mb: buffer size for Obspy reads
        :param single_threaded_access: By default, data are read via unthreaded MPI-processes.
               This can be relaxed for threaded GUI applications, though data access will still
               remain single-threaded.
        """
        self.logger = logger
        self.asdf_source = asdf_source
        self._unique_coordinates = None
        self._earth_radius = 6371  # km

        # Instantiate implementation class
        self.fds = _FederatedASDFDataSetImpl(asdf_source, logger=logger,
                                             single_item_read_limit_in_mb=single_item_read_limit_in_mb,
                                             single_threaded_access=single_threaded_access)

        # Populate coordinates
        self._unique_coordinates = defaultdict(list)

        rtps_dict = defaultdict()
        for ds_dict in self.fds.asdf_station_coordinates:
            for key in list(ds_dict.keys()):
                self._unique_coordinates[key] = [ds_dict[key][0], ds_dict[key][1]]

                rtps_dict[key] = [self._earth_radius,
                                  np.radians(90 - ds_dict[key][1]),
                                  np.radians(ds_dict[key][0])]
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
        return self._unique_coordinates

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

    def get_global_time_range(self, network, station, location=None, channel=None):
        """
        :param network: network code
        :param station: station code
        :param location: location code (optional)
        :param channel: channel code (optional)
        :return: tuple containing min and max times as UTCDateTime objects. If no matching records are found
                 min is set to 2100-01-01T00:00:00.000000Z and max is set to 1900-01-01T00:00:00.000000Z
        """

        return self.fds.get_global_time_range(network, station, location=location, channel=channel)

    # end func

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
                      endtime, trace_count_threshold=200):
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
        :return: an obspy.Stream containing waveform data over the time-rage provided
        """
        s = self.fds.get_waveforms(network, station, location, channel, starttime,
                                   endtime, trace_count_threshold)
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
        @param min_gap_length: minimum length of gap; smaller gaps in data are ignored
        @return:
        """
        return self.fds.find_gaps(network, station, location, channel, start_date_ts, end_date_ts, min_gap_length)
    # end func
# end class

if __name__ == "__main__":
    """
    How to Run Example::

        python ASDFdatabase/FederatedASDFDataSet.py /Datasets/asdf_file_index.txt

    Upon success, a db file will be created: /Datasets/f374ca9e7dd8abd2a1d58575e0d55520f30ffc23.db
    """
    import sys
    from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet

    if len(sys.argv) < 2:
        print("******** USAGE: python3 %s %s **********"% (sys.argv[0], "asdf_file_list_txt"))
        sys.exit(1)

    asdf_file_list = sys.argv[1]
    ds = FederatedASDFDataSet(asdf_file_list)