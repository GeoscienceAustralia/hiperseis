#!/usr/bin/env python
"""Helper class to iterate over station events in h5 file generated by rf library,
   without loading all the traces into memory. This is a scalable solution for very large files.
"""

import logging

import h5py
from obspyh5 import dataset2trace
from rf import RFStream

logging.basicConfig()

# pylint: disable=invalid-name, logging-format-interpolation

class IterRfH5StationEvents(object):
    """Helper class to iterate over stations in h5 file generated by extract_event_traces.py and pass
       them to RF generator. This class avoids having to load the whole file up front via obspy which
       is slow and not scalable.

       This class yields 3-channel traces per station per event.

       Data yielded per station can easily be many MB in size.
    """

    def __init__(self, h5_filename, memmap=False):
        """Initializer

        :param h5_filename: Name of file containing event seismograms in HDF5 format, indexed in \
            seismic.stream_io.EVENTIO_H5INDEX format
        :type h5_filename: str or pathlib.Path
        :param memmap: If True, memmap the open file. Can improve tractability of handling very large files.
        :type memmap: bool
        """
        self.h5_filename = h5_filename
        # self.num_components = 3
        self.memmap_input = memmap

    def _open_source_file(self):
        if self.memmap_input:
            try:
                return h5py.File(self.h5_filename, 'r', driver='core', backing_store=False)
            except OSError as e:
                logger = logging.getLogger(__name__)
                logger.error("Failure to memmap input file with error:\n{}\nReverting to default driver."
                             .format(str(e)))
        # end if
        return h5py.File(self.h5_filename, 'r')
        # end if

    def __iter__(self):
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.INFO)
        logger.info("Scanning station groups from file {}".format(self.h5_filename))
        with self._open_source_file() as f:
            wf_data = f['waveforms']
            num_stations = len(wf_data)
            count = 0
            event_count = 0

            for station_id in wf_data:
                count += 1
                logger.info("Station {} {}/{}".format(station_id, count, num_stations))
                station_data = wf_data[station_id]
                station_stream3c = []
                for event_time in station_data:
                    event_traces = station_data[event_time]
                    # if len(event_traces) != self.num_components:
                    #     logging.warning("Incorrect number of traces ({}) for stn {} event {}, skipping"
                    #                     .format(len(event_traces), station_id, event_time))
                    #     continue
                    traces = []

                    for trace_id in event_traces:
                        trace = dataset2trace(event_traces[trace_id])
                        traces.append(trace)

                    event_count += 1
                    station_stream3c.append(RFStream(traces=traces).sort())
                # end for

                # Yield the results with 3-channel trace triplets grouped together in RFStream instances.
                yield station_id, station_stream3c
            # end for
        # end with
        logger.info("Yielded {} event traces to process".format(event_count))
    # end func
