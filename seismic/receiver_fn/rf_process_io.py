#!/usr/bin/env python
"""Helper for asynchronous writing of RF data to file.
"""

import logging

import h5py
from rf import RFStream

logging.basicConfig()

# pylint: disable=invalid-name, logging-format-interpolation


def async_write(rfstream_queue, outfile_name, max_buffered=100, metadata=''):
    """Monitors asynchronous queue for data, removes from queue to buffer, then
       flushes buffer intermittently and when queue termination signal is put.

       When None is received on the queue, this is taken as the signal to terminate
       monitoring the queue.

    :param rfstream_queue: Queue into which RFStreams are pushed for writing to file.
    :type rfstream_queue: multiprocessing.Manager.Queue
    :param outfile_name: Name of file into which queued RFStream results are periodically written.
    :type outfile_name: str or Path
    :param max_buffered: Maximum number of RFStreams to buffer before flushing to file, defaults to 100
    :type max_buffered: int, optional
    :param metadata: Metadata string to write to root attribute
    :type metadata: str, optional
    """
    buffered_streams = []
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    logger.info("Starting async write thread")
    with h5py.File(outfile_name, mode='w') as h5f:
        h5f.attrs['metadata'] = metadata
    # end with
    while True:
        # Passing None into the queue is taken as signal to flush buffer and terminate thread
        rfstream = rfstream_queue.get()

        terminating = (rfstream is None)
        if not terminating:
            buffered_streams.append(rfstream)
        else:
            logger.info("Flushing result buffer...")

        if len(buffered_streams) >= max_buffered or terminating:
            stream = RFStream()
            for rf in buffered_streams:
                stream.extend(rf)
            stream.write(outfile_name, 'H5', mode='a')
            logger.info("Flushed {} streams to output file {}".format(len(buffered_streams), outfile_name))

            while buffered_streams:
                buffered_streams.pop()
                rfstream_queue.task_done()

        if terminating:
            rfstream_queue.task_done()
            break

    logger.info("Terminating async write thread")
# end func
