#!/usr/bin/env python

import logging

from rf import RFStream

logging.basicConfig()

# pylint: disable=invalid-name, logging-format-interpolation

def async_write(rfstream_queue, outfile_name, max_buffered=100):
    """Monitors asynchronous queue for data, removes from queue to buffer, then
       flushes buffer intermittently and when queue termination signal is put.

       When None is received on the queue, this is taken as the signal to terminate
       monitoring the queue.

    :param rfstream_queue: [description]
    :type rfstream_queue: multiprocessing.Manager.Queue
    :param outfile_name: [description]
    :type outfile_name: str or Path
    :param max_buffered: [description], defaults to 100
    :type max_buffered: int, optional
    """
    buffered_streams = []
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    logger.info("Starting async write thread")
    first_write = True
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
            if first_write:
                mode = 'w'
                first_write = False
            else:
                mode = 'a'
            stream.write(outfile_name, 'H5', mode=mode)
            logger.info("Flushed {} streams to output file {}".format(len(buffered_streams), outfile_name))

            while buffered_streams:
                buffered_streams.pop()
                rfstream_queue.task_done()

        if terminating:
            rfstream_queue.task_done()
            break

    logger.info("Terminating async write thread")
# end func
