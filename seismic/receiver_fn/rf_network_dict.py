#!/usr/bin/env python
# coding: utf-8
"""Class encapsulating a collection of receiver functions for stations of one network.
"""

from collections import defaultdict

# pylint: disable=invalid-name


class NetworkRFDict:
    """
    Collection of RFs for a given network indexed by station code, channel code.
    Note that location codes are not taken into account.
    """
    def __init__(self, rf_stream):
        """Initialize from rf.RFStream

        :param rf_stream: RFStream data
        :type rf_stream: rf.RFStream
        """
        self.db = defaultdict(lambda: defaultdict(list))
        self.network = None
        self.rotation = None
        self.phase = None
        for tr in rf_stream:
            net, sta, _, cha = tr.id.split('.')
            self.db[sta][cha].append(tr)
            if self.network is not None:
                assert net == self.network
            else:
                self.network = net
            # end if
            rotation = tr.stats.rotation
            phase = tr.stats.phase
            if self.rotation is not None: assert rotation == self.rotation
            else: self.rotation = rotation

            if self.phase is not None: assert phase == self.phase
            else: self.phase = phase
        # end for
    # end func

    def __iter__(self):
        return iter(self.db.items())
    # end func

    def __getitem__(self, key):
        return self.db[key]
    # end func

    def __len__(self):
        return len(self.db)
    # end func

    def keys(self):
        """Accessor for the top level keys (station codes) of the network in
        an iterable container.

        :return: Iterable of top level keys to the dictionary
        :rtype: Python iterable
        """
        return self.db.keys()
    # end func

# end class
