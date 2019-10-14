#!/usr/bin/env python
# coding: utf-8
"""Class encapsulating a collection of receiver functions for stations of one network.
"""

from collections import defaultdict

# pylint: disable=invalid-name


class NetworkRFDict:
    """Collection of RFs for a given network indexed by station code, channel code.
    """
    def __init__(self, rf_stream):
        """Initialize from rf.RFStream

        :param rf_stream: RFStream data
        :type rf_stream: rf.RFStream
        """
        self.db = defaultdict(lambda: defaultdict(list))
        self.network = None
        self.rotation = None
        for s in rf_stream:
            net, sta, _, cha = s.id.split('.')
            self.db[sta][cha].append(s)
            if self.network:
                assert net == self.network
            else:
                self.network = net
            # end if
            rotation = s.stats.rotation
            if self.rotation:
                assert rotation == self.rotation
            else:
                self.rotation = rotation
            # end if
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

# end class
