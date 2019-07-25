#!/usr/bin/env python
"""Common format for Pandas Dataframe column ordering of Network, Station and Channel metadata.
"""

from collections import OrderedDict

import numpy as np
import pandas as pd

# pylint: disable=invalid-name

TABLE_SCHEMA = OrderedDict((
    ('NetworkCode', str),
    ('StationCode', str),
    ('Latitude', np.float64),
    ('Longitude', np.float64),
    ('Elevation', np.float64),
    ('StationStart', np.datetime64),
    ('StationEnd', np.datetime64),
    ('ChannelCode', str),
    ('ChannelStart', np.datetime64),
    ('ChannelEnd', np.datetime64)))

TABLE_COLUMNS = TABLE_SCHEMA.keys()

# See https://pandas.pydata.org/pandas-docs/stable/timeseries.html#timestamp-limitations
PANDAS_MAX_TIMESTAMP = str(pd.Timestamp.max)[0:19]

DEFAULT_START_TIMESTAMP = pd.Timestamp("1964-01-01 00:00:00")
# Subtract 1 hr as a workaround for pd.Timestamp.max still being treated as NaT in some versions of Pandas
DEFAULT_END_TIMESTAMP = pd.Timestamp.max - pd.Timedelta(1.0, 'h')
