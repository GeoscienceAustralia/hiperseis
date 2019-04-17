#!/usr/bin/env python

# import os

from seismic.receiver_fn.hk_stacking import RF_c_sippl
# HK stacking code useless until source code for `hk2d` is found.
#from . import HK_stacking

RF_sippl_symbols = dir(RF_c_sippl)
print("{} symbols found in {}".format(len(RF_sippl_symbols), RF_c_sippl.__file__))
