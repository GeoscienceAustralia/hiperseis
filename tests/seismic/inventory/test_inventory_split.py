#!/usr/bin/env python

import os

from seismic.inventory.inventory_split import inventory_split


def test_inventory_split_main(tmp_path):
    """High level test of main entry point for inventory split.
    """
    self_path = os.path.dirname(os.path.abspath(__file__))
    test_data_path = os.path.join(self_path, 'data')
    stxml_file = os.path.join(test_data_path, 'test_inventory_merge_20190528T114248_noresp.xml')
    output_path = str(tmp_path)
    inventory_split(stxml_file, output_path)


# TODO:
#   * fix this test so that it actually uses test files with multiple networks (currently only AU)
#     to properly test splitting.
