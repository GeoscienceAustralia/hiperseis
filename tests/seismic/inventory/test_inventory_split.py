#!/usr/bin/env python

import os
import glob

from obspy import read_inventory

from seismic.inventory.inventory_split import inventory_split


def test_inventory_split_main(tmp_path):
    """High level test of main entry point for inventory split.
    """
    self_path = os.path.dirname(os.path.abspath(__file__))
    test_data_path = os.path.join(self_path, 'data')
    stxml_file = os.path.join(test_data_path, 'INVENTORY_20190528T163156_MERGED.xml')
    output_path = str(tmp_path)
    inventory_split(stxml_file, output_path)
    assert not os.path.exists(os.path.splitext(stxml_file)[0] + ".pkl")
    found_files = glob.glob(os.path.join(output_path, 'network_*.xml'))
    assert len(found_files) > 1
    # Confirm that each file can be read by obspy
    for f in found_files:
        inv = read_inventory(f)
        assert len(inv.networks) == 1


if __name__ == "__main__":
    import tempfile
    test_inventory_split_main(tempfile.mkdtemp())
