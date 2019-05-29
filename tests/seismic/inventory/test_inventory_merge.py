#!/usr/bin/env python

import os

from seismic.inventory.inventory_merge import inventory_merge


def test_inventory_merge_main(tmp_path):
    """High level test of main entry point for inventory merge.
    """
    self_path = os.path.dirname(os.path.abspath(__file__))
    test_data_path = os.path.join(self_path, 'data')
    iris_file = os.path.join(test_data_path, 'IRIS-ALL_tiny.xml')
    custom_file = os.path.join(test_data_path, 'INVENTORY_20190528T163156.xml')
    outfile_name = "INVENTORY_MERGED.xml"
    output_path = os.path.join(str(tmp_path), outfile_name)
    assert not os.path.exists(output_path)
    inventory_merge(iris_file, custom_file, output_path, test_mode=True)
    assert not os.path.exists(os.path.splitext(custom_file)[0] + ".pkl")
    assert os.path.exists(output_path)

# TODO: Add additional test functions that:
#   * test individual functions from inventory_merge module
#   * test the correctness of the content in the output_path file


if __name__ == "__main__":
    import tempfile
    test_inventory_merge_main(tempfile.mkdtemp())
