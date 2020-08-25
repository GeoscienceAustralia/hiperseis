import os
import json

import pytest
import numpy as np

from seismic.receiver_fn.moho_config import ConfigConstants as cc
import seismic.receiver_fn.moho_workflow as mw
from tests.conftest import TESTS

@pytest.fixture()
def data_dir():
    return data_dir = os.path.join(TESTS, 'test_seismic', 'receiver_fn', 'data')

@pytest.fixture()
def config(data_dir):
    cfg_file = os.path.join(data_dir, 'config.json')
    with open(cfg_file, 'r') as f:
        return json.load(f)

def test_moho_workflow(tmpdir, data_dir, config):
    ccp_1 = os.path.join(data_dir, 'ccp1.csv')
    ccp_2 = os.path.join(data_dir, 'ccp2.csv')

    test1_out = os.path.join(tmpdir, 'test1')
    config[cc.OUTPUT_DIR] = os.path.join(test1_out)
    for params in config[cc.METHODS]:
        if params[cc.NAME] == 'ccp_1':
            params[cc.DATA] = ccp_1
        if params[cc.NAME] == 'ccp_2':
            params[cc.DATA] = ccp_2

    test_conf1 = os.path.join(tmpdir, 'conf1.json') 
    with open(test_conf1, 'w') as f:
        json.dump(config, f)

    mw.run_workflow(test_conf1)
    
    expected_output = os.path.join(data_dir, 'expected_output')
    expected_grid = os.path.join(expected_output, 'moho_grid.csv')
    with open(expected_grid, 'r') as f:
        expected_grid = np.loadtxt(f, delimiter=',', skiprows=2)
    expected_grad = os.path.join(expected_output, 'moho_gradient.csv')
    with open(expected_grad, 'r') as f:
        expected_grad = np.loadtxt(f, delimiter=',', skiprows=2)

    test1_grid = os.path.join(test1_out, cc.MOHO_GRID)
    test1_grad = os.path.join(test1_out, cc.MOHO_GRAD)

    with open(test1_grid, 'r') as f:
        nx = int(f.readline())
        ny = int(f.readline())
        moho_grid = np.loadtxt(f, delimiter=',')

    assert ny == 9
    assert nx == 17
    np.testing.assert_equal(moho_grid, expected_grid)

    with open(test1_grad, 'r') as f:
        nx = int(f.readline())
        ny = int(f.readline())
        moho_grad = np.loadtxt(f, delimiter=',')

    assert ny == 9
    assert nx == 17
    np.testing.assert_equal(moho_grad, expected_grad)
    

