import os
import json

import pytest
import numpy as np

from seismic.receiver_fn.moho_config import ConfigConstants as cc
from seismic.receiver_fn.pointsets2grid import _bounds, _grid
import seismic.receiver_fn.moho_workflow as mw
from tests.conftest import TESTS

@pytest.fixture()
def data_dir():
    return os.path.join(TESTS, 'test_seismic', 'receiver_fn', 'data')

@pytest.fixture()
def config(data_dir):
    cfg_file = os.path.join(data_dir, 'config.json')
    with open(cfg_file, 'r') as f:
        return json.load(f)

@pytest.mark.skip(reason="The moho workflow needs to be sorted out..")
def test_moho_workflow(tmpdir, data_dir, config):
    """
    Test full workflow and that outputs are equal to expected
    """
    ccp_1 = os.path.join(data_dir, 'ccp1.csv')
    ccp_2 = os.path.join(data_dir, 'ccp2.csv')

    test1_out = os.path.join(tmpdir, 'test1')
    config[cc.OUTPUT_DIR] = os.path.join(test1_out)
    # config[cc.OUTPUT_DIR] = '/home/bren/hs_new'
    for params in config[cc.METHODS]:
        if params[cc.NAME] == 'ccp1':
            params[cc.DATA] = ccp_1
        if params[cc.NAME] == 'ccp2':
            params[cc.DATA] = ccp_2

    test_conf1 = os.path.join(tmpdir, 'conf1.json') 
    with open(test_conf1, 'w') as f:
        json.dump(config, f)

    mw.run_workflow(test_conf1)
    
    expected_output = os.path.join(data_dir, 'expected_output')
    expected_grid = os.path.join(expected_output, 'grid.csv')
    with open(expected_grid, 'r') as f:
        expected_grid = np.loadtxt(f, delimiter=',', skiprows=2)
    expected_gradient = os.path.join(expected_output, 'gradient.csv')
    with open(expected_gradient, 'r') as f:
        expected_gradient = np.loadtxt(f, delimiter=',', skiprows=2)

    test1_grid = os.path.join(test1_out, cc.MOHO_GRID)

    with open(test1_grid, 'r') as f:
        nx = int(f.readline())
        ny = int(f.readline())
        moho_grid = np.loadtxt(f, delimiter=',')

    assert ny == 9
    assert nx == 17
    # np.testing.assert_equal(moho_grid, expected_grid)
    np.testing.assert_allclose(moho_grid, expected_grid, rtol=0.01)

    test1_grad = os.path.join(test1_out, cc.MOHO_GRAD)

    with open(test1_grad, 'r') as f:
        nx = int(f.readline())
        ny = int(f.readline())
        moho_gradient = np.loadtxt(f, delimiter=',')

    assert ny == 9
    assert nx == 17
    # np.testing.assert_equal(moho_grad, expected_grad)
    np.testing.assert_allclose(moho_gradient, expected_gradient, rtol=0.5)

    # Test other outputs exist
    assert os.path.exists(os.path.join(test1_out, cc.MOHO_PLOT + '.png'))
    assert os.path.exists(os.path.join(test1_out, cc.GMT_DIR, cc.MOHO_GRID_GMT))
    assert os.path.exists(os.path.join(test1_out, cc.GMT_DIR, cc.MOHO_GRAD_GMT))
    assert os.path.exists(os.path.join(test1_out, cc.GIS_DIR, cc.MOHO_GRID_GIS))
    assert os.path.exists(os.path.join(test1_out, cc.GIS_DIR, cc.MOHO_GRAD_GIS))
    for params in config[cc.METHODS]:
        assert os.path.exists(os.path.join(
            test1_out, cc.GMT_DIR, f'{params[cc.NAME]}{cc.LOCATIONS_GMT}'))
        assert os.path.exists(os.path.join(
            test1_out, cc.GIS_DIR, f'{params[cc.NAME]}{cc.LOCATIONS_GIS}.shp'))


def test_grid():
    n_x, x_grid, n_y, y_grid = _grid(np.array((130, 0)), np.array((131, 1)), 0.25)
    assert n_x == 5
    assert n_y == 5
    np.testing.assert_equal(x_grid[0], [130., 130.25, 130.5, 130.75, 131])
    np.testing.assert_equal(x_grid[1], [130., 130.25, 130.5, 130.75, 131])
    np.testing.assert_equal(x_grid[2], [130., 130.25, 130.5, 130.75, 131])
    np.testing.assert_equal(x_grid[3], [130., 130.25, 130.5, 130.75, 131])
    np.testing.assert_equal(x_grid[4], [130., 130.25, 130.5, 130.75, 131])
    np.testing.assert_equal(y_grid[0], [0., 0., 0., 0., 0.])
    np.testing.assert_equal(y_grid[1], [0.25, 0.25, 0.25, 0.25, 0.25])
    np.testing.assert_equal(y_grid[2], [0.5, 0.5, 0.5, 0.5, 0.5])
    np.testing.assert_equal(y_grid[3], [0.75, 0.75, 0.75, 0.75, 0.75])
    np.testing.assert_equal(y_grid[4], [1., 1., 1., 1., 1.])

    n_x, x_grid, n_y, y_grid = _grid(np.array((130, 0)), np.array((131, 1)), 0.5)
    assert n_x == 3
    assert n_y == 3
    np.testing.assert_equal(x_grid[0], [130., 130.5, 131.])
    np.testing.assert_equal(x_grid[1], [130., 130.5, 131.])
    np.testing.assert_equal(x_grid[2], [130., 130.5, 131.])
    np.testing.assert_equal(y_grid[0], [0., 0., 0.])
    np.testing.assert_equal(y_grid[1], [0.5, 0.5, 0.5])
    np.testing.assert_equal(y_grid[2], [1., 1., 1.])

    with pytest.raises(ValueError):
        _grid(np.array((130, 0)), np.array((129, 1)), 0.25)

    with pytest.raises(ValueError):
        _grid(np.array((130, 0)), np.array((131, -1)), 0.25)


def test_bounds():
    bb_min, bb_max = _bounds([(130, 0), (131, 1)], [(140, 10), (141, 11)], [10, 1, 20, 2])
    np.testing.assert_equal(bb_min, np.array([10, 1]))
    np.testing.assert_equal(bb_max, np.array([20, 2]))

    bb_min, bb_max = _bounds([(130, 0), (131, 1)], [(140, 10), (141, 11)], None)
    np.testing.assert_equal(bb_min, np.array([130, 0]))
    np.testing.assert_equal(bb_max, np.array([141, 11]))
