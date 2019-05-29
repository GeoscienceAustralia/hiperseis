#!/usr/bin/env python

import pytest

from tests.mocks.inventory.mock_fdsn_xml import MockIrisResponse

@pytest.fixture(scope='session')
def iris_mocker():
    m = MockIrisResponse()
    m.start()
    return m
