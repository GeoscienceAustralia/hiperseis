import os
import pytest
from iloc_rstt.iloc_event import ILocCatalog, DELTA


@pytest.fixture
def analyst_event(test_dir):
    return os.path.join(test_dir, 'mocks', 'events',
                        'analyst_event_samples', '772009.xml')


def test_stations_in_range(analyst_event):
    cat = ILocCatalog(analyst_event)
    assert len(cat.events) == 1  # there is only one event in the xml
    cat.update()
    assert len(cat.iloc_events) == 1  # should create only one iloc event

    for e in cat.iloc_events:
        assert all(e.new_stations[DELTA] < 120)

