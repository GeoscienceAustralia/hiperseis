import pytest
import datetime
import yaml
from seismic import config

conf = {
    'inputs':
    [{'files': [{'file': 'tests/mocks/data/ev0_6.a01.gse2'},
                {'file': 'tests/mocks/data/ev0_6.a02.gse2'}],
      'name': 'my miniseed files', 'type': 'miniseed'},
     {'name': 'my events',
      'events': [{'event_id': 'ga2017abcdefg'},
                 {'event_id': 'ga2017hijklmn'}],
      'type': 'events'},
     {'times': {'end_time': datetime.datetime(2017, 3, 28, 17, 18, 30),
                'start_time': datetime.datetime(2017, 3, 28, 16, 18, 30)},
      'name': 'my time range', 'type': 'time'}]
    }


def test_conflict(random_filename):
    yaml_file = random_filename(ext='.yaml')
    with open(yaml_file, 'w') as outfile:
        yaml.dump(conf, outfile, default_flow_style=False)
    with pytest.raises(config.ConfigException, message='Too many inputs '
                                                       'provided'):
        config.Config(yaml_file)
