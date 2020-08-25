import click
import numpy as np

from seismic.receiver_fn.moho_config import ConfigConstants as cc

# Map spellings of network codes in hand-digitised data to 
# their FDS equivalent
NETWORK_CODE_MAPPINGS = {'AQT': 'AQ',
                         'OA2': 'OA'}

# Characters on station names to remove
SPECIAL_CHARS = ['*', '.']


def add_weights(data_file, weights=None):
    col_names = ['Sta', 'Lon', 'Lat', 'Depth']
    data = np.genfromtxt(data_file, delimiter=',', dtype=None, encoding=None,
                         names=col_names, usecols=(0, 1, 2, 3))
    if weights is not None:
        pass # TODO: weights from file or dict
    else:
        weights = np.ones_like(data['Depth'])
    data = np.array((data['Sta'], data['Lon'], data['Lat'], data['Depth'], weights)).T
    np.savetxt(data_file, data, fmt=['%s', '%s', '%s', '%s', '%s'], delimiter=',',
               header=','.join(col_names + ['Weight']))

@click.command()
@click.argument('data_file', type=click.Path(exists=True, dir_okay=False), required=True)
def main(data_file):
    add_weights(data_file)


if __name__ == '__main__':
    main()
