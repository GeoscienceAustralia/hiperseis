import os

import click
import numpy as np

from seismic.receiver_fn.moho_config import ConfigConstants as cc

# Map spellings of network codes in hand-digitised data to 
# their FDS equivalent
NETWORK_CODE_MAPPINGS = {'AQT': 'AQ',
                         'OA2': 'OA'}

# Characters on station names to remove
SPECIAL_CHARS = ['*', '.']


def add_weights(data_file, sample_weights=None):
    col_names = ['Sta', 'Lon', 'Lat', 'Depth']
    data = np.genfromtxt(data_file, delimiter=',', dtype=None, encoding=None,
                         names=col_names, usecols=(0, 1, 2, 3))

    weights = np.ones_like(data['Depth'])
    if sample_weights is not None:
        if isinstance(sample_weights, (str, os.PathLike)) and os.path.exists(sample_weights):
            with open(sample_weights, 'r') as fr:
                lines = fr.readlines()
            sta_weight = zip(*[l.split() for l in lines])
        elif isinstance(sample_weights, dict):
            sta_weight = sample_weights.items()

        for s, w in sta_weight:
            weights[np.where(data['Sta'] == s)] = w
    else:
        print("'sample_weights' not recognised. Must be a file of format 'STATION WEIGHT' or a "
              "dict of '{STATION: WEIGHT}'. All weights will be set to 1.0")
            
    data = np.array((data['Sta'], data['Lon'], data['Lat'], data['Depth'], weights)).T
    np.savetxt(data_file, data, fmt=['%s', '%s', '%s', '%s', '%s'], delimiter=',',
               header=','.join(col_names + ['Weight']))

@click.command()
@click.argument('data_file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.option('--weights-file', type=click.Path(exists=True, dir_okay=False), required=False)
def main(data_file, weights_file=None):
    add_weights(data_file, weights_file)


if __name__ == '__main__':
    main()
