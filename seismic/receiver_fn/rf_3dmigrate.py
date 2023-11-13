#!/bin/env python
"""
Description:
    Implements migration algorithm as described in Frassetto et al. (2010):
      Improved imaging with phase-weighted common conversion point stacks
      of receiver functions (GJI)

References:

CreationDate:   3/15/18
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     3/15/18   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import sys
import logging
import click
from seismic.receiver_fn.rf_ccp_util import Migrator

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
log = logging.getLogger('rf_3dmigrate')

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('rf-h5-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output_h5_file', type=click.Path(exists=False, dir_okay=False), required=True)
@click.option('--dz', type=float, default=0.1, show_default=True,
              help='Depth-step (km)')
@click.option('--max-depth', type=click.FloatRange(0, 750), default=150, show_default=True,
              help='Maximum depth (km) of profile')
@click.option('--fmin', type=float, default=None, show_default=True,
              help="Lowest frequency for bandpass filter; default is None."
                   "If only --fmin in provided, a highpass filter is aplied.")
@click.option('--fmax', type=float, default=None, show_default=True,
              help="Highest frequency for bandpass filter; default is None."
                   "If only --fmax is provided, a lowpass filter is applied.")
@click.option('--min-slope-ratio', type=float, default=-1, show_default=True,
              help='Apply filtering to the RFs based on the "slope_ratio" metric '
                   'that indicates robustness of P-arrival. Typically, a minimum '
                   'slope-ratio of 5 is able to pick out strong arrivals. The '
                   'default value of -1 does not apply this filter')
def main(rf_h5_file, output_h5_file, dz, max_depth, fmin, fmax, min_slope_ratio):
    """Perform 3D migration of RFs
    RF_H5_FILE : Path to RFs in H5 format
    OUTPUT_H5_FILE: H5 output file name

    Example usage:
        mpirun -np 48 python rf_3dmigrate.py OA-ZRT-R-cleaned.h5 mig.h5 --min-slope-ratio 5 --fmin 0.1
    """
    log.setLevel(logging.DEBUG)

    m = Migrator(rf_filename=rf_h5_file, dz=dz, max_depth=max_depth,
                 min_slope_ratio=min_slope_ratio, logger=log)
    m.process_streams(output_h5_file, fmin=fmin, fmax=fmax)
# end

if __name__ == "__main__":
    # call main function
    main()  # pylint: disable=no-value-for-parameter
