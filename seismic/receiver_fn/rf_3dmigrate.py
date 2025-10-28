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
from seismic.receiver_fn.rf_ccp_util import ANTVolume
from seismic.misc import setup_logger

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('rf-h5-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output_h5_file', type=click.Path(exists=False, dir_okay=False), required=True)
@click.option('--relax-sanity-checks', is_flag=True, default=False, show_default=True,
              help='RF traces with amplitudes > 1.0 or troughs around onset time are dropped by default. '
                   'This option allows RF traces with amplitudes > 1.0 to pass through')
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
@click.option('--ant-model', type=click.Path(exists=True, dir_okay=False),
              default=None, show_default=True,
              help='ANT model in .txt format to extract Vs from.')
@click.option('--ant-model-max-depth', type=float,
              default=100, show_default=True,
              help='Maximum depth (km) up to which velocities from the ANT model are to be extracted. '
                   'Has no impact if --ant-model is not speficied.')
def main(rf_h5_file, output_h5_file, relax_sanity_checks, dz, max_depth, fmin, fmax,
         min_slope_ratio, ant_model, ant_model_max_depth):
    """Perform 3D migration of RFs
    RF_H5_FILE : Path to RFs in H5 format
    OUTPUT_H5_FILE: H5 output file name

    Example usage:
        mpirun -np 48 python rf_3dmigrate.py OA-ZRT-R-cleaned.h5 mig.h5 --min-slope-ratio 5 --fmin 0.1
    """
    log = setup_logger('__func__')

    am = None
    if(ant_model):
        log.info('Loading ANT model: {}'.format(ant_model))
        am = ANTVolume(ant_model, ant_model_max_depth)
    # end if

    m = Migrator(rf_filename=rf_h5_file, dz=dz, max_depth=max_depth,
                 min_slope_ratio=min_slope_ratio, ant_model=am,
                 logger=log)
    m.process_streams(output_h5_file, relax_sanity_checks=relax_sanity_checks, fmin=fmin, fmax=fmax)
# end

if __name__ == "__main__":
    # call main function
    main()  # pylint: disable=no-value-for-parameter
