#!/usr/bin/env python
"""
Description:
    Plots CCP vertical and depth slices

References:

CreationDate:   03/09/21
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     03/09/21   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import os, sys
import logging
import numpy as np
import click
import json
import matplotlib.pyplot as plt
from PIL.PngImagePlugin import PngImageFile, PngInfo

from seismic.receiver_fn.rf_ccp_util import CCPVolume, CCP_VerticalProfile, CCP_DepthProfile, Gravity
from utils.sgrid import SGrid

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.group(context_settings=CONTEXT_SETTINGS)
def groups():
  pass
# end func

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
log = logging.getLogger('plot_ccp')

def read_profile_def(fname, ccpVolume, type):
    if(type != 'vertical' and type != 'depth'):
        assert 0, 'Profile-type {} not supported. Aborting..'.format(type)
    # end if

    profiles = []

    def validate_start_end(line):
        start, end = None, None
        try:
            start, end = [item.strip() for item in line.split('-')]
            for item in [start, end]:
                if (item not in ccpVolume._meta.keys()):
                    raise ValueError('{} not found in {}. '
                                     'Available stations are: {}. '
                                     'Aborting..'.format(item, ccpVolume._fn_list,
                                                         list(ccpVolume._meta.keys())))
                # end if
            # end for
        except Exception as e:
            print(str(e))
            assert 0, 'Invalid entry in profile-definition file. The format should be ' \
                      'NET.STA.LOC-NET.STA.LOC. Aborting..'
        # end try

        return start, end
    # end func

    fh = open(fname, 'r')
    for iline, line in enumerate(fh.readlines()):
        if(not line.strip()): continue

        if(type=='vertical'):
            start, end = validate_start_end(line)
            profiles.append([start, end])
        else:
            if(iline == 0):
                start, end = validate_start_end(line)
                profiles.append([start, end])
            else:
                try:
                    depth = np.float(line)
                    profiles.append(depth)
                except Exception as e:
                    print(str(e))
                    assert 0, 'Expected a single float value for depth. Aborting..'
                # end try
            # end if
        # end if
    # end for

    if(type=='depth' and len(profiles)<2):
        assert 0, 'Invalid profile-definition file {}. First line should be NET.STA.LOC-NET.STA.LOC ' \
                  'followed by depth(s) (km) in each following line. Aborting..'.format(fname)
    # end if

    return profiles
# end func

@click.command(name='vertical', context_settings=CONTEXT_SETTINGS)
@click.argument('ccp-h5-volume', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('profile-def', type=click.Path(exists=True, dir_okay=False), required=True)
@click.option('--gravity-grid', type=click.Path(exists=True, dir_okay=False), default=None, show_default=True,
              help='Gravity grid in tif/ers format. If provided, a gravity line-plot is added to each'
                   ' vertical profile. Note that this parameter is ignored if the output-format is set '
                   'to txt')
@click.option('--mt-sgrid', type=click.Path(exists=True, dir_okay=False), default=None, show_default=True,
              help='MT grid in sgrid format. If provided, an additional column with resistivity values are '
                   'added in the output text file. Note that this parameter is ignored if the output-format '
                   'is set to png')
@click.option('--mt-utm-zone', type=str, default=None, show_default=True,
              help='UTM zone for mt-grid e.g. 53S. This parameter is ignored if mt-grid is not provided')
@click.option('--dx', type=float, default=5, show_default=True,
              help='Horizontal distance-step (km) between start and end locations of a profile')
@click.option('--dz', type=float, default=0.5, show_default=True,
              help='Depth-step (km)')
@click.option('--max-depth', type=click.FloatRange(0, 750), default=100, show_default=True,
              help='Maximum depth (km) of profile')
@click.option('--swath-width', type=float, default=40, show_default=True,
              help='CCP amplitudes are averaged over grid elements across a swath defined by '
                   '±swath-width/2 (km) on either side of a profile. A larger swath-width results in '
                   'a smoother image -- ideally, swath-width should be guided by areal distribution '
                   'of stations around a profile')
@click.option('--ds', type=float, default=10, show_default=True,
              help='Horizontal distance-step (km) along swath, from -swath-width/2 to +swath-width/2 ')
@click.option('--extend', type=float, default=50, show_default=True,
              help='The length (km) by which to extend either end of a profile, since CCP coverage '
                   'extends laterally with depth')
@click.option('--cell-radius', type=float, default=20, show_default=True,
              help='CCP amplitudes at each grid element are computed through a combination of '
                   'inverse-distance- and instantaneous-phase-weighting applied to CCP values '
                   'that fall within a disk, defined by cell-radius (km) and dz. Cell-radius '
                   'should be guided by dx, ds, and areal distribution of stations in the CCP '
                   'volume')
@click.option('--idw-exponent', type=float, default=2, show_default=True,
              help='Exponent used in inverse-distance-weighted interpolation, which determines the '
                   'relative contribution of near and far values. A larger exponent diminishes the '
                   'contribution of faraway values')
@click.option('--pw-exponent', type=float, default=1, show_default=True,
              help='Exponent used in instantaneous phase-weighting of CCP amplitudes')
@click.option('--amplitude-min', type=float, default=-0.2, show_default=True,
              help='Minimum amplitude for colorbar normalization')
@click.option('--amplitude-max', type=float, default=0.2, show_default=True,
              help='Maximum amplitude for colorbar normalization')
@click.option('--max-station-dist', type=float, default=10, show_default=True,
              help='Stations closer than max-station-dist (km) to a profile are labelled on plots')
@click.option('--output-folder', type=str, default='./', show_default=True,
              help='Output folder')
@click.option('--output-format', type=click.Choice(['png', 'txt']), default='png', show_default=True,
              help='Output format')
@click.option('--plot-aspect', type=click.Choice(['auto', 'equal']), default='auto', show_default=True,
              help='Aspect ratio of CCP plot. Default is "auto", which applies a vertical exaggeration.')
def vertical(ccp_h5_volume, profile_def, gravity_grid, mt_sgrid, mt_utm_zone, dx, dz, max_depth, swath_width, ds,
             extend, cell_radius, idw_exponent, pw_exponent, amplitude_min, amplitude_max, max_station_dist,
             output_folder, output_format, plot_aspect):
    """ Plot CCP vertical profile \n
    CCP_H5_VOLUME: Path to CCP volume in H5 format (output of rf_3dmigrate.py) or a text file containing paths to CCP volumes in H5 format. The latter is useful for generating profiles spanning multiple neighboring networks \n
    PROFILE_DEF: text file containing start and end locations of each vertical profile as\n
                 NET.STA.LOC-NET.STA.LOC in each line\n\n

    Example usage:\n
        python plot_ccp.py vertical OA_ccp_volume.h5 slice_def.txt\n
    """
    log.setLevel(logging.DEBUG)

    # create output folder if necessary
    if not os.path.exists(output_folder): os.makedirs(output_folder)

    # sanity checks
    if (mt_sgrid and not mt_utm_zone): assert 0, 'UTM zone for {} not specified. ' \
                                                'See help for --mt-utm-zone'.format(mt_sgrid)

    log.info('Loading CCP volume {}..'.format(ccp_h5_volume))
    vol = CCPVolume(ccp_h5_volume)

    log.info('Loading profile definition file {}..'.format(profile_def))
    profiles = read_profile_def(profile_def, vol, 'vertical')

    gravity = None
    if(gravity_grid):
        log.info('Loading gravity grid {}..'.format(gravity_grid))
        gravity = Gravity(gravity_grid)
    # end if

    mt = None
    if(mt_sgrid):
        log.info('Loading MT SGrid {}..'.format(mt_sgrid))
        mt = SGrid(mt_sgrid, mt_utm_zone)
    # end if

    for s, e in profiles:
        log.info('Processing {}..'.format('-'.join([s, e])))

        vprof = CCP_VerticalProfile(vol, s, e, dx=dx, dz=dz, max_depth=max_depth,
                                    swath_width=swath_width, ds=ds, extend=extend,
                                    cell_radius=cell_radius, idw_exponent=idw_exponent,
                                    pw_exponent=pw_exponent, max_station_dist=max_station_dist)

        fname = os.path.join(output_folder, '{}-{}.{}'.format(s, e, output_format))

        if(output_format.lower() == 'txt'):
            ggx, ggd = np.meshgrid(vprof._gx, vprof._gd)

            grid = np.zeros((vprof._grid.shape[0], 5))
            grid[:, 0] = vprof._grid[:, 1]
            grid[:, 1] = vprof._grid[:, 2]
            grid[:, 2] = ggx.flatten()
            grid[:, 3] = ggd.flatten()
            grid[:, 4] = vprof._grid_vals.flatten()

            if(mt):
                grid = np.hstack((grid, np.atleast_2d(np.zeros(grid.shape[0])).T))

                grid[:, 5] = mt.get_values(grid[:, 0], grid[:, 1], grid[:, 3]*1e3, p=3)

                """
                fig, ax = plt.subplots(1)
                fig.set_size_inches(25,10)
                fig.set_dpi(300)

                vals = np.reshape(np.array(grid[:, 5]), (vprof._gd.shape[0], vprof._gx.shape[0]))
                ax.pcolor(np.log1p(vals), cmap='seismic_r')
                plt.savefig(fname+'.png')
                """

                np.savetxt(fname, grid, header='lon lat distance(km) depth(km) amplitude(arb. units) resistivity(ohm m)',
                           fmt='%10.10f')
            else:
                np.savetxt(fname, grid, header='lon lat distance(km) depth(km) amplitude(arb. units)', fmt='%10.10f')
            # end if
        else:
            fig, gax, ax = None, None, None
            if(gravity):
                fig, (gax, ax) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1, 5]}, sharex=True)
            else:
                fig, ax = plt.subplots()
            # end if

            fig.set_size_inches((20, 10))
            fig.set_dpi(300)

            ax.set_aspect(plot_aspect)

            # plot profile
            vprof.plot(ax, amp_min=amplitude_min, amp_max=amplitude_max, gax=gax, gravity=gravity)

            fig.savefig(fname)

            if(output_format.lower() == 'png'):
                # write meta-data for translating digitization
                d = {'lon1': vprof._lon1, 'lat1': vprof._lat1,
                     'lon2': vprof._lon2, 'lat2': vprof._lat2,
                     'x1':np.min(vprof._gx), 'y1':np.min(vprof._gd),
                     'x2':np.max(vprof._gx), 'y2':np.max(vprof._gd)}
                sd = json.dumps(d)

                img = PngImageFile(fname)
                meta = PngInfo()
                meta.add_text('profile_meta', sd)

                img.save(fname, pnginfo=meta)
                img.close()
            # end if
        # end if
    # end for
# end

@click.command(name='depth', context_settings=CONTEXT_SETTINGS)
@click.argument('ccp-h5-volume', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('profile-def', type=click.Path(exists=True, dir_okay=False), required=True)
@click.option('--mt-sgrid', type=click.Path(exists=True, dir_okay=False), default=None, show_default=True,
              help='MT grid in sgrid format. If provided, an additional column with resistivity values are '
                   'added in the output text file. Note that this parameter is ignored if the output-format '
                   'is set to png')
@click.option('--mt-utm-zone', type=str, default=None, show_default=True,
              help='UTM zone for mt-grid e.g. 53S. This parameter is ignored if mt-grid is not provided')
@click.option('--dx', type=float, default=5, show_default=True,
              help='Longitudinal distance-step (km) between start and end locations')
@click.option('--dy', type=float, default=5, show_default=True,
              help='Latitudinal distance-step (km) between start and end locations')
@click.option('--dz', type=float, default=1, show_default=True,
              help='Thickness (km) of layer at "depth" over which CCP values are averaged')
@click.option('--extend', type=float, default=50, show_default=True,
              help='The length (km) by which to extend either end of a profile, since CCP coverage '
                   'extends laterally with depth')
@click.option('--cell-radius', type=float, default=40, show_default=True,
              help='CCP amplitudes at each grid element are computed through a combination of '
                   'inverse-distance- and instantaneous-phase-weighting applied to CCP values '
                   'that fall within a disk, defined by cell-radius (km) and dz. Cell-radius '
                   'should be guided by dx and areal distribution of stations in the CCP '
                   'volume')
@click.option('--idw-exponent', type=float, default=2, show_default=True,
              help='Exponent used in inverse-distance-weighted interpolation, which determines the '
                   'relative contribution of near and far values. A larger exponent diminishes the '
                   'contribution of faraway values')
@click.option('--pw-exponent', type=float, default=1, show_default=True,
              help='Exponent used in instantaneous phase-weighting of CCP amplitudes')
@click.option('--amplitude-min', type=float, default=-0.2, show_default=True,
              help='Minimum amplitude for colorbar normalization')
@click.option('--amplitude-max', type=float, default=0.2, show_default=True,
              help='Maximum amplitude for colorbar normalization')
@click.option('--output-folder', type=str, default='./', show_default=True,
              help='Output folder')
@click.option('--output-format', type=click.Choice(['png', 'txt']), default='png', show_default=True,
              help='Output format')
def depth(ccp_h5_volume, profile_def, mt_sgrid, mt_utm_zone, dx, dy, dz, extend, cell_radius,
          idw_exponent, pw_exponent, amplitude_min, amplitude_max, output_folder, output_format):
    """ Plot CCP depth profile \n
    CCP_H5_VOLUME: Path to CCP volume in H5 format (output of rf_3dmigrate.py) or a text file containing paths to CCP volumes in H5 format. The latter is useful for generating profiles spanning multiple neighboring networks \n
    PROFILE_DEF: text file containing start and end locations, defining the rectangular region spanned by each
                 depth profile, in the first line as NET.STA.LOC-NET.STA.LOC, followed by depth(s) (km) in
                 subsequent lines \n\n

    Example usage:\n
        python plot_ccp.py depth OA_ccp_volume.h5 slice_def.txt\n
    """
    log.setLevel(logging.DEBUG)

    # create output folder if necessary
    if not os.path.exists(output_folder): os.makedirs(output_folder)

    # sanity checks
    if (mt_sgrid and not mt_utm_zone): assert 0, 'UTM zone for {} not specified. ' \
                                                'See help for --mt-utm-zone'.format(mt_sgrid)

    log.info('Loading CCP volume {}..'.format(ccp_h5_volume))
    vol = CCPVolume(ccp_h5_volume)

    log.info('Loading profile definition file {}..'.format(profile_def))
    profiles = read_profile_def(profile_def, vol, 'depth')

    mt = None
    if(mt_sgrid):
        log.info('Loading MT SGrid {}..'.format(mt_sgrid))
        mt = SGrid(mt_sgrid, mt_utm_zone)
    # end if

    s, e = profiles[0]
    for depth in profiles[1:]:
        log.info('Processing depth {}..'.format(depth))

        dprof = CCP_DepthProfile(vol, s, e, depth, dx=dx, dy=dy, dz=dz, extend=extend,
                                 cell_radius=cell_radius, idw_exponent=idw_exponent,
                                 pw_exponent=pw_exponent)

        fname = os.path.join(output_folder, '{}-{}-depth-{}.{}'.format(s, e, int(depth), output_format))

        if(output_format.lower() == 'txt'):
            ggx, ggy = np.meshgrid(dprof._gx, dprof._gy)

            grid = np.zeros((dprof._grid.shape[0], 5))
            grid[:, 0] = dprof._grid[:, 1]
            grid[:, 1] = dprof._grid[:, 2]
            grid[:, 2] = ggx.flatten()
            grid[:, 3] = ggy.flatten()
            grid[:, 4] = dprof._grid_vals.flatten()

            if(mt):
                grid = np.hstack((grid, np.atleast_2d(np.zeros(grid.shape[0])).T))

                grid[:, 5] = mt.get_values(grid[:, 0], grid[:, 1],
                                           np.ones(len(dprof._grid))*depth*1e3, p=3)

                """
                fig, ax = plt.subplots(1)
                fig.set_size_inches(25,10)
                fig.set_dpi(300)

                vals = np.reshape(np.array(grid[:, 5]), (dprof._gy.shape[0], dprof._gx.shape[0]))
                ax.pcolor(np.log1p(vals), cmap='seismic_r')
                plt.savefig(fname+'.png')
                """

                np.savetxt(fname, grid, header='lon lat x(km) y(km) amplitude(arb. units) resistivity(ohm m)',
                           fmt='%10.10f')
            else:
                np.savetxt(fname, grid, header='lon lat x(km) y(km) amplitude(arb. units)', fmt='%10.10f')
            # end if
        else:
            fig, ax = plt.subplots()
            fig.set_size_inches((10, 10))
            fig.set_dpi(300)

            # plot profile
            dprof.plot(ax, amp_min=amplitude_min, amp_max=amplitude_max)

            fig.savefig(fname)
        # end if
    # end for
# end

groups.add_command(vertical)
groups.add_command(depth)

if __name__ == "__main__":
    groups()
# end if
