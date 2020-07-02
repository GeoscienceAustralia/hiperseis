#!/usr/bin/env python
"""Plot vespagrams from receiver function data
"""

import os
import click
from past.builtins import xrange

import numpy as np
import matplotlib
# comment out the line below if you want to visualise the figure
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.tri as tri

# from obspy.core.utcdatetime import UTCDateTime as UTC

# Here are the libraries to deal with RFSTREAM, it uses obspy classes for event and station
import rf

# from seismic.units_utils import KM_PER_DEG

# pylint: disable=invalid-name


# Definition of the simple 1D Earth model, remember each interface will give one Ps conversion
# you can add shallow layer to study Ps conversions in sedimentary basins
z = np.array([0, 35, 35, 165])
vp = np.array([5.8, 6.5, 8.04, 8.17])
vs = np.array([3.5, 3.9, 4.48, 4.51])
simple_model = rf.simple_model.SimpleModel(z, vp, vs)

# For sake of simplicity we will generate separate model for sedimentary basin
zb = np.array([0, 5, 5, 300])
vp = np.array([4.5, 4.5, 6.5, 8.17])
vs = np.array([2.5, 2.5, 3.5, 4.51])
basin_model = rf.simple_model.SimpleModel(zb, vp, vs)


# its not advised to use standard models because they create large number of conversions at each layer
#simple_model=rf.simple_model.load_model(fname='iasp91')


#-------------Main---------------------------------

@click.command()
@click.argument('input-h5-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output-pdf-file', type=click.Path(exists=False), required=True)
def main(input_h5_file, output_pdf_file):
    ''' This program composes vespagrams to identify RF converted phases and their multiples
        please refer to Tian et al. GRL 2005 VOL. 32, L08301, doi:10.1029/2004GL021885 for good examples

        input  - H5 file with receiver functions
        output - PDF files to print

        Dependencies - rf and obspy packages beside other standard python packages
        The image is composed using triangulation. It gives good results but block median or mean must
        be implemented at some stage to reduce size of PDF.
    '''

    stream = rf.read_rf(input_h5_file, 'H5')
    rf_type = 'LQT-Q '
    filter_type = 'bandpass'
    freqmin = 0.03
    freqmax = 0.5

    # we use a zero-phase-shift band-pass filter using 2 corners. This is done in two runs forward and backward,
    # so we end up with 4 corners de facto.
    # Lets assume we have LQT orientation
    selected_stream = stream.select(component='Q').filter(filter_type, freqmin=freqmin, freqmax=freqmax,
                                                          corners=2, zerophase=True).interpolate(10)


    # if none lets try ZRT
    if not selected_stream:
        selected_stream = stream.select(component='R').filter(filter_type, freqmin=freqmin, freqmax=freqmax,
                                                              corners=2, zerophase=True).interpolate(10)
        rf_type = 'ZRT-R '
    # end if

    if not selected_stream:
        print("Tried Q and R components but neither found, quitting...")
        exit(-1)
    # end if

    station_list = []

    # here we collect station names but maybe ID is more appropriate in case of having the same station names
    # in different deployments

    for tr in selected_stream:
        station_list.append(tr.stats.station)
        net = tr.stats.network
    # end for

    pdf = PdfPages(output_pdf_file)
    case_description = rf_type + filter_type + ' ' + str(freqmin) + '-' + str(freqmax) + ' Hz'
    pdf.attach_note(case_description)
    d = pdf.infodict()
    d['Title'] = rf_type + 'RF vespagrams of ' + net + ' network'
    d['Keywords'] = case_description

    station_list = np.unique(np.array(station_list))
    print("Gathered ", len(station_list), " stations")

    # Define layout of the page outer_grid
    columns = 3
    rows = 2
    frame = 0
    figure = 1

    # ------------------------------------------
    # Main loop here over all stations
    for i, station in enumerate(station_list):
        if frame == 0:
            printed = False
            fig = plt.figure(figsize=(11.69, 8.27))
            outer_grid = gridspec.GridSpec(columns, rows, wspace=0.2, hspace=0.2)
        # end if

        print("Station ", station, i + 1, " of ", station_list.shape[0])
        traces = selected_stream.select(station=station)
        print('Contains: ', len(traces), ' events')

        # we choose short RF to simplify and speed up the processing
        # from -5 to 20 seconds and slowness range from 5 to 9 s/deg
        # its enough to see multiples and possible LAB conversion at ~19 sec (~160km)

        traces = traces.trim2(-5, 20, 'onset')
        moved = []
        slow = []

        for tr in traces:
            tr.normalize()

            # This 'if' block is designed to check correct data placement on vespagram to
            # trace the logic (debugging purposes)
            DEBUG_PLACEMENT = False
            if DEBUG_PLACEMENT and (tr.stats.slowness > 6.) and (tr.stats.slowness < 7.):
                print('altered')
                data = tr.data.copy()
                print(data.shape, tr.stats.delta)
                # 500 below corresponds to 0 with sampling rate of 100Hz
                # TODO: !Change these array indices to be computed, not magic numbers!
                data[500:800] = 1.
                moved.append(data)
            else:
                moved.append(tr.data.copy()/np.max(np.abs(tr.data)))
                slow.append(tr.stats.slowness)
            # end if
        # end for

        print("Slowness min and max: ", np.min(slow), np.max(slow))
        slow.append(np.min(slow) - 0.1)
        moved.append(np.zeros(traces[0].data.shape))
        slow.append(np.max(slow) + 0.1)
        moved.append(np.zeros(traces[0].data.shape))

        slow = np.array(slow)
        idx = np.argsort(slow)
        moved = np.nan_to_num(np.array(moved))
        # moved = np.array(moved)
        slow = slow[idx]
        moved = moved[idx, :]
        z = moved.copy()

        # Some simple stacking to reduce data size on the image, this block can be safely commented out
        idx = []
        idx.append(True) # first row with zeroes
        slo_cum = 0.
        elements = 1
        for j in xrange(1, slow.shape[0] - 2):
            if np.abs(slow[j+1] - slow[j]) < 0.1 and slo_cum < 0.2:
                slow[j+1] = (slow[j] + slow[j+1])/2.
                moved[j, :] = moved[j, :]*elements
                moved[j+1, :] = np.sum(moved[j:j+2, :], axis=0)/(elements + 1)
                elements = elements + 1
                idx.append(False)
                slo_cum = slo_cum + np.abs(slow[j+1] - slow[j])
            else:
                idx.append(True)
                slo_cum = 0
                elements = 1
            # end if
        # end for

        idx.append(True) # before last
        idx.append(True) # last row with zeroes
        idx = np.array(idx)
        print(idx.shape, slow.shape, moved.shape)

        slow = slow[idx]
        moved = moved[idx, :]
        z = moved.copy()
        # ------------------------------ end of stacking ------------------------

        # print('minmax',np.min(z),np.max(z))
        x = np.array(list(range(moved.shape[1])))*traces[0].stats.delta - 5.
        x = np.repeat([x], moved.shape[0], axis=0)
        y = np.ones((moved.shape[0], moved.shape[1]))

        phase_Ps = []
        phase_Pms = []
        phase_PpPmS = []
        phase_PpSmS = []
        phase_slow = []
         # basin part
        phase_Pbs = []
        phase_PpPbs = []

        for j in xrange(slow.shape[0]):
            y[j, :] = y[j, :]*slow[j]
            phase_Ps.append(simple_model.calculate_delay_times(slow[j], phase='PS'))
            phase_Pms.append(simple_model.calculate_delay_times(slow[j], phase='PmS'))
            phase_PpPmS.append(simple_model.calculate_delay_times(slow[j], phase='PpPmS'))
            phase_PpSmS.append(simple_model.calculate_delay_times(slow[j], phase='PpSmS'))
            phase_slow.append(np.ones(phase_Ps[-1].shape[0])*slow[j])

            # basin, we will use reflection at the top layer only
            if zb.size > 0:
                phase_Pbs.append(basin_model.calculate_delay_times(slow[j], phase='PS'))
                phase_PpPbs.append(basin_model.calculate_delay_times(slow[j], phase='PpPmS'))
            # end if
        # end for

        xi = np.linspace(-5, 20, 200)
        yi = np.linspace(0, 9, 400)

        # Gridding the data using triangulation. standard gridding doesn't work well here
        triang = tri.Triangulation(x.flatten(), y.flatten())
        interpolator = tri.LinearTriInterpolator(triang, z.flatten())
        xi, yi = np.meshgrid(xi, yi)
        zi = interpolator(xi, yi)

        # Define two plots as inner_grid to place them inside one cell of outer_grid
        inner_grid = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer_grid[frame], wspace=0.5, hspace=0.)
        ax1 = plt.Subplot(fig, inner_grid[0])
        ax2 = plt.Subplot(fig, inner_grid[1], sharex=ax1)
        lim = np.max(np.abs(zi[zi < 0])*0.5)
        levels = np.linspace(-lim, lim, 15)
        # print("Levels ",-lim,lim)
        cmap = plt.cm.jet
        cs = ax1.contourf(xi, yi, zi, levels=levels, extend='both', cmap=cmap)
        cs.cmap.set_under('k')
        cs.cmap.set_over('k')
        ax1.set_ylim(5, 9)
        ax1.set_xlim(-5, 20)
        ax1.plot(phase_Ps, slow, color='black')      # direct conversion, positive amplitude
        ax1.plot(phase_PpPmS, slow, color='crimson') # multiples,         positive amplitude
        ax1.plot(phase_PpSmS, slow, color='purple')  # multiples,         negative amplitude

        ax1.annotate('Pms', xy=(phase_Ps[-1][0], 9.1), xycoords='data', ha='center', va='bottom',
                     rotation=0., annotation_clip=False, fontsize=7)
        ax1.annotate('Ps LAB', xy=(phase_Ps[-1][-1], 9.1), xycoords='data', ha='center', va='bottom',
                     rotation=0., annotation_clip=False, fontsize=7)

        if phase_Pbs:
            ax1.annotate('Pbs', xy=(phase_Pbs[-1][0], 9.1), xycoords='data', ha='center', va='bottom',
                         rotation=0., annotation_clip=False, fontsize=7)
            ax1.plot(phase_Pbs, slow, color='black')
            ax1.plot(phase_PpPbs, slow, color='crimson')
        # end if

        ax1.spines['bottom'].set_visible(False)
        ax1.tick_params(labelbottom='off')
        ax1.spines['bottom'].set_visible(False)
        ax1.yaxis.tick_right()
        ax1.yaxis.set_label_position("right")
        xlabels = ax1.get_xticklabels()
        ylabels = ax1.get_yticklabels()
        for label in xlabels:
            label.set_rotation(90)
            label.set_fontsize(7)
        # end for
        for label in ylabels:
            label.set_rotation(90)
            label.set_fontsize(7)
        # end for

        ax1.annotate(station, xy=(-0.08, 0), ha='left', va='center',
                     xycoords='axes fraction', textcoords='offset points', rotation=90.)

        start, end = ax1.get_ylim()
        ax1.yaxis.set_ticks(np.arange(start+1, end+1, 1))

        cs = ax2.contourf(xi, -1.*yi, zi, levels=levels, extend='both', cmap=cmap)
        cs.cmap.set_under('k')
        cs.cmap.set_over('k')
        ax2.spines['top'].set_visible(False)
        ax2.set_ylim(-9, -5)
        ax2.set_xlim(-5, 20)
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position("right")
        ylabels = ax2.get_yticklabels()
        ax2.plot(phase_Ps, -slow, color='black')
        ax2.plot(phase_PpPmS, -slow, color='crimson')
        ax2.plot(phase_PpSmS, -slow, color='purple')

        ax2.annotate('+PpPms', xy=(phase_PpPmS[-1][0], -9.1), xycoords='data', ha='center', va='top',
                     rotation=0., annotation_clip=False, fontsize=7, color='crimson')
        ax2.annotate('-PpSms', xy=(phase_PpSmS[-1][0], -9.1), xycoords='data', ha='center', va='top',
                     rotation=0., annotation_clip=False, fontsize=7, color='purple')

        if phase_PpPbs:
            ax2.annotate('+PpPbs', xy=(phase_PpPbs[-1][0], -9.1), xycoords='data', ha='center', va='top',
                         rotation=0., annotation_clip=False, fontsize=7, color='crimson')

            ax2.plot(phase_PpPbs, -slow, color='crimson')
            ax2.plot(phase_Pbs, -slow, color='black')
        # end if

        for label in ylabels:
            label.set_rotation(90)
            label.set_fontsize(7)
        # end for

        if frame > 3:
            xlabels = ax2.get_xticklabels()
            for label in xlabels:
                label.set_rotation(90)
                label.set_fontsize(7)
            ax2.set_xlabel('Time (sec.)')
        else:
            ax2.set_xticklabels([])
        # end if

        if (frame % 2) != 0:
            ax2.annotate('Slowness s/deg', xy=(1.2, 1), ha='left', va='center',
                         xycoords='axes fraction', textcoords='offset points', rotation=90.)
        # end if

        start, end = ax2.get_ylim()
        ax2.yaxis.set_ticks(np.arange(start, end, 1))
        traces.moveout()
        x = np.array(list(range(traces[0].data.shape[0])))*traces[0].stats.delta - 5.
        y = traces.stack()

        # Some amplitude scaling to have nice plot
        y = y[0].data/1.5 - 5.
        ax2.plot(x, y, clip_on=False, linewidth=3, color='white')
        ax2.plot(x, y, clip_on=False, linewidth=1)

        fig.add_subplot(ax1)
        fig.add_subplot(ax2)

        frame = frame + 1
        print('frame', frame)
        if frame >= rows*columns:
            cb_ax = fig.add_axes([0.25, 0.98, 0.5, 0.02])
            labels = fig.colorbar(cs, cax=cb_ax, ticks=[np.min(zi), 0, np.max(zi)], orientation='horizontal',
                                  extend='neither', extendfrac=0.00001, extendrect=True, drawedges=False)
            # labels.set_ticks([np.min(zi), 0, np.max(zi)])
            # labels.set_ticklabels(['-', '0', '+'])
            cb_ax.set_xticks([np.min(zi), 0, np.max(zi)])
            cb_ax.set_xticklabels(['-', '0', '+'])
            # labels.ax.set_yticklabels(['-', '0', '+'])
            pdf.savefig()
            figure += 1
            frame = 0
            printed = True
            # plt.show()
            plt.close()
        # end if

    # end for

    if not printed:
        cb_ax = fig.add_axes([0.25, 0.95, 0.5, 0.02])
        labels = fig.colorbar(cs, cax=cb_ax, ticks=[-1, 0, 1], orientation='horizontal', extend='neither',
                              extendfrac=0.00001, extendrect=True, drawedges=False)
        labels.ax.set_xticklabels(['-', '0', '+', ''])
        pdf.savefig()
        plt.close()
    # end if

    pdf.close()
# end main


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
# end if
