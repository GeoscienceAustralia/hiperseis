#!/usr/bin/env python
"""Plot inversion results from Bodin code
"""

import os

import click
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
from scipy import signal

paper_size_A4_landscape = (11.69, 8.27)  # inches
fig_size = (0.9*paper_size_A4_landscape[0], 0.9*paper_size_A4_landscape[1])


def _load_convergence_sequence(filename, downsampling=1):
    d = pd.read_csv(filename, header=None, skiprows=1)
    d = d[0].values
    # We downsample without filtering for simplicity, since only long term trend is of much interest.
    return d[::downsampling]
# end func


def plot_bodin_inversion(pdf, data_dir, rf_waveform, station=''):
    """Standardized plotting of the results of Bodin RF inversion code.
    """
    convergence_ds = 1000  # downsampling rate for convergence sequences

    rf_dat = np.loadtxt(rf_waveform)
    synth_dat = np.loadtxt(os.path.join(data_dir, 'data_best.out'))

    xar_obs = rf_dat[:,0]
    yar_obs = rf_dat[:,1]
    xar_synth = synth_dat[:,0]
    yar_synth = synth_dat[:,1]

    # Plot actual vs predicted RF
    station_id = ' (' + station + ')' if station else ''
    plt.figure(figsize=fig_size)
    plt.plot(xar_obs, yar_obs, 'k', label='observed', alpha=0.8)
    plt.plot(xar_synth, yar_synth, 'r', label='synthetic (from rank 0)', alpha=0.8)
    plt.legend()
    plt.xlabel('Time (sec)')
    plt.ylabel('RF amplitude')
    plt.title('Observed vs synthesized RF from earth model' + station_id)
    aligner_color = "#a0a0a080"
    plt.gca().xaxis.grid(True, color=aligner_color, linestyle=':')
    plt.gca().xaxis.set_minor_locator(MultipleLocator(1))

    pdf.savefig(dpi=300, papertype='a4', orientation='landscape')
    plt.close()

    # read in other model results

    with open(os.path.join(data_dir, 'Posterior.out'), 'r') as posterior:
        post_dat = posterior.readlines()
    prof, disd, d_max = post_dat[0].strip('\n').split(None)
    disd = int(disd)
    beta_min, beta_max, disv, width = post_dat[1].strip('\n').split(None)
    disv = int(disv)
    post_rest = np.reshape(np.array([float(x.strip('\n')) for x in post_dat[2:]]), (disd, disv))

    c_misfit = _load_convergence_sequence(os.path.join(data_dir, 'Convergence_misfit.out'), convergence_ds)
    c_layers = _load_convergence_sequence(os.path.join(data_dir, 'Convergence_nb_layers.out'), convergence_ds)
    csig = _load_convergence_sequence(os.path.join(data_dir, 'Convergence_sigma.out'), convergence_ds)

    layers_dat = np.loadtxt(os.path.join(data_dir, 'NB_layers.out'), dtype=int)
    sigmar = np.loadtxt(os.path.join(data_dir, 'Sigma.out'))
    ave = np.loadtxt(os.path.join(data_dir, 'Average.out'))
    cpar = np.loadtxt(os.path.join(data_dir, 'Change_points.out'))

    # Main plot
    count = disd*disv
    P = post_rest
    maxx = np.argmax(P, axis=1)
    plt.figure(figsize=fig_size)
    x = [float(beta_min), float(beta_max)]
    y = [0, float(prof)]

    maxx = (maxx - 0.5) * ((float(beta_max) - float(beta_min)) / float(disv)) + float(beta_min)

    v = 31
    n = (v * float(disd) / float(prof)) + 0.5

    plt.subplot(131)
    vmax = np.max(P[2:, 2:])
    # Recommended linear colormaps are plasma, inferno or viridis
    plt.imshow(P, cmap='plasma', extent=[x[0], x[1], y[1], y[0]], vmin=0, vmax=vmax, aspect='auto')
    plt.xlabel('Vs (km/s)')
    plt.ylabel('Depth (km)')
    x_range = [float(beta_min), float(beta_max)]
    plt.xlim(x_range)
    y_range = [0, float(prof)]
    plt.ylim(y_range)
    plt.title('Ensemble PDF' + station_id)
    plt.gca().invert_yaxis()
    plt.gca().yaxis.set_minor_locator(MultipleLocator(2))

    plt.subplot(132)
    plt.plot([float(beta_min), (float(beta_max) - 2 * float(width))], [0, float(d_max)], 'k', alpha=0.8)
    plt.plot([float(beta_min) + 2 * float(width), float(beta_max)], [0, float(d_max)], 'k', alpha=0.8)
    plt.plot(maxx, ave[:,0], 'k', label='Most probable', alpha=0.8)
    plt.plot(ave[:,1], ave[:,0], 'r', label='Mean', alpha=0.8)
    plt.xlabel('Vs (km/s)')
    plt.xlim(x_range)
    plt.ylim(y_range)
    plt.legend()
    plt.title('Velocity profile')
    plt.gca().invert_yaxis()
    plt.gca().yaxis.set_minor_locator(MultipleLocator(2))

    plt.subplot(133)
    plt.plot(cpar[:,1] / max(cpar[:,1]), cpar[:,0], 'k')
    plt.fill_betweenx(cpar[:,0], cpar[:,1] / max(cpar[:,1]), 0, color='darkgrey')
    plt.xlim([0, 1])
    plt.ylim(y_range)
    # Prevalance of change points at this depth among the model ensemble (confidence in delta-V)
    plt.title('Prevalence of change points')
    plt.xlabel('P(transition)')
    plt.gca().invert_yaxis()
    plt.gca().yaxis.set_minor_locator(MultipleLocator(2))

    pdf.savefig(dpi=300, papertype='a4', orientation='landscape')
    plt.close()

    # Histograms
    plt.figure(figsize=fig_size)

    plt.subplot(211)
    plt.bar(np.array(range(len(layers_dat))), height=layers_dat/float(layers_dat.sum()), width=1)
    plt.xlim([1, len(layers_dat)])
    plt.xlabel('# layers')
    plt.ylabel('P(# layers)')
    plt.title('Layer likelihood distribution' + station_id)

    plt.subplot(212)
    plt.bar(sigmar[:,0], height=sigmar[:,1]/float(sigmar[:,1].sum()),
            width=(sigmar[1,0] - sigmar[0,0]))
    plt.xlim([sigmar[0,0], sigmar[-1,0]])
    plt.xlabel(r'$\sigma$ for RF')
    plt.ylabel(r'p($\sigma$)')

    pdf.savefig(dpi=300, papertype='a4', orientation='landscape')
    plt.close()

    # Convergence plots...
    plt.figure(figsize=fig_size)

    plt.subplot(311)
    it_num = np.arange(len(c_misfit))*convergence_ds
    plt.semilogy(it_num, c_misfit, label='Mean', alpha=0.8)
    plt.ylabel('Misfit')  # unit = ??
    plt.legend()
    plt.grid(color=aligner_color, linestyle=':')

    plt.subplot(312)
    it_num = np.arange(len(c_layers))*convergence_ds
    plt.plot(it_num, c_layers, label='Mean', alpha=0.8)
    plt.ylabel('# layers')
    plt.grid(color=aligner_color, linestyle=':')

    plt.subplot(313)
    it_num = np.arange(len(csig))*convergence_ds
    plt.semilogy(it_num, csig, label='Mean', alpha=0.8)
    plt.xlabel('Iteration number')
    plt.ylabel(r'$\sigma$')
    plt.grid(color=aligner_color, linestyle=':')

    plt.suptitle('Convergence history' + station_id, y=0.92)

    pdf.savefig(dpi=300, papertype='a4', orientation='landscape')
    plt.close()

# end func


# -------------Main---------------------------------

@click.command()
@click.option('--rf-waveform', type=click.Path(exists=True, dir_okay=False), required=True,
              help=r'Input Receiver Function waveform used for inversion (.dat file)')
@click.option('--input-folder', type=click.Path(exists=True, dir_okay=True, file_okay=False), required=True,
              help=r'Folder containing the input files for generating PDF report. Same as output folder of inversion.')
@click.option('--station-id', type=str, help='Station identification string', default='', show_default=True)
@click.argument('pdf-output-folder', type=click.Path(exists=True, dir_okay=True, file_okay=False), required=True)
def main(rf_waveform, input_folder, pdf_output_folder, station_id=''):
    station_nodot = station_id.replace('.', '_') + '_' if station_id else ''
    output_file = os.path.join(pdf_output_folder, station_nodot + 'RF_inversion_result.pdf')
    with PdfPages(output_file) as pdf:
        plot_bodin_inversion(pdf, input_folder, rf_waveform, station=station_id)
    # end with
# end func


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
# end if
