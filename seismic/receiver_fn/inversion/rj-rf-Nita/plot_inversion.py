#!/usr/bin/env python
"""Plot inversion results from Bodin code
"""

import os

import click
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

paper_size_A4_landscape = (11.69, 8.27)  # inches
fig_size = (0.9*paper_size_A4_landscape[0], 0.9*paper_size_A4_landscape[1])


def plot_bodin_inversion(data_dir, rf_waveform='RF_obs.dat', pdf_outpath='.'):
    """TODO
    """
    output_file = os.path.join(pdf_outpath, 'RF_inversion_result.pdf')

    with PdfPages(output_file) as pdf:

        with open(rf_waveform, 'r') as rfin:
            rf_dat = rfin.readlines()

        with open(os.path.join(data_dir, 'data_best.out'), 'r') as rf_synth:
            synth_dat = rf_synth.readlines()

        xar_obs = []
        xar_synth = []
        yar_obs = []
        yar_synth = []

        for j in range(len(rf_dat)):
            x, y = rf_dat[j].strip('\n').split(None)
            xar_obs.append(float(x))
            yar_obs.append(float(y))
            a, b = synth_dat[j].strip('\n').split(None)
            xar_synth.append(float(a))
            yar_synth.append(float(b))

        plt.figure(figsize=fig_size)
        plt.plot(xar_obs, yar_obs, 'k', label='observed', alpha=0.8)
        plt.plot(xar_synth, yar_synth, 'r', label='synthetic', alpha=0.8)
        plt.legend()
        plt.xlabel('Time (sec)')
        plt.ylabel('RF amplitude')
        plt.title('Observed vs synthesized RF from earth model')
        aligner_color = "#a0a0a080"
        plt.gca().xaxis.grid(True, color=aligner_color, linestyle=':')

        pdf.savefig(dpi=300, papertype='a4', orientation='landscape')
        plt.close()

        # read in stuff

        with open(os.path.join(data_dir, 'Posterior.out'), 'r') as posterior:
            post_dat = posterior.readlines()
        prof, disd, d_max = post_dat[0].strip('\n').split(None)
        beta_min, beta_max, disv, width = post_dat[1].strip('\n').split(None)
        post_rest = post_dat[2:]

        with open(os.path.join(data_dir, 'Convergence_misfit.out'), 'r') as conv_misfit:
            misf_dat = conv_misfit.readlines()
        c_misfit = np.zeros([2, len(misf_dat)])
        for y in range(len(misf_dat)):
            c_misfit[0][y], c_misfit[1][y] = misf_dat[y].strip('\n').split(None)

        with open(os.path.join(data_dir, 'Convergence_nb_layers.out'), 'r') as conv_layers:
            convlay_dat = conv_layers.readlines()
        c_layers = np.zeros([2, len(convlay_dat)])
        for a in range(len(convlay_dat)):
            c_layers[0][a], c_layers[1][a] = convlay_dat[a].strip('\n').split(None)

        with open(os.path.join(data_dir, 'Convergence_sigma.out'), 'r') as conv_sigma:
            csigma_dat = conv_sigma.readlines()
        csig = np.zeros([2, len(csigma_dat)])
        for b in range(len(csigma_dat)):
            csig[0][b], csig[1][b] = csigma_dat[b].strip('\n').split(None)

        with open(os.path.join(data_dir, 'NB_layers.out'), 'r') as layers:
            layers_dat = layers.readlines()

        with open(os.path.join(data_dir, 'Sigma.out'), 'r') as sigma:
            sigmadat = sigma.readlines()
        sigmar = np.zeros([2, len(sigmadat)])
        for w in range(len(sigmadat)):
            sigmar[0][w], sigmar[1][w] = sigmadat[w].strip('\n').split(None)

        with open(os.path.join(data_dir, 'Average.out'), 'r') as average:
            ave_data = average.readlines()
        ave = np.zeros([2, len(ave_data)])
        for q in range(len(ave_data)):
            ave[0][q], ave[1][q] = ave_data[q].strip('\n').split(None)

        with open(os.path.join(data_dir, 'Change_points.out'), 'r') as cp:
            cp_data = cp.readlines()
        cpar = np.zeros([2, len(cp_data)])
        for t in range(len(cp_data)):
            cpar[0][t], cpar[1][t] = cp_data[t].strip('\n').split(None)

        # Main plot
        count = 0
        P = np.zeros([int(disd), int(disv)])
        maxx = np.zeros([int(disd), 1])
        for i in range(int(disd)):
            for j in range(int(disv)):
                P[i][j] = float((post_rest)[count].strip('\n'))
                count += 1
            maxx[i] = P[i].argmax()

        plt.figure(figsize=fig_size)
        x = [float(beta_min), float(beta_max)]
        y = [0, float(prof)]

        maxx = (maxx - 0.5) * ((float(beta_max) - float(beta_min)) / float(disv)) + float(beta_min)

        v = 31
        n = (v * float(disd) / float(prof)) + 0.5

        plt.subplot(131)
        vmax = np.max(P[P != P[0][0]])
        # Recommended linear colormaps are plasma, inferno or viridis
        plt.imshow(P, cmap='plasma', extent=[x[0], x[1], y[1], y[0]], vmin=0, vmax=vmax, aspect='auto')
        plt.xlabel('Vs (km/s)')
        plt.ylabel('Depth (km)')
        x_range = [float(beta_min), float(beta_max)]
        plt.xlim(x_range)
        y_range = [0, float(prof)]
        plt.ylim(y_range)
        plt.title('Candidate solutions')
        plt.gca().invert_yaxis()

        plt.subplot(132)
        plt.plot([float(beta_min), (float(beta_max) - 2 * float(width))], [0, float(d_max)], 'k', alpha=0.8)
        plt.plot([float(beta_min) + 2 * float(width), float(beta_max)], [0, float(d_max)], 'k', alpha=0.8)
        plt.plot(maxx, ave[0], 'k', label='Most probable', alpha=0.8)
        plt.plot(ave[1], ave[0], 'r', label='Mean', alpha=0.8)
        plt.xlabel('Vs (km/s)')
        plt.xlim(x_range)
        plt.ylim(y_range)
        plt.legend()
        plt.title('Velocity profile')
        plt.gca().invert_yaxis()

        plt.subplot(133)
        plt.plot(cpar[1] / max(cpar[1]), cpar[0], 'k')
        plt.fill_betweenx(cpar[0], cpar[1] / max(cpar[1]), 0, color='darkgrey')
        plt.xlim([0, 1])
        plt.ylim(y_range)
        plt.title('Confidence')
        plt.xlabel('P(transition)')
        plt.gca().invert_yaxis()

        pdf.savefig(dpi=300, papertype='a4', orientation='landscape')
        plt.close()

        # Histograms
        plt.figure(figsize=fig_size)

        plt.subplot(211)
        plt.bar(np.array(range(len(layers_dat))), height=(np.array(layers_dat, dtype='int')) /
                                                         sum(np.array(layers_dat, dtype='float')), width=1)
        plt.xlim([1, len(layers_dat)])
        plt.xlabel('# layers')
        plt.ylabel('P(# layers)')
        plt.title('Layer likelihood distribution')

        plt.subplot(212)
        plt.bar(sigmar[0], height=(np.array(sigmar[1], dtype='float') / sum(np.array(sigmar[1], dtype='float'))),
                width=(sigmar[0][1] - sigmar[0][0]))
        plt.xlim([sigmar[0][0], sigmar[0][-1]])
        plt.xlabel(r'$\sigma$ for RF')
        plt.ylabel(r'p($\sigma$)')

        pdf.savefig(dpi=300, papertype='a4', orientation='landscape')
        plt.close()

        # Convergence plots...
        plt.figure(figsize=fig_size)

        plt.subplot(311)
        plt.plot(c_misfit[0][1:], label='variable 0', alpha=0.8)
        plt.plot(c_misfit[1][1:], label='variable 1', alpha=0.8)
        # plt.plot([c_misfit[0][0],0],[c_misfit[0][0],c_misfit[1][0]],'r')
        plt.ylabel('Misfit')  # unit = ??
        plt.legend()
        plt.grid(color=aligner_color, linestyle=':')

        plt.subplot(312)
        plt.plot(c_layers[0][1:], label='variable 0', alpha=0.8)
        plt.plot(c_layers[1][1:], label='variable 0', alpha=0.8)
        plt.ylabel('# layers')
        plt.grid(color=aligner_color, linestyle=':')

        plt.subplot(313)
        plt.plot(csig[0][1:], label='variable 0', alpha=0.8)
        plt.plot(csig[1][1:], label='variable 0', alpha=0.8)
        plt.xlabel('Iteration number')
        plt.ylabel(r'$\sigma$')
        plt.grid(color=aligner_color, linestyle=':')

        plt.suptitle('Convergence history')

        pdf.savefig(dpi=300, papertype='a4', orientation='landscape')
        plt.close()
    # end with pdf
# end func


# -------------Main---------------------------------

@click.command()
@click.option('--rf-waveform', type=click.Path(exists=True, dir_okay=False), required=True,
              help=r'Input Receiver Function waveform used for inversion (.dat file)')
@click.option('--input-folder', type=click.Path(exists=True, dir_okay=True, file_okay=False), required=True,
              help=r'Folder containing the input files for generating PDF report. Same as output folder of inversion.')
@click.argument('pdf-output-folder', type=click.Path(exists=True, dir_okay=True, file_okay=False), required=True)
def main(rf_waveform, input_folder, pdf_output_folder):
  plot_bodin_inversion(input_folder, rf_waveform, pdf_output_folder)
# end func


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
# end if
