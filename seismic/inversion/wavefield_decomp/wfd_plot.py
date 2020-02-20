#!/usr/bin/env python
# coding: utf-8
"""
Plotting helper functions for Wavefield Decomposition module.
"""

import numpy as np

import matplotlib.pyplot as plt

# pylint: disable=invalid-name


def plot_Esu_space(H, k, Esu, title=None, savefile_name=None, show=True, c_range=(None, None)):
    """
    Plot SU energy as function of H-k.

    :param H: Grid of H coordinates corresponding to Esu point values.
    :param k: Grid of k coordinates corresponding to Esu point values.
    :param Esu: Energy values on H-k grid.
    :param title: Plot title [OPTIONAL]
    :param savefile_name: Output file in which to save plot [OPTIONAL]
    :param show: If True, display the image and block until it is closed.
    :param c_range: Custom range of Esu to contour (min, max values)
    :return: None
    """
    colmap = 'plasma'
    plt.figure(figsize=(16, 12))
    min_E, max_E = (np.nanmin(Esu) if c_range[0] is None else c_range[0],
                    np.nanmax(Esu) if c_range[1] is None else c_range[1])
    plt.contourf(k, H, Esu, levels=np.linspace(min_E, max_E, 51), cmap=colmap)
    plt.colorbar()
    plt.contour(k, H, Esu, levels=np.linspace(min_E, max_E, 11), colors='k', linewidths=1, antialiased=True)
    plt.xlabel('Crustal $\kappa$', fontsize=14)
    plt.ylabel('Crustal $H$ (km)', fontsize=14)
    plt.tick_params(right=True, labelright=True, which='both')
    plt.tick_params(top=True, labeltop=True, which='both')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.minorticks_on()
    plt.xlim(np.min(k), np.max(k))
    plt.ylim(np.min(H), np.max(H))
    plt.grid(linestyle=':', color="#80808080")
    if title is not None:
        plt.title(title, fontsize=20, y=1.05)
    # end if
    if savefile_name is not None:
        plt.savefig(savefile_name, dpi=300)
    # end if
    if show:
        plt.show()
    # end if
    plt.close()
# end func
