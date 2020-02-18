#!/usr/bin/env python
# coding: utf-8
"""
Plotting helper functions for Wavefield Decomposition module.
"""

import numpy as np

import matplotlib.pyplot as plt

# pylint: disable=invalid-name


def plot_Esu_space(H, k, Esu, title=None, savefile_name=None, show=True):
    colmap = 'plasma'
    plt.figure(figsize=(16, 12))
    plt.contourf(k, H, Esu, levels=50, cmap=colmap)
    plt.colorbar()
    plt.contour(k, H, Esu, levels=10, colors='k', linewidths=1, antialiased=True)
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
