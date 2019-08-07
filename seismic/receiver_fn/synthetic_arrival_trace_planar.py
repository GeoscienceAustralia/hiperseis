#!/usr/bin/env python
"""Use a basic planar, 2-layer model of only the crust and the Moho to generate
synthetic arrival traces for known model characteristics. Intended to be used
for model validation.
"""

import scipy


moho_depth = 40  # km
# Top layer is layer zero, layer below Moho is layer 1.
V_p = [6.2, 6.4]  # km/s, one value per layer
V_s = [3.8, 4.0]  # km/s, one value per layer
f_s = 50.0  # Hz
src_pulse = square_wave_from_scipy
