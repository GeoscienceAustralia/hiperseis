#!/usr/bin/env python
# coding: utf-8

class LayerProps():
    """
    Helper class to contain layer bulk material properties
    """

    def __init__(self, vp, vs, rho, thickness):
        """
        Constructor for given properties

        :param vp: P-wave body wave velocity
        :param vs: S-wave body wave velocity
        :param rho: Bulk material density
        :param thickness: 1D (vertical) thickness of the layer.
        """
        self._Vp = vp
        self._Vs = vs
        self._rho = rho
        self._H = thickness  # H value here is thickness of the individual layer, NOT depth relative to surface
    # end func

    @property
    def Vp(self):
        return self._Vp
    # end func

    @property
    def Vs(self):
        return self._Vs
    # end func

    @property
    def rho(self):
        return self._rho
    # end func

    @property
    def H(self):
        return self._H
    # end func

    def __repr__(self):
        return '(' + ', '.join(['Vp=' + str(self._Vp),
                                'Vs=' + str(self._Vs),
                                'ρ=' + str(self._rho),
                                'H=' + str(self._H)]) + ')'
    # end func

# end class

