#!/usr/bin/env python
# coding: utf-8
"""Helper classes to encapsulate model properties.
"""

# pylint: disable=invalid-name


class LayerProps():
    """
    Helper class to contain layer bulk material properties
    """

    def __init__(self, vp, vs, rho, thickness):
        """
        Constructor for given properties

        :param vp: P-wave body wave velocity
        :type vp: float
        :param vs: S-wave body wave velocity
        :type vs: float
        :param rho: Bulk material density
        :type rho: float
        :param thickness: 1D (vertical) thickness of the layer.
        :type thickness: float
        """
        self._Vp = vp
        self._Vs = vs
        self._rho = rho
        self._H = thickness  # H value here is thickness of the individual layer, NOT depth relative to surface
    # end func

    @property
    def Vp(self):
        """Get P-wave body wave velocity
        """
        return self._Vp
    # end func

    @property
    def Vs(self):
        """Get S-wave body wave velocity
        """
        return self._Vs
    # end func

    @property
    def rho(self):
        """Get bulk material density
        """
        return self._rho
    # end func

    @property
    def H(self):
        """Get layer thickness
        """
        return self._H
    # end func

    def __repr__(self):
        return '(' + ', '.join(['Vp=' + str(self._Vp),
                                'Vs=' + str(self._Vs),
                                'ρ=' + str(self._rho),
                                'H=' + str(self._H)]) + ')'
    # end func

# end class
