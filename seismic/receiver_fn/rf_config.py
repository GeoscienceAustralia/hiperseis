"""
Description:
    Parser for config file used in generating RFs

References:

CreationDate:   13/12/23
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     13/12/23   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

from collections import defaultdict
import copy
import json

DEFAULT_RF_TYPE = 'prf'  # from ['prf', 'srf']
DEFAULT_RESAMPLE_RATE_HZ = 10.0
DEFAULT_FILTER_BAND_HZ = (0.02, 1.00)
DEFAULT_TAPER_LIMIT = 0.05
DEFAULT_TRIM_START_TIME_SEC = -50.0
DEFAULT_TRIM_END_TIME_SEC = 150.0
DEFAULT_ROTATION_TYPE = 'zrt'   # from ['zrt', 'lqt']
DEFAULT_DECONV_DOMAIN = 'time'  # from ['time', 'freq', 'iter']
DEFAULT_GAUSS_WIDTH = 1.0
DEFAULT_WATER_LEVEL = 0.01
DEFAULT_SPIKING = 0.5
DEFAULT_NORMALIZE = False
DEFAULT_PLOT_DIR = None

class RFConfig:
    def __init__(self, config_file: str):
        """
        Config file consists of 3 sub-dictionaries. One named "filtering" for
        input stream filtering settings, one named "processing" for RF processing
        settings, and one named "correction" for rotating/swapping/negating channel
        data for one or more named stations with potential orientation discrepancies.
        Each of these sub-dicts is described below::

        "filtering":  # Filtering settings
        {
          "resample_rate": float # Resampling rate in Hz
          "taper_limit": float   # Fraction of signal to taper at each end, between 0 and 0.5
          "filter_band": (float, float) # Filter pass band (Hz). Not required for freq-domain deconvolution.
        }

        "processing":  # RF processing settings
        {
          "rf_type": str # Choice of ['prf, srf']
          "custom_preproc":
          {
            "import": 'import custom symbols',  # statement to import required symbols
            "func": 'preproc functor'  # expression to get handle to custom preprocessing functor
            "args": {}  # additional kwargs to pass to func
          }
          "trim_start_time": float # Trace trim start time in sec, relative to onset
          "trim_end_time": float # Trace trim end time in sec, relative to onset
          "rotation_type": str # Choice of ['zrt', 'lqt']. Rotational coordinate system
                               # for aligning ZNE trace components with incident wave direction
          "deconv_domain": str # Choice of ['time', 'freq', 'iter']. Whether to perform deconvolution
                               # in time or freq domain, or iterative technique
          "gauss_width": float # Gaussian freq domain filter width. Only required for freq-domain deconvolution
          "water_level": float # Water-level for freq domain spectrum. Only required for freq-domain deconvolution
          "spiking": float # Spiking factor (noise suppression), only required for time-domain deconvolution
          "normalize": bool # Whether to normalize RF amplitude
        }

        "correction": # corrections to be applied to data for named stations prior to RF computation
        {
          "plot_dir": str # path to folder where plots related to orientation corrections are to be saved
          "swap_ne": list # list of NET.STA.LOC for which N and E channels are to be swapped, e.g ["OA.BL27."],
          "rotate": list # list of NET.STA.LOC that are to be rotated to maximize P-arrival energy on \
                           the primary RF component, e.g ["OA.BL27."]
          "negate": list # list of NET.STA.LOC.CHA that are to be negated, e.g ["OA.BL27..HHZ"]
          "recompute_inclinations": # list of NET.STA.LOC for which the inclination angles in the LQT reference
                                      frame are to be recomputed by minimizing the energy in the L component.
                                      Note that this parameter is only applicable for S-RFs.
        }

        @param config_file: name of config file
        """
        # define default config
        self.default_config = \
            {
                "filtering":
                    {
                        "resample_rate": DEFAULT_RESAMPLE_RATE_HZ,
                        "taper_limit": DEFAULT_TAPER_LIMIT,
                        "filter_band": DEFAULT_FILTER_BAND_HZ,
                    },

                "processing":
                    {
                        "rf_type": DEFAULT_RF_TYPE,
                        "trim_start_time": DEFAULT_TRIM_START_TIME_SEC,
                        "trim_end_time": DEFAULT_TRIM_END_TIME_SEC,
                        "rotation_type": DEFAULT_ROTATION_TYPE,
                        "deconv_domain": DEFAULT_DECONV_DOMAIN,
                        "gauss_width": DEFAULT_GAUSS_WIDTH,
                        "water_level": DEFAULT_WATER_LEVEL,
                        "spiking": DEFAULT_SPIKING,
                        "normalize": DEFAULT_NORMALIZE,
                    },

                "correction":
                    {
                        "plot_dir": DEFAULT_PLOT_DIR,
                        "swap_ne": [],
                        "rotate": [],
                        "negate": [],
                        "recompute_inclinations": [],
                    }
            }

        # load user config
        self.config_file = config_file
        if(self.config_file is not None):
            try:
                with open(self.config_file, 'r') as cf:
                    self.config = json.load(cf)
                # end with
            except Exception as e:
                raise RuntimeError('Failed to read config file: {}. Aborting..'.format(self.config_file))
            # end try
        else:
            self.config = {}
        # end if

        self.config_filtering = {}
        self.config_processing = {}
        self.config_correction = {}

        self._merge()
        self._validate()
    # end func

    def _merge(self):
        """
        Merges default config with the user-provided config
        @return:
        """
        def selective_merge(base_obj, delta_obj):
            if not isinstance(base_obj, dict):
                return delta_obj
            # end if
            common_keys = set(base_obj).intersection(delta_obj)
            new_keys = set(delta_obj).difference(common_keys)
            for k in common_keys:
                base_obj[k] = selective_merge(base_obj[k], delta_obj[k])
            for k in new_keys:
                base_obj[k] = delta_obj[k]
            return base_obj
        # end func

        self.config = selective_merge(copy.deepcopy(self.default_config), self.config)

        self.config_filtering = self.config.get("filtering")
        self.config_processing = self.config.get("processing")
        self.config_correction = self.config.get("correction")
    # end func

    def _validate(self):
        # validate <filtering> block
        cf_keys = {'resample_rate', 'taper_limit', 'filter_band'}
        if(not set(self.config_filtering.keys()).issubset(cf_keys)):
            raise ValueError('Invalid key(s) found in <filtering> block in the config file. '
                             'Valid keys are: {}'.format(cf_keys))
        # end if

        # validate <processing> block
        cp_keys = {'rf_type', 'custom_preproc', 'trim_start_time', 'trim_end_time',
                   'rotation_type', 'deconv_domain', 'gauss_width', 'water_level',
                   'spiking', 'normalize'}
        if(not set(self.config_processing.keys()).issubset(cp_keys)):
            raise ValueError('Invalid key(s) found in <processing> block in the config file. '
                             'Valid keys are: {}'.format(cp_keys))
        # end if
        rf_type = self.config_processing.get('rf_type', DEFAULT_RF_TYPE)
        rotation_type = self.config_processing.get('rotation_type', DEFAULT_ROTATION_TYPE)
        deconv_domain = self.config_processing.get('deconv_domain', DEFAULT_DECONV_DOMAIN)

        if(rf_type not in ['prf', 'srf']):
            raise ValueError("Invalid rf_type. Must be one of ['prf', 'srf']")
        if(rotation_type.lower() not in ['zrt', 'lqt']):
            raise ValueError("Invalid rotation_type. Must be one of ['zrt', 'lqt']")
        if(deconv_domain not in ['time', 'freq', 'iter']):
            raise ValueError("Invalid deconv_domain. Must be one of ['time', 'freq', 'iter']")

        if(rf_type is 'srf' and rotation_type.lower() is not 'lqt'):
            raise ValueError('Rotation-type for srf must be lqt.')
        # end if

        # validate <correction> block
        cc_keys = {'plot_dir', 'swap_ne', 'rotate', 'negate', 'recompute_inclinations'}
        if (not set(self.config_correction.keys()).issubset(cc_keys)):
            raise ValueError('Invalid key(s) found in <correction> block in the config file. '
                             'Valid keys are: {}'.format(cc_keys))
        # end if
    # end func
# end class

class Corrections:
    def __init__(self, config: RFConfig, hdf_keys: list):
        self.config = config

        config_correction = self.config.config_correction

        self._swap_ne_list = config_correction.setdefault('swap_ne', [])
        self._rotate_list = config_correction.setdefault('rotate', [])
        self._negate_list = config_correction.setdefault('negate', [])
        self._recompute_inclinations_list = config_correction.setdefault('recompute_inclinations', [])
        self._plot_dir = config_correction.setdefault('plot_dir', None)
        self._corrections = defaultdict(lambda: None)

        if (type(self._swap_ne_list) != list): assert 0, \
            'config:correction: Expected a list of net.sta.loc for key: swap_ne'
        if (type(self._rotate_list) != list): assert 0, \
            'config:correction: Expected a list of net.sta.loc for key: rotate'
        if (type(self._negate_list) != list): assert 0, \
            'config:correction: Expected a list of net.sta.loc.cha for key: negate'
        if (type(self._recompute_inclinations_list) != list): assert 0, \
            'config:correction: Expected a list of net.sta.loc.cha for key: recompute_inclinations'

        # NE channel swaps
        try:
            for item in self._swap_ne_list:
                net, sta, loc = item.split('.')
                found = False
                for hdf_item in hdf_keys:
                    if ('.'.join([net, sta, loc]) == hdf_item):
                        if (not self._corrections['swap_ne']): self._corrections['swap_ne'] = defaultdict(bool)
                        self._corrections['swap_ne'][hdf_item] = True
                        found = True
                        break
                    # end if
                # end for
                if (not found):
                    assert 0, 'Station {} in config:correction:swap_ne block not ' \
                              'found in input dataset'.format(item)
                # end if
            # end for
        except Exception as e:
            print(e)
            assert 0, 'An invalid net.sta.loc item was found in config:correction:rotate block'
        # end try

        # Rotations
        try:
            # sanity check
            if(len(self._rotate_list) and self.config.config_processing['rf_type'] == 'srf'):
                raise RuntimeError('Rotations are only applicable for P-RFs.')
            # end if

            for item in self._rotate_list:
                net, sta, loc = item.split('.')
                found = False
                for hdf_item in hdf_keys:
                    if ('.'.join([net, sta, loc]) == hdf_item):
                        if (not self._corrections['rotate']): self._corrections['rotate'] = defaultdict(bool)
                        self._corrections['rotate'][hdf_item] = True
                        found = True
                        break
                    # end if
                # end for
                if (not found):
                    assert 0, 'Station {} in config:correction:rotate block not ' \
                              'found in input dataset'.format(item)
                # end if
            # end for
        except Exception as e:
            print(e)
            assert 0, 'An invalid net.sta.loc item was found in config:correction:rotate block'
        # end try

        # Negations
        try:
            for item in self._negate_list:
                net, sta, loc, cha = item.split('.')
                found = False
                for hdf_item in hdf_keys:
                    if ('.'.join([net, sta, loc]) == hdf_item):
                        if (not self._corrections['negate']): self._corrections['negate'] = defaultdict(list)
                        self._corrections['negate'][hdf_item].append = cha[-1]
                        found = True
                        break
                    # end if
                # end for
                if (not found):
                    assert 0, 'Station {} in config:correction:negate block not ' \
                              'found in input dataset'.format(item)
                # end if
            # end for
        except Exception as e:
            print(e)
            assert 0, 'An invalid net.sta.loc.cha item was found in config:correction:negate block'
        # end try

        # Recompute inclinations
        try:
            # sanity check
            if(len(self._recompute_inclinations_list) and \
                    self.config.config_processing['rf_type'] == 'prf'):
                raise RuntimeError('Inclinations can be recomputed for S-RFs only')
            # end if

            for item in self._recompute_inclinations_list:
                net, sta, loc, cha = item.split('.')
                found = False
                for hdf_item in hdf_keys:
                    if ('.'.join([net, sta, loc]) == hdf_item):
                        if (not self._corrections['recompute_inclinations']):
                            self._corrections['recompute_inclinations'] = defaultdict(bool)
                        self._corrections['recompute_inclinations'][hdf_item] = True
                        found = True
                        break
                    # end if
                # end for
                if (not found):
                    assert 0, 'Station {} in config:correction:negate block not ' \
                              'found in input dataset'.format(item)
                # end if
            # end for
        except Exception as e:
            print(e)
            assert 0, 'An invalid net.sta.loc.cha item was found in config:correction:recompute_inclinations block'
        # end try
    # end func

    def needsChannelSwap(self, netstaloc):
        if(self._corrections['swap_ne']):
            return self._corrections['swap_ne'][netstaloc]
        else:
            return False
        # end if
    # end func

    def needsRotation(self, netstaloc):
        if(self._corrections['rotate']):
            return self._corrections['rotate'][netstaloc]
        else:
            return False
        # end if
    # end func

    def needsNegation(self, netstaloc):
        if(self._corrections['negate']):
            result = bool(len(self._corrections['negate'][netstaloc]))

            if(result): return self._corrections['negate'][netstaloc]
            else: return False
        else:
            return False
    # end func

    def needsInclinationRecomputed(self, netstaloc):
        if(self._corrections['recompute_inclinations']):
            return self._corrections['recompute_inclinations'][netstaloc]
        else:
            return False
    # end func

    def needsCorrections(self, netstaloc):
        return self.needsChannelSwap(netstaloc) or \
               self.needsRotation(netstaloc) or \
               self.needsNegation(netstaloc) or \
               self.needsInclinationRecomputed(netstaloc)
    # end func

    @property
    def plot_dir(self):
        return self._plot_dir
    # end func
# end class

if __name__ == "__main__":
    rfc = RFConfig('config.json')
    pass
# end if