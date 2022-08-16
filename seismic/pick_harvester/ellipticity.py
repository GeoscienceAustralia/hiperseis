from collections import defaultdict
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import traceback
import os

class Phase:
    def __init__(self, ecdist, depth, tau0, tau1, tau2):
        self.ecdist = ecdist
        self.depth = depth
        self.tau0 = tau0
        self.tau1 = tau1
        self.tau2 = tau2

        self.tau0io = RegularGridInterpolator((self.ecdist, self.depth), self.tau0, method='linear',
                                              bounds_error=False, fill_value=0.)
        self.tau1io = RegularGridInterpolator((self.ecdist, self.depth), self.tau1, method='linear',
                                              bounds_error=False, fill_value=0.)
        self.tau2io = RegularGridInterpolator((self.ecdist, self.depth), self.tau2, method='linear',
                                              bounds_error=False, fill_value=0.)
    # end func
# end class

class Ellipticity:
    def __init__(self):
        self._elcordir_tbl = os.path.dirname(os.path.abspath(__file__)) + '/ellip/elcordir.tbl'
        self.depth = np.array([0.0, 100.0, 200.0, 300.0, 500.0, 700.0], dtype='f4')
        self.phases = defaultdict(list)
        self.phase_alias = { b'Pg': b'Pup',
                             b'Sg': b'Sup',
                             b'pPg': b'pP',
                             b'sPg': b'sP',
                             b'pSg': b'pS',
                             b'sSg': b'sS',
                             b'Pb': b'P',
                             b'Sb': b'S',
                             b'pPb': b'pP',
                             b'sPb': b'sP',
                             b'pSb': b'pS',
                             b'sSb': b'sS',
                             b'Pn': b'P',
                             b'Sn': b'S',
                             b'pPn': b'pP',
                             b'sPn': b'sP',
                             b'pSn': b'pS',
                             b'sSn': b'sS',
                             b'SPn': b'SP',
                             b'SPb': b'SP',
                             b'SPg': b'SP',
                             b'SnP': b'SP',
                             b'PSn': b'PS',
                             b'PnPn': b'PP',
                             b'SnSn': b'SS',
                             b'p': b'Pup',
                             b's': b'Sup',
                             b'Pdif': b'Pdiff',
                             b'Sdif': b'Sdiff' }
        self._parse_table()
    # end func

    def _parse_table(self):
        line = open(self._elcordir_tbl).read()
        ll = len(line)

        indices = np.arange(0, ll+80, 80)
        rows = list()

        for i in range(len(indices)-1):
            row = line[indices[i]:indices[i+1]].split()
            rows.append(row)
        #end for

        phase = None
        ecdist = list()
        tau = 0
        tau0 = list()
        tau1 = list()
        tau2 = list()
        for irow, row in enumerate(rows):
            if len(row) == 4:
                if irow:
                    ecdist = np.array(ecdist).astype('f4')
                    tau0 = np.array(tau0).astype('f4')
                    tau1 = np.array(tau1).astype('f4')
                    tau2 = np.array(tau2).astype('f4')
                    self.phases[phase] = Phase(ecdist, self.depth, tau0, tau1, tau2)
                    ecdist = list()
                    tau0 = list()
                    tau1 = list()
                    tau2 = list()
                #end if
                phase = (row[0]).encode() # convert to byte-string
            elif len(row) == 1:
                ecdist.append(float(row[0]))
            elif len(row) == 6 and tau == 0:
                tau0.append(row)
                tau = (tau + 1) % 3
            elif len(row) == 6 and tau == 1:
                tau1.append(row)
                tau = (tau + 1) % 3
            elif len(row) == 6 and tau == 2:
                tau2.append(row)
                tau = (tau + 1) % 3
            #end if
        #end for

        for pnew, palias in self.phase_alias.items():
            if(palias in self.phases.keys()): self.phases[pnew] = self.phases[palias]
        # end for
    # end func

    def get_correction(self, phase, ecdist, depth_km, elat, azim):
        azim = np.radians(azim)
        ecolat = np.radians(90 - elat)
        s3 = np.sqrt(3.0) / 2.0

        try:
            if(type(phase) == np.ndarray):
                result = np.zeros(len(phase), dtype='f4')

                corr_count = 0
                for cphase in self.phases.keys():
                    indices = np.argwhere(cphase == phase).flatten()

                    if(len(indices)):
                        sc0 = 0.25 * (1.0 + 3.0 * np.cos(2.0 * ecolat[indices]))
                        sc1 = s3 * np.sin(2.0 * ecolat[indices])
                        sc2 = s3 * np.sin(ecolat[indices]) * np.sin(ecolat[indices])

                        tau0 = self.phases[cphase].tau0io((ecdist[indices], depth_km[indices]))
                        tau1 = self.phases[cphase].tau1io((ecdist[indices], depth_km[indices]))
                        tau2 = self.phases[cphase].tau2io((ecdist[indices], depth_km[indices]))

                        result[indices] = sc0*tau0 + sc1*np.cos(azim[indices])*tau1 + \
                                          sc2*np.cos(2.0*azim[indices])*tau2
                        corr_count += len(indices)
                    # end if
                # end for
                #if(corr_count < len(phase)): print('Warning: some phases {} not found..'.format(set(list(phase)) - set(self.phases.keys())))

                return result
            else:
                result = 0.
                sc0 = 0.25 * (1.0 + 3.0 * np.cos(2.0 * ecolat))
                sc1 = s3 * np.sin(2.0 * ecolat)
                sc2 = s3 * np.sin(ecolat) * np.sin(ecolat)

                tau0 = self.phases[phase].tau0io((ecdist, depth_km))
                tau1 = self.phases[phase].tau1io((ecdist, depth_km))
                tau2 = self.phases[phase].tau2io((ecdist, depth_km))

                result = sc0 * tau0 + sc1 * np.cos(azim) * tau1 + \
                         sc2 * np.cos(2.0 * azim) * tau2

                return result
            # end if
        except Exception as e:
            print(traceback.format_exc())
            raise ValueError('Phase not found..')
        # end try
    # end func
# end class
