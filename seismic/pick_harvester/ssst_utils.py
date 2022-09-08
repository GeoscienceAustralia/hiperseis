import numpy as np
from pyproj import Geod
from tqdm import tqdm
from collections import defaultdict
import h5py

def h5_to_named_array(hfn, key):
    h = h5py.File(hfn, 'r')

    result = np.array([])
    try:
        if (type(h[key]) == h5py.Group):
            names = []
            formats = []
            data = []
            for k in h[key].keys():
                names.append(k)
                var = h[key][k][:]
                data.append(var)
                formats.append(var.dtype.str)
            # end for

            fields = {'names': names, 'formats': formats}
            result = np.zeros(len(data[0]), dtype=fields)

            for item, k in zip(data, names):
                result[k] = item
            # end for
        else:
            result = h[key][:]
        # end if
    except Exception as e:
        print(str(e))
    # end try
    h.close()

    return result
# end func

def get_iters(hfn):
    h = h5py.File(hfn, 'r')

    iters = list(map(int, h.keys()))
    h.close()
    return iters
# end func

class SSST_Result:
    def __init__(self, h5_fn, iteration=-1):
        self.h5_fn = h5_fn
        self.iter = get_iters(self.h5_fn)[iteration]
        self.geod = Geod(a=180 / np.pi, f=0)

        self.arrivals = h5_to_named_array(self.h5_fn, '{}/arrivals'.format(self.iter))
        self.events = h5_to_named_array(self.h5_fn, '{}/events'.format(self.iter))
        self.event_id_to_idx = h5_to_named_array(self.h5_fn, '0/events/event_id_to_idx')
        self.residual = h5_to_named_array(self.h5_fn, '{}/arrivals/residual'.format(self.iter))
        self.tcorr = h5_to_named_array(self.h5_fn, '{}/arrivals/tcorr'.format(self.iter))
        self.is_AUTO_arrival = h5_to_named_array(self.h5_fn, '0/arrivals/is_AUTO_arrival')

        self.elons = self.events['lon'][self.event_id_to_idx[self.arrivals['event_id']]]
        self.elats = self.events['lat'][self.event_id_to_idx[self.arrivals['event_id']]]
        self.edepths_km = self.events['depth_km'][self.event_id_to_idx[self.arrivals['event_id']]]
        self.eorigin_ts = self.events['origin_ts'][self.event_id_to_idx[self.arrivals['event_id']]]

        self.phase = self.arrivals['phase']
        self.slons = self.arrivals['lon']
        self.slats = self.arrivals['lat']
        self.selevs_km = self.arrivals['elev_m'] / 1e3
        self.arrival_ts = self.arrivals['arrival_ts']

        _, _, self.ecdists = self.geod.inv(self.elons, self.elats, self.slons, self.slats)
        equality = h5_to_named_array(self.h5_fn, '{}/events/event_quality'.format(self.iter))

        self.equality = equality[self.event_id_to_idx[self.arrivals['event_id']]]
        self.is_P = self.phase.astype('S1') == b'P'
        self.is_S = self.phase.astype('S1') == b'S'
        self.slope_ratio = self.arrivals['quality_measure_slope']
    # end func

    def find_paired(self):
        paired = np.zeros(len(self.arrivals), dtype='?')
        has_p_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int))))
        indices = np.arange(len(self.arrivals))

        # Generate a dict that groups all P arrivals by eid, net, sta, loc
        for idx in tqdm(indices[self.is_P], desc='Grouping P-arrivals: '):
            eid, net, sta, loc = self.arrivals[idx]['event_id'], self.arrivals[idx]['net'], \
                                 self.arrivals[idx]['sta'], self.arrivals[idx]['loc']
            has_p_dict[eid][net][sta][loc] = 1
        # end for

        # Label each S arrival that has a corresponding P arrival
        for idx in tqdm(indices[self.is_S], desc='Labelling S-arrivals: '):
            eid, net, sta, loc = self.arrivals[idx]['event_id'], self.arrivals[idx]['net'], \
                                 self.arrivals[idx]['sta'], self.arrivals[idx]['loc']

            if (has_p_dict[eid][net][sta][loc]): paired[idx] = 1
        # end for

        return paired
    # end func
# end class
