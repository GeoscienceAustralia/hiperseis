import numpy as np
import traceback

class Grid:
    def __init__(self, min_lon, max_lon,
                 min_colat, max_colat,
                 nx, ny,
                 depth_km_list,
                 depth_refinement_factor=1,
                 block_index_offset=0):
        assert ((min_lon >= 0) & (min_lon <= 360) & \
                (max_lon >= 0) & (max_lon <= 360) & \
                (min_colat >= 0) & (min_colat <= 190) & \
                (max_colat >= 0) & (max_colat <= 180) & \
                (min_lon < max_lon) & (min_colat < max_colat)), \
                'Ensure lons[0, 360], colats[0, 180], min_lon < max_lon and min_lat < max_lat'
        assert ((depth_refinement_factor >= 1) & ((depth_refinement_factor & (depth_refinement_factor - 1) == 0))), \
               'Invalid z_refinement_factor; must be a power of 2 and >= 1'

        self.min_lon = min_lon
        self.max_lon = max_lon
        self.min_colat = min_colat
        self.max_colat = max_colat
        self.depth_km_list = np.array(depth_km_list)
        self.nx = int(nx)
        self.ny = int(ny)
        self.depth_refinement_factor = int(depth_refinement_factor)
        self.block_index_offset = block_index_offset

        self.dx = (self.max_lon - self.min_lon) / (self.nx - 1)
        self.dy = (self.max_colat - self.min_colat) / (self.ny - 1)

        if (self.depth_refinement_factor > 1):
            rf = self.depth_refinement_factor

            while (rf > 1):
                curr = self.depth_km_list
                refined = np.zeros((len(curr) - 1) * 2 + 1)

                refined[0::2] = curr
                refined[1::2] = curr[:-1] + (curr[1:] - curr[:-1]) / float(rf)

                self.depth_km_list = refined
                rf /= 2
            # wend
        # end if

        self.depth_spans = np.vstack([self.depth_km_list[0:-1],
                                      self.depth_km_list[1:]]).T
        self.nd = len(self.depth_km_list)
        self.ncells = (self.nx - 1) * (self.ny - 1) * (self.nd - 1)
    # end func

    def _is_within_grid(self, lon, colat, depth_km):
        return (lon >= self.min_lon) & (lon <= self.max_lon) & \
               (colat >= self.min_colat) & (colat <= self.max_colat) & \
               (depth_km >= self.depth_km_list[0]) & (depth_km <= self.depth_km_list[-1])
    # end func

    def get_block_index(self, lon, colat, depth_km):
        assert ((lon >= 0) & (lon <= 360) & \
                (colat >= 0) & (colat <= 180)) & \
                (depth_km <= self.depth_km_list[-1]), \
                'Ensure lon[0, 360], colat[0, 180], depth_km[{}, {}]'.format(self.depth_km_list[0],
                                                                         self.depth_km_list[-1])

        if (self._is_within_grid(lon, colat, depth_km)):
            k = np.argwhere((depth_km >= self.depth_spans[:, 0]) & \
                            (depth_km < self.depth_spans[:, 1]))
            if (k.size):
                k = k[0][0]
            else:
                k = self.nd - 2  # d-cell indices [0, nd-2]

            i = np.min([int((lon - self.min_lon) / self.dx), self.nx - 2])  # x-cell indices [0, nx-2]
            j = np.min([int((colat - self.min_colat) / self.dy), self.ny - 2])  # y-cell indices [0, ny-2]

            bidx = k * (self.nx - 1) * (self.ny - 1) + j * (self.nx - 1) + i + 1
            bidx += self.block_index_offset

            return bidx
        # end if

        return -1
    # end func

    def __str__(self):
        r = ['Grid: {}'.format(hex(id(self)))]
        for k in self.__dict__.keys():
            r.append('{}: {}'.format(k, self.__dict__[k]))
        # end for
        return '\n'.join(r)
    # end func
# end class

class NestedGrid:
    def __init__(self, param_file, outer_depth_refinement_factor=1, inner_depth_refinement_factor=1):
        # parse parameterization
        self.param_file = param_file

        try:
            fh = open(param_file)
            # read number of outer grid cells in x, y and d
            oncx, oncy, oncd, odx, ody = map(float, fh.readline().strip().split())
            self.o_depth_km_list = []
            for i in np.arange(int(oncd + 1)):
                self.o_depth_km_list.append(float(fh.readline().strip()))
            # end for

            # read extents of inner grid and number of cells in x, y and d
            imin_lon, imax_lon, imin_lat, imax_lat, incx, incy, incd = map(float, fh.readline().strip().split())
            imin_colat = 90 - imax_lat
            imax_colat = 90 - imin_lat

            self.i_depth_km_list = []
            for i in np.arange(int(incd + 1)):
                self.i_depth_km_list.append(float(fh.readline().strip()))
            # end for

            # instantiate inner and outer grids
            self.ig = Grid(imin_lon, imax_lon, imin_colat, imax_colat,
                           incx + 1, incy + 1, self.i_depth_km_list, inner_depth_refinement_factor)
            self.og = Grid(0, 360, 0, 180, oncx + 1, oncy + 1, self.o_depth_km_list,
                           outer_depth_refinement_factor, block_index_offset=self.ig.ncells)

            assert self.og.dx == odx, 'Prescribed x-resolution of outer grid does not match that inferred..'
            assert self.og.dy == ody, 'Prescribed y-resolution of outer grid does not match that inferred..'
        except Exception as e:
            print('Error: failed to read param file {} ({})..'.format(self.param_file,
                                                                      traceback.format_exc()))
            raise RuntimeError
        # end try
    # end func

    def get_block_index(self, lon, colat, depth_km):
        bidx = self.ig.get_block_index(lon, colat, depth_km)
        if (bidx < 0):
            bidx = self.og.get_block_index(lon, colat, depth_km)
        # end if

        return bidx
    # end func

    def is_inner_block(self, block_index:int):
        return block_index <= self.og.block_index_offset
    # end func

    def __str__(self):
        r = ['NestedGrid: {}'.format(hex(id(self)))]
        for k in self.__dict__.keys():
            r.append('{}: {}'.format(k, self.__dict__[k]))
        # end for
        return '\n'.join(r)
    # end func
# end class

