#! /usr/bin/env python
"""
Description:
    Define a 3D grid discretization of the Earth

CreationDate:
    2/09/2018

Developer:
    | fei.zhang@ga.gov.au
"""

import os
import logging
import numpy as np
import pandas as pd
# from deprecated import deprecated

log = logging.getLogger(__name__)


class Grid2:
    """
    A non-uniform Grid Model of Earth for eventsource-->station rays clustering and sorting.

    Reference to the inversion Fortran code model parameter file param1x1 in this directory:
        - Global: 5x5 degree  and a list of 1+16 depths
        - Region: 1x1 degree  and a list of 1+29 depths
    """

    # def __init__(self, nx=360, ny=180, dz=10000):
    def __init__(self, ndis=2, param_file=None):
        """
        constructor for a non-uniform Grid.
        Please modify the parameters in this method according to your desired Grid attributes.
        :param ndis: 2, 4 (define the ref discretization by 2 or 4)
        """
        self.ndis = ndis

        if param_file is None:
            print("Use the default harded-coded Grid parameters in ", __file__)
            self._define_params()
            # self.show_properties()
        else:
            # get params from a file
            self._parse_parametrisation(param_file)

            # Further re-scaling definition
            self.dx = self.dx / self.ndis
            self.dy = self.dy / self.ndis
            self.gdx = self.gdx / self.ndis
            self.gdy = self.gdy / self.ndis

            self.nx = self.nx * self.ndis
            self.ny = self.ny * self.ndis
            self.gnx = self.gnx * self.ndis
            self.gny = self.gny * self.ndis

            self.refrmeters = self._refine_depth(self.rmeters, ndis=self.ndis)

            self.refgmeters = self._refine_depth(self.gmeters, ndis=self.ndis)

            self.REGION_MAX_BN = self._get_max_region_block_number()

        return


    def _define_params(self):
        """
        define the grid parameters
        :param ndis:
        :return:
        """

        # Region bounadry definition
        self.LON = (100.0, 190.0)
        self.LAT = (-54.0, 0.0)

        # region grid size in lat-long plane
        self.nx = self.ndis * 90  # 1/4 degree cell size in LONG
        self.ny = self.ndis * 54  # 1/4 degree cell size in LAT
        self.dx = (self.LON[1] - self.LON[0]) / self.nx  # 1/ndis degree
        self.dy = (self.LAT[1] - self.LAT[0]) / self.ny  # 1/ndis degree

        # Region's depths in KM according to Fortran param1x1 (read in from file)
        self.rdepth = (0, 10, 35, 70, 110, 160, 210, 260, 310, 360, 410, 460, 510, 560,
                       610, 660, 710, 810, 910, 1010, 1110, 1250, 1400, 1600)

        self.nz = len(self.rdepth) -1

        # Region Depth in meters
        self.rmeters = 1000 * np.array(self.rdepth)

        self.refrmeters = self._refine_depth(self.rmeters, ndis=self.ndis)

        # global grid model params (outside the region)
        self.gdx = 5 * self.dx  # =360/self.gnx;  5 times as big as the regional grid size
        self.gdy = 5 * self.dy  # =180/self.gny;  5 times as big as the regional grid size

        self.gnx = int(360 / self.gdx)  # 360/2.5 =144
        self.gny = self.gnx / 2  # =77

        # Global's depths in KM according to Fortran param1x1
        self.gdepth = (0, 110, 280, 410, 660, 840, 1020, 1250, 1400, 1600, 1850, 2050, 2250, 2450, 2600, 2750, 2889)
        self.gnz = len(self.gdepth) -1

        self.gmeters = 1000 * np.array(self.gdepth)  # Region Depth in meters

        # self.gdz = 20000  # global outside the ANZ zone 20KM uniform cellsize in depth
        self.refgmeters = self._refine_depth(self.gmeters)

        self.REGION_MAX_BN = self._get_max_region_block_number()

        return

    def _parse_parametrisation(self, param_file):
        """
        Reads in the local and global grid parameters from an external
        parametrisation file (by default we load  param_file='./example_param.txt')
        """

        with open(param_file, 'r') as F:
            # Rewind to file beginning
            F.seek(0)

            ## PART 1 - GLOBAL GRID
            # Read-in the global grid (number of cells and resolution)
            gnx, gny, gnz, gdx, gdy = F.readline().strip('\n').split(' ')
            self.gnx = int(gnx)
            self.gny = int(gny)
            self.gnz = int(gnz)
            self.gdx = float(gdx)
            self.gdy = float(gdy)

            # There is currently no option to reduce the coverage of the global
            # grid, so resolution information is a little redundant. Nevertheless,
            # let's perform a sanity test.
            assert (360. / self.gdx) == self.gnx
            assert (180. / self.gdy) == self.gny

            # Load global depth nodes
            self.gdepth = np.zeros(self.gnz + 1)
            for i in range(self.gnz + 1):
                try:
                    depth = float(F.readline().strip('\n'))
                except ValueError as err:
                    print(('There appear to be less global depth levels than ' + \
                           'was indicated! The header suggested %d levels, but ' + \
                           'I can only parse %d. Please check the param file.') % (self.gnz + 1, i))
                    raise err
                self.gdepth[i] = depth

            # Global depth in meters
            self.gmeters = 1000. * self.gdepth

            ## PART 2 - LOCAL GRID
            # Read-in the local grid (extent and number of cells)
            try:
                LON_min, LON_max, LAT_min, LAT_max, nx, ny, nz = F.readline().strip('\n').split(' ')
                LON_min = float(LON_min)
                LON_max = float(LON_max)
                LAT_min = float(LAT_min)
                LAT_max = float(LAT_max)
                self.LON = (LON_min, LON_max)
                self.LAT = (LAT_min, LAT_max)
                self.nx = int(nx)
                self.ny = int(ny)
                self.nz = int(nz)
            except ValueError as err:
                print(('There appear to be more global depth levels than ' + \
                       'was indicated! The header suggested %d levels; ' + \
                       'please check the param file.') % (self.nz + 1))
                raise err

            # Local grid resolution
            self.dx = (self.LON[1] - self.LON[0]) / self.nx
            self.dy = (self.LAT[1] - self.LAT[0]) / self.ny

            # Load local depth nodes
            self.rdepth = np.zeros(self.nz + 1)
            for i in range(self.nz + 1):
                try:
                    depth = float(F.readline().strip('\n'))
                except ValueError as err:
                    print(('There appear to be less local depth levels than ' + \
                           'was indicated! The header suggested %d levels, but ' + \
                           'I can only parse %d. Please check the param file.') % (self.nz + 1, i))
                    raise err
                self.rdepth[i] = depth

            # Local depth in meters
            self.rmeters = 1000. * self.rdepth

            #print(self.rdepth)
            return

    def _refine_depth(self, dep_meters, ndis=2):
        """
        define a refined depths discretization, 1/ndis size of the depth steps of input np array dep_meters
        :param dep_meters: numpy array listing of the dpeth in inversion model
        :param ndis: refined discretization, 2, 4
        :return: an np array of depths
        """

        refdepth = np.arange((dep_meters.size - 1) * ndis + 1)
        # print(refdepth)
        for i in range(refdepth.size - 1):
            i0 = int(i / ndis)
            i1 = i % ndis
            # print(i0, i1)
            refdepth[i] = dep_meters[i0] + (dep_meters[i0 + 1] - dep_meters[i0]) * i1 * (1.0 / ndis)

        refdepth[-1] = dep_meters[-1]  # last depth is same

        # print(refdepth)  # check the smaller-step discretization of depth
        return refdepth

    def _get_max_region_block_number(self):
        """
        Compute the MAX BLOCK NUMBER for a point in the region box (with a finer grid)
        This max block number is used to offset bn for points in the global zone.
        assuming the event's max depth is 2000KM

        :return: int max_nb
        """
        # assume 2000KM max depth, compute the z-layers layer_depth
        layer_depth = 5000
        k = round(2000000.0 / layer_depth) + 1
        bn = k * self.nx * self.ny + 1  # must make sure this number is big enough, to separate region|global

        log.info("The REG_MAX_BN = %s", bn)

        return bn

    def is_point_in_region(self, lat, lon):
        """
        test if the event or station point is in the region box?
        
        :param lat: latitude
        :param lon: longitude
        :return: boolean T or F
        """

        # if (abs(lat) > 90):
        #    log.error("wrong lattitude value %s", lat)

        y = (lat + 90) % 180 - 90.0  # convert lat into y which will be in [0,180) then to (-90,90])

        x = lon % 360  # convert longitude into x which must be in [0,360)

        if (lat <= self.LAT[1] and lat >= self.LAT[0] and x <= self.LON[1] and x >= self.LON[0]):
            log.debug("(%s, %s) in the region", lat, lon)
            return True
        else:
            log.debug("(%s, %s) NOT in the region ", lat, lon)
            return False

    def find_block_number(self, lat, lon, z):
        """
        find the index-number of an events/station in a non-uniform grid
        each spatial point (lat, lon, z) is mapped to a uniq block_number.

        :param lat: latitude (-90,90)
        :param lon: longitude (0,180)
        :param z: depth
        :return: int block number AND (xc,yc, zcm)
        """
        x = lon % 360  # convert lon into x which must be in [0,360)
        y = (lat + 90.0) % 180  # convert lat into y which will be in [0,180)

        if self.is_point_in_region(lat, lon) is True:

            i = round(x / self.dx) + 1
            j = round(y / self.dy) + 1

            # k = round(z / self.dz) + 1
            (k, zcm) = self.get_depth_index(z, self.refrmeters)

            block_number = (k - 1) * self.nx * self.ny + (j - 1) * self.nx + i

            xc = ((i - 1) + 0.5) * self.dx  # cell block center longitude in deg
            yc = ((j - 1) + 0.5) * self.dy  # cell block centre lattitude in deg

            # xc = (i + 0.5) * self.dx  # cell block center longitude in deg
            # yc = (j + 0.5) * self.dy  # cell block centre lattitude in deg
            zc = zcm  # cell block centre depth in meters

        else:  # point NOT in ANZ region, but in the global

            i = round(x / self.gdx) + 1
            j = round(y / self.gdy) + 1

            # k = round(z / self.gdz) + 1 # this is a constant steps
            (k, zcm) = self.get_depth_index(z, self.refgmeters)
            g_block_number = (k - 1) * self.gnx * self.gny + (j - 1) * self.gnx + i
            block_number = g_block_number + self.REGION_MAX_BN
            # global block_number is offset by the regional max block number, to ensure uniqueness

            xc = ((i - 1) + 0.5) * self.gdx  # cell block center longitude in deg
            yc = ((j - 1) + 0.5) * self.gdy  # cell block centre lattitude in deg
            # xc = (i  + 0.5) * self.gdx  # cell block center longitude in deg
            # yc = (j  + 0.5) * self.gdy  # cell block centre lattitude in deg
            zc = zcm  # cell block centre depth in meters

        yc = yc % 180 - 90.0  # Lattitude from (0,180) back to [-90,90) could be 91.25 => -88.75
        xc = xc % 360  # map the longitude xc to the [0,360)

        # return int(block_number)
        return (int(block_number), xc, yc, zc)  # return block_number and the block's centre coordinate xc,yc,zc)

    def get_depth_index(self, z, dep_meters):
        """
        given a point with depth z meters, and a np-array of refined depth in meters,
        find the index which z fall into.

        :param z: an event depth in meters
        :param dep_meters: an array of numbers corresponding to a refined depth discretization.
        :return:  depth index and the cell block centre depth in metres.
        """

        diff_refdepth = dep_meters - z

        # print (diff_refdepth)
        kindex = np.where(diff_refdepth > 0)[0]

        if kindex.size > 0:
            # print ("the max kindex is %s with positive diff value %s" % (kindex[0], diff_refdepth[kindex[0]]))
            k = kindex[0]
        else:
            k = diff_refdepth.size - 1  # the last index

        zcm = (dep_meters[k] + dep_meters[k - 1]) / 2.0

        assert k >= 1  # k should be 1, 2,...

        return (k, zcm)

    def __str__(self):
        """
        String representaiton of the object
        https://stackoverflow.com/questions/4932438/how-to-create-a-custom-string-representation-for-a-class-object
        :return:
        """
        print("****** Show all attributes of the Grid instance:  ******* ")
        mykeys = self.__dict__.keys()

        sorted_keys = sorted(mykeys)

        mylist=[]
        for k in sorted_keys:
            #print('%-10s' % (str(k)) + ':', my_grid.__dict__[k])
            mylist.append("%s : %s "%(str(k), my_grid.__dict__[k]))

        return str(mylist)

    def show_properties(self):
        """
        print the properties of the grid definition

        """

        print(str(self))  # print a string representation of this object

        print("**********************************")
        print("The refined discretization divided by = ", self.ndis)
        print("Reginal grid prop= ", self.dx, self.dy, self.nx, self.ny)
        print("Gloabl grid prop= ", self.gdx, self.gdy, self.gnx, self.gny)

        print("The MAX BLOCK Number inside ANZ region", self.REGION_MAX_BN)

        print("The refined regional depth steps and size:", self.refrmeters, self.refrmeters.size)
        print("The refined global depth steps and size:", self.refgmeters, self.refgmeters.size)

        dmeters = 10000  # an event depth in meters
        dmeters = 10001  # +1 in event depth to cross cell boundary|
        dmeters = 0

        k, zcm = self.get_depth_index(dmeters, self.refrmeters)

        print("Regional: The k index = %s" % k)
        print("Regional: The zcm  = %s" % zcm)

        k, zcm = self.get_depth_index(dmeters, self.refgmeters)

        print("Global: The k index = %s" % k)
        print("Global: The zcm  = %s" % zcm)

        # Test some special border points
        print("Global:", self.find_block_number(-90, 0, 1))
        print("Global:", self.find_block_number(-89, 0, 1))

        print("regional:", self.find_block_number(-1, 101, 1))
        print("regional:", self.find_block_number(-2, 102, 1))

        return

    def generate_latlong_grid(self, depthmeters=0.0):
        """
        create a csv file containing: (block_number, lat,long, depthm=0, xc,yc,zc)

        :return: Pandas data frame
        """

        # how many points per degree in lat and lon respectively
        per_degree = 10  # 5 means 0.2 degree cell size
        latgrid = np.linspace(-90, 90, num=180 * per_degree + 1)
        longrid = np.linspace(0, 360, num=360 * per_degree + 1)

        mygrid = []
        for xin in longrid:
            for yin in latgrid:
                (bn, xc, yc, zc) = self.find_block_number(yin, xin, depthmeters)
                arow = (bn, xc, yc, zc, xin, yin, depthmeters)
                mygrid.append(arow)

        mypdf = pd.DataFrame(mygrid, columns=['blockn', 'xc', 'yc', 'zc', 'long', 'lat', 'depthmeter'])

        # mypdf.to_csv("/tmp/whole_grid.csv", index=False) # write CSV without the sequence index column 01,2,3,

        final_pdf = mypdf.drop_duplicates(subset=['blockn'], keep='first', inplace=False)

        return final_pdf

    def generate_grid3D(self):
        """
        Loop over all the reference depths, both regional and global depth list.

        :return: pandas dataframe
        """

        pdf3d = self.generate_latlong_grid()

        for adepth in self.refrmeters[1:10]:
            print("Making grid at depth:", adepth)
            apdf = self.generate_latlong_grid(depthmeters=adepth)
            pdf3d = pdf3d.append(apdf)
            print("The size of the csv=", pdf3d.shape)

        for adepth in self.refgmeters[1:]:
            print("Making grid at depth:", adepth)
            apdf = self.generate_latlong_grid(depthmeters=adepth)
            pdf3d = pdf3d.append(apdf)
            print("The size of the csv=", pdf3d.shape)

        print("The final size of the pandas dataframe to be saved into CSV file: ", pdf3d.shape)

        pdf3d.to_csv("/tmp/cluster_grid_3d.csv", index=False)

        return pdf3d

# @deprecated(reason="Use the class Grid2")
class UniformGrid:
    """
    This is the original (simple) uniform grid model for Earth.

    It is replaced by the more generic non-uniform model defined above as class Grid2
    """

    def __init__(self, nx=360, ny=180, dz=10000):
        """
        constructor for Grid. Define 3D cell properties
        :param nx: integer number of cells in longitude (360)
        :param ny: lattitude cells (180)
        :param dz: depth in meteres (default=10KM)
        """
        self.nx = nx
        self.ny = ny
        self.dx = 360.0 / nx
        self.dy = 180.0 / ny
        self.dz = dz

    def find_block_number(self, lat, lon, z):
        """
        find the 3D-block number in this grid

        :param lat: lattitude
        :param lon: longitude
        :param z: elevation
        :return: int block number
        """
        # y = 90. - lat  # original wrong
        y = (lat + 90) % 180  # y will be in (0,180)
        x = lon % 360
        i = round(x / self.dx) + 1
        j = round(y / self.dy) + 1
        k = round(z / self.dz) + 1
        block_number = (k - 1) * self.nx * self.ny + (j - 1) * self.nx + i

        return int(block_number)


# ========================================================
# quick test run to see if the Grid definition is right??
# ========================================================
if __name__ == "__main__":

    # my_grid = Grid2(ndis=4)
    # my_grid = Grid2(param_file='/g/data/ha3/fxz547/Githubz/passive-seismic/seismic/traveltime/example_param.txt')
    my_grid = Grid2(param_file='/g/data/ha3/fxz547/Githubz/passive-seismic/seismic/traveltime/param1x1')
    # my_grid = Grid2(param_file='/g/data/ha3/fxz547/Githubz/passive-seismic/seismic/traveltime/param2x2')
    # my_grid = Grid2()
    my_grid.show_properties()

    # Generate the whole grid:
    # mypdf = my_grid.generate_latlong_grid()
    # print (mypdf.head())

    # my3dgrid = my_grid.generate_grid3D()  # this may take a long time
    # print("Final size of the 3D grid", my3dgrid.shape)
