# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 14:15:44 2018

@author: u81234
"""
import numpy as np


class grid3():
    def __init__(self):
        self.parse_parametrisation()
        
#        self.dz = 5000  # assumed uniform cellsize in depth inside the ANZ zone 5KM
#        self.refrmeters = self._refine_depth(self.rmeters, ndis=self.ndis)
#        
#        # self.gdz = 20000  # global outside the ANZ zone 20KM uniform cellsize in depth
#        self.refgmeters = self._refine_depth(self.gmeters)
#
#        self.REGION_MAX_BN = self._get_max_region_block_number()
#
#        return
    
    def parse_parametrisation(self, filepath='./param'):
        """
        Reads in the local and global grid parameters from an external 
        parametrisation file (by default we load './param')
        """
        with open(filepath, 'r') as F:
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
            self.gdepth = np.zeros(self.gnz+1)
            for i in range(self.gnz+1):
                try:
                    depth = float(F.readline().strip('\n'))
                except ValueError as err:
                    print(('There appear to be less global depth levels than ' + \
                           'was indicated! The header suggested %d levels, but ' + \
                           'I can only parse %d. Please check the param file.') % (self.gnz+1, i))
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
                       'please check the param file.') % (self.nz+1))
                raise err
            
            # Local grid resolution
            self.dx = (self.LON[1] - self.LON[0]) / self.nx
            self.dy = (self.LAT[1] - self.LAT[0]) / self.ny
            
            # Load local depth nodes
            self.rdepth = np.zeros(self.nz+1)
            for i in range(self.nz+1):
                try:
                    depth = float(F.readline().strip('\n'))
                except ValueError as err:
                    print(('There appear to be less local depth levels than ' + \
                           'was indicated! The header suggested %d levels, but ' + \
                           'I can only parse %d. Please check the param file.') % (self.nz+1, i))
                    raise err
                self.rdepth[i] = depth
            
            # Local depth in meters
            self.rmeters = 1000. * self.rdepth
    

if __name__=='__main__':
    G = grid3()
    L = G.__dict__.keys()
    L.sort()
    for k in L:
        print('%-10s' % (str(k)) + ':', G.__dict__[k])
    