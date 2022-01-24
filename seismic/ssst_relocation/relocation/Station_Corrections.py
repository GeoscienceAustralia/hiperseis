# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 14:08:43 2021

@author: U37509
"""

import numpy as np

class Station(object): 
    def __init__(self, name, picks):
        """
        Class to facilitate computation of source specific station terms 
        for a single station. Creates a Station_phase_attr object for each
        phase associated with arrivals at the station.
        
        
        Parameters
        ----------
        name : string
            Name of station.
            
        picks : numpy.ndarray
            Structured array with the following data type.
            dtype = [('event_id', 'U20'), ('pick_id', 'U30'), ('stat', 'U10'), 
                     ('net', 'U5'), ('cha', 'U10'), ('elon', 'single'), 
                     ('ecolat', 'single'), ('edepth', 'single'), 
                     ('origin_time', 'double'), ('mag', 'half'), 
                     ('slon', 'single'), ('scolat', 'single'), 
                     ('selev', 'half'), ('phase', 'U8'), 
                     ('arrival_time', 'double'), ('ptt', 'single'), 
                     ('tcor', 'half'), ('residual', 'half'), ('snr', 'half'), 
                     ('qualityMeasureCWT', 'half'), ('domFreq', 'half'), 
                     ('qualityMeasureSlope', 'half'), ('bandIndex', 'uint8'), 
                     ('nSigma', 'uint8')]
            
        
        """
        self.name = name
        self.picks = picks
        self.phases = list(np.unique(picks['phase']))
        if 'X' in self.phases:
            self.phases.remove('X')
        #end if
        if 'Px' in self.phases:
            self.phases.remove('Px')
        #end if
        if 'Sx' in self.phases:
            self.phases.remove('Sx')
        #end if
        for phase in self.phases: setattr(self, phase, Station_phase_attr())
    #end func
#end class
    
class Station_phase_attr(): 
    def __init__(self):
        """
        This is designed to be an attribute to the Station object. One is 
        created for each phase present in the arrival at the station.
        
        During SSST algorithm, (lon, colat, depth) of each pick is placed into
        the 'picks' attribute, and (lon, colat, depth, residual) of each pick
        for which the origin quality is good are saved in the 'resids' 
        attribute.
        
        
        """
        self.resids = []
        self.picks = []
    #end func
#end class
    
def ang_dist(lon1, colat1, lon2, colat2, units='degrees'):
    """
    Function to calculate the angular distance from (lon1, colat1) to 
    (lon2, colat2).
    
    
    Parameters
    ----------
    lon1 : float
        Longitude of point 1.
        
    colat1 : float
        Colatitude of point 1.
        
    lon2 : float
        Longitude of point 2.
    
    colat2 : float
        Colatitude of point 2.
        
    units : string (optional)
        'degrees' or 'radians' describing the units in which lon1, lat1, lon2,
        lat2 are input. Default is 'degrees'.
        
        
    Returns
    -------
    value : float
        Angular distance from (lon1, colat1) to (lon2, colat2).
        
        
    """
    if units == 'degrees':
        lon1 = lon1*np.pi/180
        colat1 = colat1*np.pi/180
        lon2 = lon2*np.pi/180
        colat2 = colat2*np.pi/180
    else:
        colat1 = np.pi/2 - colat1
        colat2 = np.pi/2 - colat2
    #end if
    value = 2*np.arcsin(np.sqrt(np.sin((colat1 - colat2)/2)**2 + \
                                np.sin((lon1 - lon2)/2)**2* \
                                np.sin(colat1)*np.sin(colat2)))
    if units == 'degrees': value = value*180/np.pi
    return value
#end func

def calculate_station_corrections(statnames, picks, rank, config,
                                  unstable_events=list()):
    """
    Function to compute time corrections for picks.
    For each station in 'statnames', a 'Station' object is created and all
    corresponding picks are added to it. Then, for each phase associated with
    arrivals at the station, a SST or SSST correction is computed, and added
    to the list 'correction_list'.
    
    
    Parameters
    ----------
    statnames : list
        List of station names.
        
    picks : numpy.ndarray
        Structured array with the following data type.
        dtype = [('event_id', 'U20'), ('pick_id', 'U30'), ('stat', 'U10'), 
                 ('net', 'U5'), ('cha', 'U10'), ('elon', 'single'), 
                 ('ecolat', 'single'), ('edepth', 'single'), 
                 ('origin_time', 'double'), ('mag', 'half'), 
                 ('slon', 'single'), ('scolat', 'single'), ('selev', 'half'), 
                 ('phase', 'U8'), ('arrival_time', 'double'), 
                 ('ptt', 'single'), ('tcor', 'half'), ('residual', 'half'), 
                 ('snr', 'half'), ('qualityMeasureCWT', 'half'), 
                 ('domFreq', 'half'), ('qualityMeasureSlope', 'half'), 
                 ('bandIndex', 'uint8'), ('nSigma', 'uint8')]
    
    rank : integer
        Rank of the in-use processor.
        
    config : configparser.SectionProxy object
        Information from config file.
       
    unstable_events : list
        List of event IDs which have poor origin quality. If an event moves
        more than a threshold distance between iterations of the relocation 
        algorithm, it is classed as unstable.
        
    
    Returns
    -------
    picks : numpy.ndarray
        Input pick array with added travel time corrections.
        
        
    """
    method = config['correction_method']
    thr = float(config['corr_thr_dist_deg'])
    
    statlist = [Station(stat, picks[picks['stat'] == stat]) \
                for stat in statnames]
    correction_list = list()
    deg2km = 6371*np.pi/180
    
    for stat in statlist:
        for pick in stat.picks:
            phase = pick['phase']
            event = pick['event_id']
            if phase not in ['Px', 'Sx', 'X']:
                obj = getattr(stat, phase) 
                if event not in unstable_events:
                    obj.resids.append([pick['elon'], pick['ecolat'], 
                                       pick['edepth'], pick['residual']])
                    obj.picks.append([pick['pick_id'], pick['elon'], 
                                      pick['ecolat'], pick['edepth']])
                else:
                    obj.picks.append([pick['pick_id'], pick['elon'], 
                                      pick['ecolat'], pick['edepth']])
                #end if
                setattr(stat, phase, obj)
            #end if
        #end for
        for phase in stat.phases:
            obj = getattr(stat, phase)
            obj.resids = np.array(obj.resids)
            obj.picks = np.array(obj.picks)
            [x1, y1, z1] = np.transpose(obj.picks[:, 1:4]).astype(float)
            z1 = z1/(1e3*deg2km)
            pickIds = obj.picks[:, 0]
            [x2, y2, z2] = np.transpose(obj.resids[:, 0:3])
            z2 = z2/(1e3*deg2km)
            f = obj.resids[:, 3]
            if method == 'SSST':
                obj.method = 'SSST'
                sph_corr = sphere_corr(x1, y1, z1, x2, y2, z2, f, thr)
                obj.correction = [(pickIds[i], sph_corr[i]) \
                                  for i in range(len(pickIds))]
            elif method == 'SST':
                obj.method = 'SST'
                med_corr = np.median(f)
                obj.correction = [(pickIds[i], med_corr) \
                                  for i in range(len(pickIds))]
            else:
                obj.method = 'None'
                obj.correction = [(pickIds[i], 0.0) \
                                  for i in range(len(pickIds))]
            #end if
            correction_list.extend(obj.correction)
            setattr(stat, phase, obj)
        #end for
    #end for
    correction_arr = np.array(correction_list,
                              dtype=[('pick_id', 'U30'), ('tcor', 'float16')])
    correction_arr = np.sort(correction_arr, order='pick_id')
    picks = np.sort(picks, order='pick_id')
    ind = np.isin(picks['pick_id'], correction_arr['pick_id'])
    picks['tcor'][ind] = correction_arr['tcor']
    return picks
#end func
    
def sphere_corr(x1, y1, z1, x2, y2, z2, f, thr, minpoints=5):
    """
    Computes a time correction for each event at location (x1, y1, z1) based on
    the travel time residuals for neighbouring events at (x2, y2, z2), if 
    angular distance between epicentres is less than 'thr' and if the minimum
    number of points within the boundary sphere is more than 'minpoints'.
    
    
    Parameters
    ----------
    x1 : numpy.ndarray
        Longitude.
        
    y1 : numpy.ndarray
        Colatitude.
        
    z1 : numpy.ndarray
        Depth.
        
    x2 : numpy.ndarray
        Longitude.
        
    y2 : numpy.ndarray
        Colatitude.
        
    z2 : numpy.ndarray
        Depth.
        
    f : numpy.ndarray
        Travel time residual.
        
    thr : float
        Maximum angular distance (degrees) between events to be considered 
        'neighbours'.
        
    minpoints : integer
        Minimum number of points in neighbourhood of an event to compute a 
        time correction.
        
    
    Returns
    -------
    medians : numpy.ndarray
        Median travel time residual for events within 'thr' distance of each
        event at points (x1, y1, z1).
        
        
    """
    
    r = (ang_dist(np.expand_dims(x1, 1), np.expand_dims(y1, 1), 
                  np.expand_dims(x2, 1).T, np.expand_dims(y2, 1).T))
    
    indices = [np.argwhere(r[i, :] < thr) for i in range(len(f))]
    samples = [f[indices[i]] for i in range(len(f))]
    medians = np.array([np.median(samples[i]) if len(samples[i]) >= minpoints \
                        else 0 for i in range(len(samples))])
    return medians
#end func