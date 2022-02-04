"""
Description
-----------
This module is used by the event relocation and phase redefinition algorithm 
for travel time prediction.

Developer: Lachlan Adams 
Contact: lachlan.adams@ga.gov.au or lachlan.adams.1996@outlook.com

"""

import glob, os
import numpy as np
from scipy import interpolate

class tt_table_object():
    def __init__(self, phase, ecdists, depths, tt, dtdd, dtdh, nanval=-999.0):
        """
        Object designed to allow interpolation of travel time tables for a 
        particular phase for prediction of travel times for seismic waves.
        The object will contain epicentral distance an depth value arrays, the
        travel time ('tt') table and its first vertical (dtdh) and horizontal 
        (dtdd) derivatives, and a 2D piecewise cubic bezier polynomial 
        interpolating function for each of the 'tt', 'dtdd', and 'dtdh' tables.
        
        
        Parameters
        ----------
        phase : string
        
        ecdists : numpy.ndarray
            Array of epicentral distance samples for tables.
            
        depths : numpy.ndarray
            Array of depth samples for tables.
            
        tt : numpy.ndarray
            Array of predicted travel time values.
            
        dtdd : numpy.ndarray
            Horizontal derivative of 'tt' table.
            
        dtdh : numpy.ndarray
            Vertical derivative of 'dtdh' table.
            
        nanval : float
            Value in tables representing NaN (not a number).
            
            
        """
        
        
        self.ecdists = ecdists
        self.depths = depths
        self.tt = tt
        self.dtdh = dtdh
        self.dtdd = dtdd
        self.name = phase
        self.max_dist = np.max(ecdists)
        self.max_depth = np.max(depths)
        
        self.interp_tt(nanval)
        self.interp_dtdd(nanval)
        self.interp_dtdh(nanval)
    #end func
    
    def interp_tt(self, nanval):
        """
        Interpolate the travel time table.
        
        
        """
        values = np.reshape(np.transpose(self.tt), self.tt.size)
        x, y = np.meshgrid(self.ecdists, self.depths)
        X = np.reshape(x, len(values))
        Y = np.reshape(y, len(values))
        points = np.transpose(np.array([X, Y]))
        points = points[values != nanval]
        values = values[values != nanval]
        self.tt_interp = interpolate.CloughTocher2DInterpolator(points, values)
        self.tt = np.transpose(self.tt_interp(x, y))
    #end func
    
    def interp_dtdh(self, nanval):
        """
        Interpolate the first vertical derivative of the travel time table.
        
        
        """
        values = np.reshape(np.transpose(self.dtdh), self.dtdh.size)
        x, y = np.meshgrid(self.ecdists, self.depths)
        X = np.reshape(x, len(values))
        Y = np.reshape(y, len(values))
        points = np.transpose(np.array([X, Y]))
        points = points[values != nanval]
        values = values[values != nanval]
        self.dtdh_interp = interpolate.CloughTocher2DInterpolator(points, 
                                                                  values)
        self.dtdh = np.transpose(self.dtdh_interp(x, y))
    #end func
    
    def interp_dtdd(self, nanval):
        """
        Interpolate the first horizontal derivative of the travel time table.
        
        
        """
        values = np.reshape(np.transpose(self.dtdd), self.dtdd.size)
        x, y = np.meshgrid(self.ecdists, self.depths)
        X = np.reshape(x, len(values))
        Y = np.reshape(y, len(values))
        points = np.transpose(np.array([X, Y]))
        points = points[values != nanval]
        values = values[values != nanval]
        self.dtdd_interp = interpolate.CloughTocher2DInterpolator(points, 
                                                                  values)
        self.dtdd = np.transpose(self.dtdd_interp(x, y))
    #end func
    
class ellipcorr_object():
    def __init__(self, ecdist, depth, tau0, tau1, tau2):
        """
        Object used to interpolate and store coefficient tables for computation 
        of ellipticity corrections. If more than 3 epicentral distance samples
        are available, a rectilinear bivariate cubic spline is used to 
        interpolate the table, or if not, bilinear interpolation is used.
        
        
        Parameters
        ----------
        ecdist : numpy.ndarray
            Epicentral distance samples for coefficient tables.
            
        depth : numpy.ndarray
            Depth samples for coefficient tables.
            
        tau0 : numpy.ndarray
            tau0 value table.
            
        tau1 : numpy.ndarray
            tau0 value table.
            
        tau2 : numpy.ndarray
            tau0 value table.
            
            
        """
        self.ecdist = ecdist
        self.depth = depth
        self.tau0 = tau0
        self.tau1 = tau1
        self.tau2 = tau2
        if len(ecdist) < 4:
            self.method = 'linear'
            self.tau0_interp = \
                interpolate.RegularGridInterpolator((ecdist, depth), tau0, 
                                                    method='linear',
                                                    bounds_error=False, 
                                                    fill_value=0)
            self.tau1_interp = \
                interpolate.RegularGridInterpolator((ecdist, depth), tau1, 
                                                    method='linear',
                                                    bounds_error=False, 
                                                    fill_value=0)
            self.tau2_interp = \
                interpolate.RegularGridInterpolator((ecdist, depth), tau2, 
                                                    method='linear',
                                                    bounds_error=False, 
                                                    fill_value=0)
        else:
            self.method = 'cubic'
            self.tau0_interp = \
                interpolate.RectBivariateSpline(ecdist, depth, tau0)
            self.tau1_interp = \
                interpolate.RectBivariateSpline(ecdist, depth, tau1)
            self.tau2_interp = \
                interpolate.RectBivariateSpline(ecdist, depth, tau2)
    #end func
    
    def compute_values(self, ecdist, depth):
        """
        Function to compute ellipticity correction coefficients using 
        interpolating functions. If given coordinates are outside coordinate 
        range, values of 0.0 are returned.
        
        
        Parameters
        ----------
        ecdist : float, or numpy.ndarray
            Epicentral distance values.
            
        depth : float, or numpy.ndarray
            Depth values. Must have same length as 'ecdist'.
            
        
        Returns
        -------
        tau0 : numpy.ndarray
        
        tau1 : numpy.ndarray
        
        tau2 : numpy.ndarray
        
        
        """
        if self.method == 'linear':
            tau0 = self.tau0_interp((ecdist, depth))
            tau1 = self.tau1_interp((ecdist, depth))
            tau2 = self.tau2_interp((ecdist, depth))
        elif self.method == 'cubic':
            in_bbox1 = np.logical_and(depth <= np.max(self.depth),
                                      depth >= np.min(self.depth))
            in_bbox2 = np.logical_and(ecdist <= np.max(self.ecdist),
                                      ecdist >= np.min(self.ecdist))
            in_bbox = np.logical_and(in_bbox1, in_bbox2)
            tau0 = self.tau0_interp(ecdist, depth, grid=False)
            tau1 = self.tau1_interp(ecdist, depth, grid=False)
            tau2 = self.tau2_interp(ecdist, depth, grid=False)
            tau0[~in_bbox] = 0.0
            tau1[~in_bbox] = 0.0
            tau2[~in_bbox] = 0.0
        #end if
        return tau0, tau1, tau2
    #end func
#end class
    
def read_tt_table(input_path, phase):
    """
    Reads pre-computed travel time tables. For each file, the first row 
    contains coordinates of epicentral distance samples, the second contains 
    depth samples, and each subsequent row contains travel time values. There 
    is one row of values in the table for each epicentral distance sample.
    
    
    Parameters
    ----------
    input_path : string
        Directory containing .tt, .dtdd, and .dtdh files for the required 
        phase.
        
    phase : string
    
    
    Returns
    -------
    ecdists : numpy.ndarray
        Array of epicentral distance values.
        
    depths : numpy.ndarray
        Array of depth values.
        
    tt : numpy.ndarray
        Table of predicted travel time values.
        
    dtdd : numpy.ndarray
        Table of horizontal derivatives of the predicted travel time values.
        
    dtdh : numpy.ndarray
        Table of vertical derivatives of the predicted travel time values.
        
    
    """
    
    
    ttfile = os.path.join(input_path, str(phase + '.tt'))
    dtddfile = os.path.join(input_path, str(phase + '.dtdd'))
    dtdhfile = os.path.join(input_path, str(phase + '.dtdh'))
    
    with open(ttfile, 'r') as file:
        rows = file.read().splitlines()
        ecdists = np.array(rows[0].split()).astype(float)
        depths = np.array(rows[1].split()).astype(float)
        tt = np.array([row.split() for row in rows[2:]]).astype(float)
    #end with
    with open(dtddfile, 'r') as file:
        rows = file.read().splitlines()
        dtdd = np.array([row.split() for row in rows[2:]]).astype(float)
    #end with
    with open(dtdhfile, 'r') as file:
        rows = file.read().splitlines()
        dtdh = np.array([row.split() for row in rows[2:]]).astype(float)
    #end with
    
    return ecdists, depths, tt, dtdd, dtdh
#end func
    
def read_ellipcorr_table(filename, depth=np.array([0.0, 100.0, 200.0, 300.0, 
                                                   500.0, 700.0])):
    """
    Reads ellipticity correction coefficients from file in the format used by 
    the ellipcorr.f Fortran subroutine written by Brian Kennett.
    File contains a single string. When partitioned into strings of length 80
    the format is as follows.
    
    [[phase, number of distance samples, minimum distance, maximum distance],
     [distance value],
     [tau0 values],
     [tau1 values],
     [tau2 values]]
    
    Extra rows are added for each distance value. There are six items in each 
    tau0, tau1, and tau2 row which correspond to six hypocentre depth values 
    (0, 100, 200, 300, 500, 700).
    
    
    Parameters
    ----------
    filename : string
    
    depth : numpy.ndarray
        The depth values used to create the coefficient tables. These values 
        are fixed.
        
    
    Returns
    -------
    dct : dictionary
        Dictionary containing an "ellipcorr_object" for each available phase.
        This object contains tables of tau0, tau1, and tau2 values, and 
        epicentral distance and depth coordinates for these values, as well as 
        interpolating functions of the tables.
        
    phase_list : list of string
        List of phases for which ellipticity corrections are available.
        
        
    """
    
    
    with open(filename, 'r') as file:
        string = file.read()
    #end with

    l = len(string)
    ind = np.arange(0, l+80, 80)
    rows = list()
    phase_list = list()

    for i in range(len(ind)-1):
        row = string[ind[i]:ind[i+1]].split()
        rows.append(row)
        if len(row) == 4:
            phase_list.append(row[0])
        #end if
    #end for

    dct = {}
    first = True
    phase = None
    ecdist = list()
    tau = 0
    tau0 = list()
    tau1 = list()
    tau2 = list()
    for row in rows:
        if len(row) == 4:
            if not first:
                ecdist = np.array(ecdist).astype(float)
                tau0 = np.array(tau0).astype(float)
                tau1 = np.array(tau1).astype(float)
                tau2 = np.array(tau2).astype(float)
                dct[phase] = ellipcorr_object(ecdist, depth, tau0, tau1, tau2)
                ecdist = list()
                tau0 = list()
                tau1 = list()
                tau2 = list()
            else:
                first = False
            #end if
            phase = row[0]
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
    
    return dct, phase_list
#end func

def process_tt_tables(tt_table_path):
    """
    Processes a list of travel time tables so that the tables are interpolated
    prior to obtaining predicted travel times from them. The 'NaN' value used
    in the tables is assumed to be -999.0.
    
    
    Parameters
    ----------
    tt_table_path : string
        Directory in which travel time tables ('.tt') are stored.
        
    
    Returns
    -------
    dct : dictionary
        Dictionary of (phase: tt_table_object) pairs. Phases are the 
        names of the travel time tables.
        
    phase_list : list
        List of strings representing names of phases that travel time tables 
        are available for.
        
        
    """
    filelist = [os.path.basename(file) for file in \
                glob.glob('%s/*.tt'%tt_table_path)]
    phase_list = np.array([file.split('.')[0] for file in filelist])
    dct = {}
    for phase in phase_list:
        ecdists, depths, tt, dtdd, dtdh = read_tt_table(tt_table_path, phase)
        dct[phase] = tt_table_object(phase, ecdists, depths, tt, dtdd, dtdh, 
                                     nanval=-999.0)
    #end for
    return dct, phase_list
#end func
    
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

def azimuth(lon1, colat1, lon2, colat2, units='degrees'):
    """
    Function to calculate the azimuth from (lon1, lat1) to (lon2, lat2).
    
    
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
    azim : float
        Azimuth from (lon1, colat1) to (lon2, colat2).
        
        
    """
    if units == 'degrees':
        degrad = np.pi/180.0
        a = lon1*degrad
        b = colat1*degrad
        x = lon2*degrad
        y = colat2*degrad
    else:
        a = lon1
        b = colat1
        x = lon2
        y = colat2
    #end if
    azim = np.arctan(np.sin(x - a)/(np.sin(b)/np.tan(y) - \
                                    np.cos(b)*np.cos(x - a)))
    
    offset1 = np.pi*(colat1 < colat2)
    offset2 = np.pi*(np.logical_and(colat1 == colat2, b > np.pi/2.0))
    
    azim = (azim + offset1 + offset2 + 2*np.pi) % (2*np.pi)
    if units == 'degrees': 
        azim = azim/degrad
    return azim
#end func
    
def IsP(phase):
    """
    Returns 1 if a wave is a P wave, -1 if S wave, or 0 if undetermined.
    
    
    Parameters
    ----------
    phase : string
    
    
    Returns
    -------
    isp : integer
    
    
    """
    isp = 0
    if phase == 'p':
        isp = 1
    elif phase == 's':
        isp = -1
    else:
        l = len(phase)
        for i in np.arange(l, 0, -1):
            char = phase[i-1]
            if char == 'P': 
                isp = 1
                break
            elif char == 'S' or char == 'L':
                isp = -1
                break
            #end if
        #end for
    #end if
    return isp
#end func
    
def elev_corr(wt, phase, ecdist, edepth, selev, dct):
    """
    Computes an elevation correction based on the surface wave velocity for 
    P or S waves. 
    
    
    Parameters
    ----------
    wt : integer, or numpy.ndarray
        Type of wave. 1 if P, -1 if S.
        
    phase : string, or list
        Name of phase.
        
    ecdist : float, or numpy.ndarray
        Angular distance from source to receiver.
        
    edepth : float, or numpy.ndarray
        Event depth (km).
        
    selev : float, or numpy.ndarray
        Station elevation (m).
    
    dct : dictionary
        Dictionary of phase_TT_table objects containing travel time tables.
        
    
    Returns
    -------
    elev_corr : float, or numpy.ndarray
        Time correction (s).
        
        
    """
    deg2km = 6371*np.pi/180
    PSurfVel = 5.8
    SSurfVel = 3.46
    if wt == 1:
        surfvel = PSurfVel
    elif wt == -1:
        surfvel = SSurfVel
    else:
        surfvel = 0.0
    #end if
    if surfvel != 0.0:
        if type(phase) == str:
            dt_dd = dct[phase].dtdd_interp(ecdist, edepth)
        elif type(phase) == list:
            dt_dd = np.array([dct[ph].dtdd_interp(ecdist, edepth) \
                              for ph in phase])
        #end if
        elev_corr = (surfvel*dt_dd/deg2km)**2
        if type(elev_corr) == float:
            if elev_corr > 1: elev_corr = 1.0/elev_corr
        else:
            elev_corr[elev_corr > 1] = 1.0/elev_corr[elev_corr > 1]
        #end if
        elev_corr = np.sqrt(1.0 - elev_corr)
        elev_corr = elev_corr*selev/(1000*surfvel)
    else:
        if type(ecdist) == float:
            elev_corr = 0.0
        else:
            elev_corr = np.zeros(len(ecdist))
        #end if
    #end if
    return elev_corr
#end func

def ellip_corr(phases, azim, ecdist, ecolat, depth, dct):
    degrad = np.pi/180
    
    azim = azim*degrad
    ecolat = ecolat*degrad
    
    s3 = np.sqrt(3.0)/2.0
    sc0 = 0.25*(1.0+3.0*np.cos(2.0*ecolat))
    sc1 = s3*np.sin(2.0*ecolat)
    sc2 = s3*np.sin(ecolat)*np.sin(ecolat)
    
    tcor = list()
    if type(phases) == str:
        phases = [phases]
    #end if
    
    for phase in phases:
        if phase in dct.keys():
            tau0, tau1, tau2 = dct[phase].compute_values(ecdist, depth)
            tcor.append(sc0*tau0 + sc1*np.cos(azim)*tau1 + \
                        sc2*np.cos(2.0*azim)*tau2)
        else:
            phase_alias = {'Pg': 'Pup', 'Sg': 'Sup', 'pPg': 'pP', 'sPg': 'sP', 
                           'pSg': 'pS', 'sSg': 'sS', 'Pb': 'P', 'Sb': 'S', 
                           'pPb': 'pP', 'sPb': 'sP', 'pSb': 'pS', 'sSb': 'sS', 
                           'Pn': 'P', 'Sn': 'S', 'pPn': 'pP', 'sPn': 'sP', 
                           'pSn': 'pS', 'sSn': 'sS', 'SPn': 'SP', 'SPb': 'SP', 
                           'SPg': 'SP', 'SnP': 'SP', 'PSn': 'PS', 'PnPn': 'PP', 
                           'SnSn': 'SS', 'p': 'Pup', 's': 'Sup', 
                           'Pdif': 'Pdiff', 'Sdif': 'Sdiff'}
            if phase in phase_alias.keys():
                tau0, tau1, tau2 = \
                    dct[phase_alias[phase]].compute_values(ecdist, depth)
                tcor.append(sc0*tau0 + sc1*np.cos(azim)*tau1 + \
                            sc2*np.cos(2.0*azim)*tau2)
            else:
                tcor.append(np.zeros_like(ecdist))
            #end if
        #end if
    #end for
    
    tcor = np.vstack(tcor)
    return tcor
#end func    
    
def travel_time(ecdist, edepth, phase_list, dct):
    """
    Compute travel times for waves by interpolation of travel time tables.
    Uses pre-computed travel time table 2D cubic bezier polynomial interpolant
    for prediction.
    
    
    Parameters
    ----------
    ecdist : numpy.ndarray, or float
        Angular distance from source to receiver.
        
    edepth : numpy.ndarray, or float
        Source depth.
        
    phase_list : list, or string
        Phase of wave to compute travel time for.
        
    dct : dictionary
        Dictionary of phase_TT_table objects containing pre-computed 
        interpolant for the table.
        
    
    Returns
    -------
    ptt : numpy.ndarray, or float
        Interpolated predicted travel time.
        
        
    """
    if type(phase_list) == str:
        ptt = dct[phase_list].tt_interp(ecdist, edepth)
    else:
        ptt = np.array([dct[phase].tt_interp(ecdist, edepth) \
                        for phase in phase_list])
    #end if
    return ptt
#end func
    
def compute_travel_time(wt, ecdist, azim, ecolat, edepth, selev, ett, 
                        phase_list, TT_dict, ellipcorr_dict, thr):
    """
    Find the predicted travel times for a wave by checking all possible phases.
    Uses elevation corrections for all phases, and computes ellipticity 
    corrections for the phases with travel time within 'thr' seconds of
    empirical travel time.
    
    
    Parameters
    ----------
    wt : integer
        Wave type (P=1, S=-1, undefined=0)
    
    ecdist : numpy.ndarray
        Angular distance between source and receiver.
        
    azim : numpy.ndarray
        Azimuth from receiver to source.
        
    ecolat : numpy.ndarray
        Epicentral colatitude.
        
    edepth : numpy.ndarray
        Event depth (m).
        
    selev : numpy.ndarray
        Height of station above sea level (m).
        
    ett : numpy.ndarray
        Empirical travel time (seconds).
        
    phase_list : list
        List of strings, representing names of phases to check.
        
    TT_dict : dictionary
        Dictionary of phase_TT_table objects to use to interpolate travel 
        times.
        
    ellipcorr_dict : dictionary
        Dictionary of ellipcorr_object objects to use to interpolate 
        ellipticity correction table.
        
    thr : float
        Threshold time to classify candidate predicted travel times.
        
    
    Returns
    -------
    phases : numpy.ndarray
        Array representing most likely phase for input waves.
        
    bestptt : numpy.ndarray
        Most likely travel time for redefined phases.
        
    minresid : numpy.ndarray
        Travel time residual between empirical travel times and predicted 
        travel time for redefined phases.
        
        
    """
    ptt = travel_time(ecdist, edepth, phase_list, TT_dict)
    elev_corr_val = elev_corr(wt, phase_list, ecdist, edepth, selev, TT_dict)
    ellip_corr_val = ellip_corr(phase_list, azim, ecdist, ecolat, edepth, 
                                ellipcorr_dict)
    
    ptt = ptt + elev_corr_val + ellip_corr_val
    resid = ett - ptt
    
    minind = np.where(np.abs(resid) == np.nanmin(np.abs(resid), axis=0))
    minind = np.transpose(minind)[np.transpose(minind)[:, 1].argsort()]
    
    phases = np.array(['Px' if wt == 1 else 'Sx' if wt == -1 else 'X' \
                       for i in range(len(ett))], dtype='U8')
    minresid = np.zeros_like(ett)
    bestptt = np.zeros_like(ett)
    
    for i, j in minind:
        if np.abs(resid[i, j]) < thr:
            phases[j] = phase_list[i]
            minresid[j] = resid[i, j]
            bestptt[j] = ptt[i, j]
        #end if
    #end for
    
    return phases, bestptt, minresid
#end func
    
def predict_travel_times(picks, phase_list, TT_dict, ellipcorr_dict, config):
    """
    Predicts travel times and redefines phases for input picks.
    First pass - uses "use_phases" only.
    Second pass - if a match wasn't found from the above list, tries all 
    available phases.
    
    
    Parameters
    ----------
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
        
    phase_list : list
        List of phases for which travel time tables are available.
        
    TT_dict : dictionary
        Dictionary of phase_TT_table objects to use to interpolate travel 
        times.
        
    ellipcorr_dict : dictionary
        Dictionary of ellipcorr_object objects to use to interpolate 
        ellipticity correction table.
        
    config : configparser.SectionProxy object
        Information from config file.
        
    
    Returns
    picks : numpy.ndarray
        Input pick array with updated predicted travel time, phase, and 
        travel time residual values.
        
        
    """
    use_phases = config['phases'].split(', ')
    thr_p = float(config['thr_p'])
    thr_s = float(config['thr_s'])
    
    slon = picks['slon']
    scolat = picks['scolat']
    selev = picks['selev']
    elon = picks['elon']
    ecolat = picks['ecolat']
    edepth = picks['edepth']/1e3
    ett = picks['arrival_time'] - picks['origin_time']
    wave_types = np.array([IsP(pick['phase']) for pick in picks])
    ecdist = ang_dist(elon, ecolat, slon, scolat)
    azim = azimuth(elon, ecolat, slon, scolat)
    
    # First pass
    ind_p = wave_types == 1
    ind_s = wave_types == -1
    ind_x = wave_types == 0

    is_p = [IsP(phase) for phase in use_phases]
    p_phase_list = [use_phases[i] for i in range(len(use_phases)) \
                    if is_p[i] == 1]
    s_phase_list = [use_phases[i] for i in range(len(use_phases)) \
                    if is_p[i] == -1]
    
    phases_p, ptt_p, resid_p = \
        compute_travel_time(1, ecdist[ind_p], azim[ind_p], ecolat[ind_p], 
                            edepth[ind_p], selev[ind_p], ett[ind_p], 
                            p_phase_list, TT_dict, ellipcorr_dict, thr_p)
    
    phases_s, ptt_s, resid_s = \
        compute_travel_time(-1, ecdist[ind_s], azim[ind_s], ecolat[ind_s], 
                            edepth[ind_s], selev[ind_s], ett[ind_s], 
                            s_phase_list, TT_dict, ellipcorr_dict, thr_s)
    
    phases_x, ptt_x, resid_x = \
        compute_travel_time(0, ecdist[ind_x], azim[ind_x], ecolat[ind_x], 
                            edepth[ind_x], selev[ind_x], ett[ind_x], 
                            use_phases, TT_dict, ellipcorr_dict, thr_s)
        
    picks['phase'][ind_p] = phases_p
    picks['phase'][ind_s] = phases_s
    picks['phase'][ind_x] = phases_x
    
    picks['ptt'][ind_p] = ptt_p
    picks['ptt'][ind_s] = ptt_s
    picks['ptt'][ind_x] = ptt_x
    
    picks['residual'][ind_p] = resid_p
    picks['residual'][ind_s] = resid_s
    picks['residual'][ind_x] = resid_x
    
    # Second pass
    ind_p = picks['phase'] == 'Px'
    ind_s = picks['phase'] == 'Sx'
    ind_x = picks['phase'] == 'X'
    
    is_p = [IsP(phase) for phase in phase_list]
    p_phase_list = [phase_list[i] for i in range(len(phase_list)) \
                    if is_p[i] == 1 and phase_list[i] not in use_phases]
    s_phase_list = [phase_list[i] for i in range(len(phase_list)) \
                    if is_p[i] == -1 and phase_list[i] not in use_phases]
    
    phases_p, ptt_p, resid_p = \
        compute_travel_time(1, ecdist[ind_p], azim[ind_p], ecolat[ind_p], 
                            edepth[ind_p], selev[ind_p], ett[ind_p], 
                            p_phase_list, TT_dict, ellipcorr_dict, thr_p)
    
    phases_s, ptt_s, resid_s = \
        compute_travel_time(-1, ecdist[ind_s], azim[ind_s], ecolat[ind_s], 
                            edepth[ind_s], selev[ind_s], ett[ind_s], 
                            s_phase_list, TT_dict, ellipcorr_dict, thr_s)
        
    phases_x, ptt_x, resid_x = \
        compute_travel_time(0, ecdist[ind_x], azim[ind_x], ecolat[ind_x], 
                            edepth[ind_x], selev[ind_x], ett[ind_x], 
                            phase_list, TT_dict, ellipcorr_dict, thr_s)
        
    picks['phase'][ind_p] = phases_p
    picks['phase'][ind_s] = phases_s
    picks['phase'][ind_x] = phases_x
    
    picks['ptt'][ind_p] = ptt_p
    picks['ptt'][ind_s] = ptt_s
    picks['ptt'][ind_x] = ptt_x
    
    picks['residual'][ind_p] = resid_p
    picks['residual'][ind_s] = resid_s
    picks['residual'][ind_x] = resid_x
    
    return picks
#end func