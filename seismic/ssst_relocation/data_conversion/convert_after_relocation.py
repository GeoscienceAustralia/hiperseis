import argparse, os
import numpy as np

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
    Determine if a wave arrives in P phase or S phase, looking at the last lag.
    
    
    Parameters
    ----------
    phase : string
        Phase
        
    
    Returns
    -------
    isp : integer
        1 if P wave, -1 if S wave, 0 if undetermined.
        
        
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
    
def convert_pick_array_to_list(pick_array, output_path, include_tcor=False):    
    """
    Converts an array of picks in the format used by Source Specific Station 
    Term relocation method into the required output format.
    
    
    Parameters
    ----------
    pick_array : numpy.ndarray
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
    
    output_path : string
        Output directory.
        
    include_tcor : boolean
        If a column is required for time corrections, set include_tcor = True.
        
    
    Returns
    -------
    p_combined : list
        List of P wave picks with format as below.
        pick_rows = [event_id, origin_time, mag, elon, elat, edepth (km), net,
                     sta, cha, arrival_time, slon, slat, az, baz, dist, 
                     residual, snr, qualityMeasureCWT, domFreq, 
                     qualityMeasureSlope, bandIndex, nsigma]
        
    s_combined : list
        List of S wave picks with format as for p_combined.
        
    all_picks : list
        Lisf of all picks with format as for p_combined and s_combined.
        
        
    """
    p_combined = list()
    s_combined = list()
    all_picks = list()
    
    ev_id = pick_array['event_id']
    origin_time = pick_array['origin_time']
    mag = pick_array['mag']
    elon = pick_array['elon']
    elat = 90.0 - pick_array['ecolat']
    edepthkm = pick_array['edepth']/1e3
    net = pick_array['net']
    sta = pick_array['stat']
    cha = pick_array['cha']
    arrival_time = pick_array['arrival_time']
    slon = pick_array['slon']
    slat = 90.0 - pick_array['scolat']
    az = azimuth(pick_array['elon'], pick_array['ecolat'], pick_array['slon'], 
                 pick_array['scolat'])
    baz = azimuth(pick_array['slon'], pick_array['scolat'], pick_array['elon'], 
                  pick_array['ecolat'])
    dist = ang_dist(pick_array['elon'], pick_array['ecolat'], 
                    pick_array['slon'], pick_array['scolat'])
    residual = pick_array['residual']
    snr = pick_array['snr']
    qualityMeasureCWT = pick_array['qualityMeasureCWT']
    domFreq = pick_array['domFreq']
    qualityMeasureSlope = pick_array['qualityMeasureSlope']
    bandIndex = pick_array['bandIndex']
    nSigma = pick_array['nSigma']
    
    phase = pick_array['phase']
    tcor = pick_array['tcor']
    
    cha[cha == ''] = 'None'
    
    with open(os.path.join(output_path, 'out.txt'), 'w') as file:
        file.write('Converting pick array \n')
    #end with
    
    for i in range(len(pick_array)):
        if (i+1) % 10000 == 0:
            with open(os.path.join(output_path, 'out.txt'), 'a') as file:
                file.write(str('Converting pick ' + str(i+1) + ' of ' + \
                               str(len(pick_array)) + '\n'))
            #end with
        #end if
        row = [ev_id[i], origin_time[i], mag[i], elon[i], elat[i], edepthkm[i], 
               net[i], sta[i], cha[i], arrival_time[i], phase[i], slon[i], 
               slat[i], az[i], baz[i], dist[i], residual[i], snr[i], 
               qualityMeasureCWT[i], domFreq[i], qualityMeasureSlope[i], 
               bandIndex[i], nSigma[i]]
        if include_tcor:
            row.append(tcor[i])
        #end if
        all_picks.append(row)
        if phase[i] in ['P', 'Pg']:
            p_combined.append(row)
        elif phase[i] in ['S', 'Sg']:
            s_combined.append(row)
        #end if
    #end for
    return p_combined, s_combined, all_picks
#end func
    
def write_to_csv(lst, filename, include_tcor=False):
    """
    Write a list of picks to file.
    
    
    Parameters
    ----------
    lst : list
        List of picks with row format described below in 'header'.
        
    filename : string
        File name.
        
    include_tcor : boolean
        If a column is required for time corrections, set include_tcor = True.        
        
        
    """
    header = str('#eventID originTimestamp mag originLon originLat ' + \
                 'originDepthKm net sta cha pickTimestamp phase ' + \
                 'stationLon stationLat az baz distance ttResidual snr ' + \
                 'qualityMeasureCWT domFreq qualityMeasureSlope bandIndex ' + \
                 'nSigma')
    if include_tcor:
        header = header + ' tcor'
    #end if
    with open(filename, 'w') as file:
        file.write(str(header + '\n'))
        for row in lst:
            string = str(' '.join([str(item) for item in row]) + '\n')
            file.write(string)
        #end for
    #end with
#end func            

def process():
    """
    Read in pick array from output file of Source Specific Station Term
    relocation method, convert into pick list, and output picks as .txt files.
    
    Input format - Structured array with the following data type.
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
    
    Output format - List with rows as follows.
        pick_rows = [event_id, origin_time, mag, elon, elat, edepth (km), net,
                     sta, cha, arrival_time, slon, slat, az, baz, dist, 
                     residual, snr, qualityMeasureCWT, domFreq, 
                     qualityMeasureSlope, bandIndex, nsigma]
    
    
    Arguments
    ---------
    event_file : string
        Name of numpy binary file to extrack picks from.
        
    output_path : string
        Output path.
        
    include_tcor : boolean
        If a column is required for time corrections, set include_tcor = True.
        
    
    Usage
    -----
    python convert_after_relocation.py --event_file .../event_file.npy
    --output_path .../output_path/
    
    
    """
    parser = argparse.ArgumentParser(description='Catalogue converter')

    parser.add_argument("--event_file", type=str, required=True)
    parser.add_argument("--output_path", type=str, default=".")
    parser.add_argument("--include_tcor", type=bool, default=False)
    args = parser.parse_args()
    
    event_file = args.event_file
    output_path = args.output_path
    include_tcor = args.include_tcor
    
    pick_array = np.load(event_file)
    
    p_combined, s_combined, all_picks = \
        convert_pick_array_to_list(pick_array, output_path, 
                                   include_tcor=include_tcor)
    
    filename = os.path.join(output_path, 'ensemble.p.txt')
    write_to_csv(p_combined, filename, include_tcor=include_tcor)
    
    filename = os.path.join(output_path, 'ensemble.s.txt')
    write_to_csv(s_combined, filename, include_tcor=include_tcor)
    
    filename = os.path.join(output_path, 'ensemble.all.txt')
    write_to_csv(all_picks, filename, include_tcor=include_tcor)
#end func
    
if __name__ == '__main__':
    process()
#end if