# -*- coding: utf-8 -*-
"""
Created on Thu Dec 23 16:08:36 2021

@author: U37509
"""

import numpy as np
from obspy.taup import TauPyModel
import argparse, configparser, csv, glob, os

def gradient(x, y, f):
    g = np.ones_like(f)*np.nan
    gx = np.ones_like(f)*np.nan
    gy = np.ones_like(f)*np.nan
    
    for i in range(len(x)):
        for j in range(len(y)):
            if np.isnan(f[i,j]):
                continue
            #end if
            if i == 0:
                gx[i,j] = (f[i+1,j] - f[i,j])/(x[i+1] - x[i])
            elif i == len(x)-1:
                gx[i,j] = (f[i,j] - f[i-1,j])/(x[i] - x[i-1])
            else:
                if np.isnan(f[i+1,j]):
                    gx[i,j] = (f[i,j] - f[i-1,j])/(x[i] - x[i-1])
                elif np.isnan(f[i-1,j]):
                    gx[i,j] = (f[i+1,j] - f[i,j])/(x[i+1] - x[i])
                else:
                    gx[i,j] = (f[i+1,j] - f[i-1,j])/(x[i+1] - x[i-1])
                #end if
            #end if
            if j == 0:
                gy[i,j] = (f[i,j+1] - f[i,j])/(y[j+1] - y[j])
            elif j == len(y)-1:
                gy[i,j] = (f[i,j] - f[i,j-1])/(y[j] - y[j-1])
            else:
                if np.isnan(f[i,j+1]):
                    gy[i,j] = (f[i,j] - f[i,j-1])/(y[j] - y[j-1])
                elif np.isnan(f[i, j-1]):
                    gy[i,j] = (f[i,j+1] - f[i,j])/(y[j+1] - y[j])
                else:
                    gy[i,j] = (f[i,j+1] - f[i,j-1])/(y[j+1] - y[j-1])
                #end if
            #end if
        #end for
    #end for
    g = np.sqrt(gx**2 + gy**2)
    return g, gx, gy
#end func
    
def read_iloc_TT_table(filename, nanval=-999.0):
    """
    Reads a travel time table for a given phase, in the format of those 
    provided with the distribution of the iLoc earthquake location algorithm.
    
    The structure of these files are as follows.
    - Header - information about the file.
    - Delta samples - a few rows will contain the epicentral distance values
        used for the grid.
    - Depth samples - a single row contains depth samples used.
    - Travel times - a grid of travel times, with one row per epicentral 
        distance sample, and one row per depth sample.
    - dt/dd - a grid containing first horizontal derivative of the travel time
        table (\partial time/\partial distance).
    - dt/dh - a grid containing the first vertical derivative of the travel 
        time table (\partial time/\partial height).
    
    
    Parameters
    ----------
    filename : string
        Name of travel time table.
        
        
    Returns
    -------
    ecdists : numpy.ndarray
        1D array of epicentral distance grid points used for the table.
        
    depths : numpy.ndarray
        1D array of depth grid points used for the table.
        
    times : numpy.ndarray
        2D array of travel times with dimensions (distance, depth).
        
    dtdd_values : numpy.ndarray
        2D array containing values of the partial derivative of 'times'
        with respect to epicentral distance.
    
    
    """
    ecdists = []
    depths = []
    times = []
    dtdd_values = []
    dists_start = 0
    dists_end = 0
    depths_start = 0
    depths_end = 0
    times_start = 0
    times_end = 0
    dtdd_values_start = 0
    dtdd_values_end = 0
    with open('%s'%filename) as file:
        reader = csv.reader(file, delimiter=' ')
        ind = 0
        for row in reader:
            if ' '.join(row) == '# delta samples': 
                dists_start = ind + 1
            elif ' '.join(row) == '# depth samples': 
                dists_end = ind - 1
                depths_start = ind + 1
            elif ' '.join(row) == \
                '# travel times (rows - delta, columns - depth)':
                depths_end = ind - 2
                times_start = ind + 2
            elif ' '.join(row) == '# dtdd (rows - delta, columns - depth)':
                times_end = ind - 2
                dtdd_values_start = ind + 2
            elif ' '.join(row) == '# dtdh (rows - delta, columns - depth)':
                dtdd_values_end = ind - 2
            ind += 1
        #end for
    with open('%s'%filename) as file:
        reader = csv.reader(file, delimiter=' ', skipinitialspace=True)
        ind = 0
        for row in reader:
            if ind in range(dists_start, dists_end+1):
                ecdists = ecdists + [float(item) for item in row if item != '']
            elif ind in range(depths_start, depths_end+1):
                depths = depths + [float(item) for item in row if item != '']
            elif ind in range(times_start, times_end+1):
                times = times + [[float(item) for item in row if item != '']]
            elif ind in range(dtdd_values_start, dtdd_values_end+1):
                dtdd_values = dtdd_values + \
                              [[float(item) for item in row if item != '']]
            #end if
            ind += 1
        #end for
    ecdists = np.array(ecdists)
    depths = np.array(depths)
    dtdd_values = np.array(dtdd_values)
    times = np.array(times)
    times[times == nanval] = np.nan
    dtdd_values[dtdd_values == nanval] = np.nan
    return ecdists, depths, times, dtdd_values
#end func
    
def write_to_csv(ecdists, depths, tt, dtdd, dtdh, phase, output_path):
    tt[np.isnan(tt)] = -999.0
    dtdd[np.isnan(dtdd)] = -999.0
    dtdh[np.isnan(dtdh)] = -999.0
    with open(os.path.join(output_path, str(phase + '.tt')), 'w') as file:
        file.write(' '.join([("%4.6f"%val).rjust(11) for val in ecdists]) + \
                   '\n')
        file.write(' '.join([("%4.6f"%val).rjust(11) for val in depths]) + \
                   '\n')
        for row in tt:
            file.write(' '.join([("%4.6f"%val).rjust(11) for val in row]) + \
                       '\n')
        #end for
    #end with
    with open(os.path.join(output_path, str(phase + '.dtdd')), 'w') as file:
        file.write(' '.join([("%4.6f"%val).rjust(11) for val in ecdists]) + \
                   '\n')
        file.write(' '.join([("%4.6f"%val).rjust(11) for val in depths]) + \
                   '\n')
        for row in dtdd:
            file.write(' '.join([("%4.6f"%val).rjust(11) for val in row]) + \
                       '\n')
        #end for
    #end with
    with open(os.path.join(output_path, str(phase + '.dtdh')), 'w') as file:
        file.write(' '.join([("%4.6f"%val).rjust(11) for val in ecdists]) + \
                   '\n')
        file.write(' '.join([("%4.6f"%val).rjust(11) for val in depths]) + \
                   '\n')
        for row in dtdh:
            file.write(' '.join([("%4.6f"%val).rjust(11) for val in row]) + \
                       '\n')
        #end for
    #end with
#end func
    
def compute_using_taup(ecdists, depths, phase, model):
    times = np.ones((len(ecdists), len(depths)))*np.nan
    
    for i in range(len(ecdists)):
        for j in range(len(depths)):
            arrivals = model.get_travel_times(source_depth_in_km=depths[j],
                                              distance_in_degree=ecdists[i],
                                              phase_list=phase,
                                              receiver_depth_in_km=0.0)
            if arrivals != list():
                times[i, j] = arrivals[0].time
            #end if
        #end for
    #end for
    
    return times
#end func
    
def process():
    parser = argparse.ArgumentParser(description='Travel time table generator')
    
    parser.add_argument("--iloc_path", type=str, default='')
    parser.add_argument("--output_path", type=str, default='.')
    parser.add_argument("--from_iloc_tables", type=bool, default=False)
    parser.add_argument("--config_file", type=str, default='')
    parser.add_argument("--model", type=str, default='iasp91')
    
    """
    import sys
    sys.argv = ['/g/data/ha3/la8536/SSST/data_conversion/tt_table_calculator.py',
                "--iloc_path", '/g/data/ha3/la8536/SSST/TT_tables/',
                "--output_path", '/g/data/ha3/la8536/SSST/TT_tables/',
                "--from_iloc_tables", 'True']
    """
    
    args = parser.parse_args()
    output_path = args.output_path
    from_iloc_tables = args.from_iloc_tables
    
    if from_iloc_tables:
        print('Reading travel times from iLoc tables')
        iloc_path = args.iloc_path
        for file in glob.glob(os.path.join(iloc_path, '*.tab')):
            phase = os.path.basename(file).split('.')[0]
            ecdist, depth, tt, _ = \
                read_iloc_TT_table(file)
            _, dtdd, dtdh = gradient(ecdist, depth, tt)
            write_to_csv(ecdist, depth, tt, dtdd, dtdh, phase, output_path)
        #end for
    else:
        print('Computing travel times using taup')
        
        model = TauPyModel(model=args.model)
        
        config = configparser.ConfigParser()
        config.sections()
        config.read(args.config_file)
        
        sections = config._sections.keys()
        
        phases = list()
        ecdists = list()
        depths = list()
        
        for section in sections:
            phases.append(section.split())
            ecdists.append(np.array(config[section]['dist'].split()) \
                           .astype(float))
            depths.append(np.array(config[section]['depth'].split()) \
                          .astype(float))
        #end for
        
        for i in range(len(phases)):
            print('Finding travel times for', phases[i])
            phase = phases[i]
            ecdist = ecdists[i]
            depth = depths[i]
            tt = compute_using_taup(ecdist, depth, phase, model)
            _, dtdd, dtdh = gradient(ecdist, depth, tt)
            write_to_csv(ecdist, depth, tt, dtdd, dtdh, phase[0], output_path)
        #end for
    #end if
#end func

if __name__ == '__main__':
    process()
#end if