# -*- coding: utf-8 -*-
"""
Created on Thu Dec 23 11:42:59 2021

@author: U37509
"""

import argparse
import numpy as np
import os

def get_fields(arr1, fields):
    dtype = [(field, str(arr1[field].dtype)) for field in fields]
    arr2 = np.empty(len(arr1), dtype=dtype)
    for field in fields:
        arr2[field] = arr1[field]
    #end for
    return arr2
#end func

def process():
    parser = argparse.ArgumentParser(description='Results')
    
    parser.add_argument("--input_files", type=str, nargs='+', required=True)
    parser.add_argument("--output_path", type=str, default='.')
    
    """
    import sys
    sys.argv = ['/g/data/ha3/la8536/SSST/relocation/results.py', 
                "--input_file", 
                '/g/data/ha3/la8536/SSST/input_events/picks.npy',
                "--output_path", 
                '/g/data/ha3/la8536/SSST/input_events/']
    """
    
    args = parser.parse_args()
    
    input_files = args.input_files
    output_path = args.output_path
    
    for file in input_files:
        picks = np.load(file)
        
        hypocentres = get_fields(picks, ['event_id', 'elon', 'ecolat', 
                                         'edepth', 'origin_time'])
        _, inds = np.unique(hypocentres['event_id'], return_index = True)
        hypocentres = hypocentres[inds]
        
        residuals = get_fields(picks, ['phase', 'ptt', 'tcor', 'residual'])
        
        prefix1 = str(os.path.basename(file).split('.')[0] + '_')
        
        np.save(os.path.join(output_path, str(prefix1 + 'hypocentres')),
                hypocentres)
        np.save(os.path.join(output_path, str(prefix1 + 'residuals')),
                residuals)
    #end for
#end func

if __name__ == '__main__':
    process()
#end if