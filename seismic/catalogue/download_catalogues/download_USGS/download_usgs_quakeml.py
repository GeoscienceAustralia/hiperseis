"""
Description:
    This script will download a catalogue of events from the ISC database given
    various input criteria by the user. An output file will be generated for 
    each month or each day in the requested time span (depending on argument 
    "--split_by_day".
    
Developer: Lachlan Adams 
Contact: lachlan.adams@ga.gov.au or lachlan.adams.1996@outlook.com

"""

import argparse, glob, os
from ftplib import FTP
import numpy as np

def process():
    """
    Download USGS catalogue fulfilling specified criteria.
    
    
    Usage
    -----
    python download_usgs_quakeml.py --start_date YY-MM-DD --end_date YYYY-MM-DD 
    --output_path .../output_path
    
    
    """
    parser = argparse.ArgumentParser(description='ISC event downloader')

    parser.add_argument("--start_date", type=str, default=None,
                        help="Enter as '--start_date YYYY-MM-DD'")
    parser.add_argument("--end_date", type=str, default=None,
                        help="Enter as '--end_date YYYY-MM-DD'")
    parser.add_argument("--output_path", type=str, default='.')
    
    args = parser.parse_args()
    output_path = args.output_path
    start_date = args.start_date
    end_date = args.end_date
    
    start_year, start_month, start_day = \
        (int(item) for item in start_date.split('-'))
    end_year, end_month, end_day = \
        (int(item) for item in end_date.split('-'))
    
    download_events(start_year, start_month, start_day, end_year, end_month, 
                    end_day, output_path)    
#end func
    
def days_in_month(year, month):
    """
    Returns the number of days in a month, taking into account leap years.
    January, March, May, July, August, October, December = 31 Days.
    April, June, September, November = 30 Days.
    February:
        if year is not a multiple of 4 = 28 days
        if year is a multiple of 4:
            if year is a multiple of 400: 29 days
            if year is a multiple of 100 but not of 400: 28 days
            otherwise 29 days.
    
    
    Parameters
    ----------
    year : integer
    
    month : integer
    
    
    Returns
    -------
    days : integer
    
    
    """
    if month in [1, 3, 5, 7, 8, 10, 12]:
        days = 31
    elif month in [4, 6, 9, 11]:
        days = 30
    else:
        if year % 4 == 0:
            if year % 400 == 0:
                days = 29
            elif year % 100 == 0:
                days = 28
            else:
                days = 29
            #end if
        else:
            days = 28
        #end if
    #end if
    return days
#end func
    
def julian_day(year, month, day):
    return np.sum([days_in_month(int(year), int(m)) \
                   for m in range(1, int(month))]).astype(int) + int(day)
#end func
    
def download_events(start_year, start_month, start_day, end_year, end_month, 
                    end_day, path,):
    """
    Executes a ftp retrieval command to download an event catalogue from the 
    USGS database.
    File is saved in directory 'path' with the name 
    "YYYYWW_cat_quakeml.zip".
    
    
    Parameters
    ----------
    start_year : integer
    
    start_month : integer
    
    start_day : integer
        
    end_year : integer
    
    end_month : integer
    
    end_day : integer
    
    path : string
        Directory in which to save output file.
        
        
    """
    if not os.path.exists(path):
        os.mkdir(path)
    #end if
    
    julian_startday = julian_day(start_year, start_month, start_day)
    start_week = int(julian_startday/7) + 1
    julian_endday = julian_day(end_year, end_month, end_day)
    end_week = np.min([int(julian_endday/7) + 1, 52])
    
    ftp = FTP('hazards.cr.usgs.gov')
    ftp.login()
    ftp.cwd('NEICPDE/quakeml')
    filelist = ftp.nlst()
    
    for filename in filelist:
        YYYYWW = filename.split('_')[0]
        YYYY = int(YYYYWW[:4])
        WW = int(YYYYWW[4:])
        if YYYY == start_year and YYYY == end_year:
            if WW >= start_week and WW <= end_week:
                with open(os.path.join(path, filename), 'wb') as file:
                    ftp.retrbinary('RETR %s'%filename, file.write)
                #end with
            else:
                continue
            #end if
        elif YYYY == start_year:
            if WW >= start_week:
                with open(os.path.join(path, filename), 'wb') as file:
                    ftp.retrbinary('RETR %s'%filename, file.write)
                #end with
            else:
                continue
            #end if
        elif YYYY > start_year and YYYY < end_year:
                with open(os.path.join(path, filename), 'wb') as file:
                    ftp.retrbinary('RETR %s'%filename, file.write)
                #end with
        elif YYYY == end_year:
            if WW <= end_week:
                with open(os.path.join(path, filename), 'wb') as file:
                    ftp.retrbinary('RETR %s'%filename, file.write)
                #end with
            else:
                continue
            #end if
        else:
            continue
        #end if
        print('Downloaded file', filename)
    #end for
    
    try:
        ftp.quit()
    except:
        ftp.close()
    #end try
    
    for filename in glob.glob(os.path.join(path, '*.zip')):
        os.system('unzip %s -d %s'%(filename, path))
    #end for
#end func
    
if __name__ == '__main__':
    process()
#end func