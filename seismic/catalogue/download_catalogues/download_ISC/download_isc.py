"""
Description:
    This script will download a catalogue of events from the ISC database given
    various input criteria by the user. An output file will be generated for 
    each month or each day in the requested time span (depending on argument 
    "--split_by_day".
    
Developer: Lachlan Adams 
Contact: lachlan.adams@ga.gov.au or lachlan.adams.1996@outlook.com

"""

import argparse, os
import numpy as np
from obspy import UTCDateTime

def process():
    """
    Download ISC catalogue fulfilling specified criteria, splitting into a
    catalogue for each month.
    
    
    Usage
    -----
    python download_isc.py --bbox min_lon min_lat max_lon max_lat
    --start_date YY-MM-DD --start_time hh:mm:ss --end_date YYYY-MM-DD 
    --end_time hh:mm:ss --mag_range min_mag max_mag --output_path 
    .../output_path --split_by_day
    
    
    """
    parser = argparse.ArgumentParser(description='ISC event downloader')

    parser.add_argument("--bbox", nargs=4, type=float, default=None, 
                        help="Enter as '--bbox min_lon min_lat max_lon max_lat'")
    parser.add_argument("--start_date", type=str, default=None,
                        help="Enter as '--start_date YYYY-MM-DD'")
    parser.add_argument("--start_time", type=str, default=None,
                        help="Enter as '--start_time hh:mm:ss'")
    parser.add_argument("--end_date", type=str, default=None,
                        help="Enter as '--end_date YYYY-MM-DD'")
    parser.add_argument("--end_time", type=str, default=None,
                        help="Enter as '--end_time hh:mm:ss'")
    parser.add_argument("--mag_range", nargs=2, type=float, default=None,
                        help="Enter as '--mag_range min_mag max_mag'")
    parser.add_argument("--output_path", type=str, default=None)
    parser.add_argument("--split_by_day", type=bool, default=False)
    
    args = parser.parse_args()
    output_path = args.output_path
    bbox = args.bbox
    start_date = args.start_date
    start_time = args.start_time
    end_date = args.end_date
    end_time = args.end_time
    mag_range = args.mag_range
    split_by_day = args.split_by_day
    
    if output_path is None:
        output_path = input("Output path? ")
    if bbox is None:
        bbox = input("Bounding box (min_lon, min_lat, max_lon, max_lat)? ")
        min_lon, min_lat, max_lon, max_lat = \
            (float(item) for item in bbox.split(','))
    else:
        min_lon, min_lat, max_lon, max_lat = bbox
    if start_date is None:
        start_date = input("Start date (YYYY-MM-DD)? ")
    if start_time is None:
        start_time = input("Start time (hh:mm:ss)? ")
    if end_date is None:
        end_date = input("End date (YYYY-MM-DD)? ")
    if end_time is None:
        end_time = input("End time (hh:mm:ss)? ")
    if mag_range is None:
        mag_range = input("Magnitude range (min_mag, max_mag)? ")
        min_mag, max_mag = (float(item) for item in mag_range.split(','))
    else:
        min_mag, max_mag = mag_range
    
    start_year, start_month, start_day = \
        (int(item) for item in start_date.split('-'))
    end_year, end_month, end_day = \
        (int(item) for item in end_date.split('-'))
    
    if split_by_day:
        get_events_by_day(min_lon, min_lat, max_lon, max_lat, start_date, 
                          start_time, end_date, end_time, min_mag, max_mag, 
                          output_path)
    else:
        get_events_by_month(min_lon, min_lat, max_lon, max_lat, start_year, 
                            start_month, start_day, start_time, end_year, 
                            end_month, end_day, end_time, min_mag, max_mag, 
                            output_path)
    #end if
    
#end func
    
def get_events_by_day(min_lon, min_lat, max_lon, max_lat, start_date, 
                      start_time, end_date, end_time, min_mag, max_mag, 
                      output_path):
    st = UTCDateTime(str(start_date + 'T' + start_time))
    et = UTCDateTime(str(end_date + 'T' + end_time))
    
    time_range = et - st
    ints = np.hstack([np.arange(0, time_range, 86400), np.array([time_range])])
    
    days = [st + ints[i] for i in range(len(ints))]
    start_times = days[:-1]
    end_times = days[1:]
    for day in range(len(days)-1):
        start_year = start_times[day].year
        start_month = start_times[day].month
        start_day = start_times[day].day
        start_hour = start_times[day].hour
        start_minute = start_times[day].minute
        start_second = start_times[day].second
        end_year = end_times[day].year
        end_month = end_times[day].month
        end_day = end_times[day].day
        end_hour = end_times[day].hour
        end_minute = end_times[day].minute
        end_second = end_times[day].second
        start_time = '%s:%s:%s'%(str(start_hour).zfill(2), 
                                 str(start_minute).zfill(2), 
                                 str(start_second).zfill(2))
        end_time = '%s:%s:%s'%(str(end_hour).zfill(2), 
                               str(end_minute).zfill(2), 
                               str(end_second).zfill(2))
        download_events(min_lon, min_lat, max_lon, max_lat, start_year, 
                        start_month, start_day, start_time, end_year, 
                        end_month, end_day, end_time, min_mag, max_mag, 
                        output_path)
    #end for
#end func    

def get_events_by_month(min_lon, min_lat, max_lon, max_lat, start_year, 
                        start_month, start_day, start_time, end_year, 
                        end_month, end_day, end_time, min_mag, max_mag, 
                        output_path):
    if end_year == start_year:
        year = start_year
        if end_month == start_month:
            month = start_month
            download_events(min_lon, min_lat, max_lon, max_lat, year, month, 
                            start_day, start_time, year, month, end_day, 
                            end_time, min_mag, max_mag, output_path)
        else:
            months = np.arange(start_month, end_month+1)
            for month in months:
                max_days = days_in_month(year, month)
                if month == start_month:
                    download_events(min_lon, min_lat, max_lon, max_lat, year, 
                                    month, start_day, start_time, year, 
                                    month, max_days, '23:59:59', min_mag, 
                                    max_mag, output_path)
                elif month == end_month:
                    download_events(min_lon, min_lat, max_lon, max_lat, year, 
                                    month, 1, '00:00:00', year, month, end_day, 
                                    end_time, min_mag, max_mag, output_path)
                else:
                    download_events(min_lon, min_lat, max_lon, max_lat, year, 
                                    month, 1, '00:00:00', year, month, 
                                    max_days, '23:59:59', min_mag, max_mag, 
                                    output_path)
                #end if
            #end for
        #end if
    else:
        years = np.arange(start_year, end_year+1)
        for year in years:
            if year == start_year:
                months = np.arange(start_month, 13)
                for month in months:
                    max_days = days_in_month(year, month)
                    if month == start_month:
                        download_events(min_lon, min_lat, max_lon, max_lat, 
                                        year, month, start_day, start_time, 
                                        year, month, max_days, '23:59:59', 
                                        min_mag, max_mag, output_path)
                    else:
                        download_events(min_lon, min_lat, max_lon, max_lat, 
                                        year, month, 1, '00:00:00', year,
                                        month, max_days, '23:59:59', min_mag, 
                                        max_mag, output_path)
                    #end if
                #end for
            elif year == end_year:
                months = np.arange(1, end_month+1)
                for month in months:
                    max_days = days_in_month(year, month)
                    if month == end_month:
                        download_events(min_lon, min_lat, max_lon, max_lat, 
                                        year, month, 1, '00:00:00', year, 
                                        month, end_day, end_time, min_mag, 
                                        max_mag, output_path)
                    else:
                        download_events(min_lon, min_lat, max_lon, max_lat, 
                                        year, month, 1, '00:00:00', year, 
                                        month, max_days, '23:59:59', min_mag, 
                                        max_mag, output_path)
                    #end if
                #end for
            else:
                months = np.arange(1, 13)
                for month in months:
                    max_days = days_in_month(year, month)
                    download_events(min_lon, min_lat, max_lon, max_lat, year, 
                                    month, 1, '00:00:00', year, month, 
                                    max_days, '23:59:59', min_mag, max_mag,
                                    output_path)
                #end for
            #end if
        #end for
    #end if
#end func    
            
def download_events(min_lon, min_lat, max_lon, max_lat, start_year, 
                    start_month, start_day, start_time, end_year, end_month, 
                    end_day, end_time, min_mag, max_mag, path):
    """
    Executes a 'wget' command to download an event catalogue from the ISC
    database.
    File is saved in directory 'path' with the name 
    "YYYY-MM-DDThh:mm:ss_YYYY-MM-DDThh:mm:ss.xml", i.e. start date/time 
    followed by end_time.
    
    
    Parameters
    ----------
    min_lon : float
    
    min_lat : float
    
    max_lon : float
    
    max_lat : float
    
    start_year : integer
    
    start_month : integer
    
    start_day : integer
    
    start_time : string
        String in format 'hh:mm:ss'.
        
    end_year : integer
    
    end_month : integer
    
    end_day : integer
    
    end_time : string
        String in format 'hh:mm:ss'.
    
    min_mag : float
    
    max_mag : float
    
    path : string
        Directory in which to save output file.
        
        
    """
    if not os.path.exists(path):
        os.mkdir(path)
    #end if
    
    start_hour, start_minute, start_second = start_time.split(':')
    end_hour, end_minute, end_second = end_time.split(':')
    
    filename = '%s-%s-%sT%s-%s-%s_%s-%s-%sT%s-%s-%s.xml' \
                %(start_year, str(start_month).zfill(2), 
                  str(start_day).zfill(2), str(start_hour).zfill(2),
                  str(start_minute).zfill(2), str(start_second).zfill(2),
                  end_year, str(end_month).zfill(2), str(end_day).zfill(2), 
                  str(end_hour).zfill(2), str(end_minute).zfill(2), 
                  str(end_second).zfill(2))
    file = os.path.join(path, filename)
    
    wget_string = str("wget --output-document=%s "%file + \
                      "'http://www.isc.ac.uk/cgi-bin/web-db-v4?" + \
                      "out_format=QuakeML&" + \
                      "request=STNARRIVALS&" + \
                      "ttime=on&" + \
                      "iscreview=on&" + \
                      "stnsearch=GLOBAL&" + \
                      "searchshape=RECT&" + \
                      "bot_lat=%s&"%min_lat + \
                      "top_lat=%s&"%max_lat + \
                      "left_lon=%s&"%min_lon + \
                      "right_lon=%s&"%max_lon + \
                      "start_year=%s&"%start_year + \
                      "start_month=%s&"%start_month + \
                      "start_day=%s&"%start_day + \
                      "start_time=%s&"%start_time + \
                      "end_year=%s&"%end_year + \
                      "end_month=%s&"%end_month + \
                      "end_day=%s&"%end_day + \
                      "end_time=%s&"%end_time + \
                      "min_mag=%s&"%min_mag + \
                      "max_mag=%s'"%max_mag)
    os.system(wget_string)
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
    
if __name__ == '__main__':
    process()
#end func