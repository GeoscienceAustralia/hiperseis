"""
Description
-----------
This script will download a catalogue of events from the USGS database given 
various input criteria by the user. 
    
Developer: Lachlan Adams 
Contact: lachlan.adams@ga.gov.au or lachlan.adams.1996@outlook.com

"""

import argparse, os

def process():
    """
    Download USGS catalogue fulfilling specified criteria.
    
    
    Usage
    -----
    python download_usgs.py --bbox min_lon min_lat max_lon max_lat
    --start_date YY-MM-DD --start_time hh:mm:ss --end_date YYYY-MM-DD 
    --end_time hh:mm:ss --mag_range min_mag max_mag --output_path 
    .../output_path --output_format format
    
    
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
    parser.add_argument("--output_path", type=str, default='.')
    parser.add_argument("--output_format", type=str, default='csv')
    
    args = parser.parse_args()
    output_path = args.output_path
    bbox = args.bbox
    start_date = args.start_date
    start_time = args.start_time
    end_date = args.end_date
    end_time = args.end_time
    mag_range = args.mag_range
    output_format = args.output_format
    
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
    
    download_events(min_lon, min_lat, max_lon, max_lat, start_year, 
                    start_month, start_day, start_time, end_year, end_month, 
                    end_day, end_time, min_mag, max_mag, output_path, 
                    output_format=output_format)    
#end func
    
def download_events(min_lon, min_lat, max_lon, max_lat, start_year, 
                    start_month, start_day, start_time, end_year, end_month, 
                    end_day, end_time, min_mag, max_mag, path, 
                    output_format='csv'):
    """
    Executes a 'wget' command to download an event catalogue from the USGS
    database.
    File is saved in directory 'path' with the name 
    "YYYY-MM-DDThh:mm:ss_YYYY-MM-DDThh:mm:ss.csv", i.e. start date/time 
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
        
    output_format : string
        'csv', 'xml', 'text', etc.
        
        
    """
    if not os.path.exists(path):
        os.mkdir(path)
    #end if
    
    start_hour, start_minute, start_second = start_time.split(':')
    end_hour, end_minute, end_second = end_time.split(':')
    
    filename = '%s-%s-%sT%s-%s-%s_%s-%s-%sT%s-%s-%s.csv' \
                %(start_year, str(start_month).zfill(2), 
                  str(start_day).zfill(2), str(start_hour).zfill(2),
                  str(start_minute).zfill(2), str(start_second).zfill(2),
                  end_year, str(end_month).zfill(2), str(end_day).zfill(2), 
                  str(end_hour).zfill(2), str(end_minute).zfill(2), 
                  str(end_second).zfill(2))
    file = os.path.join(path, filename)
    start_time = '%s-%s-%sT%s:%s:%s' \
                %(start_year, str(start_month).zfill(2), 
                  str(start_day).zfill(2), str(start_hour).zfill(2),
                  str(start_minute).zfill(2), str(start_second).zfill(2))
    end_time = '%s-%s-%sT%s:%s:%s' \
                %(end_year, str(end_month).zfill(2), str(end_day).zfill(2), 
                  str(end_hour).zfill(2), str(end_minute).zfill(2), 
                  str(end_second).zfill(2))
    
    wget_string = str("wget --output-document=%s "%file + \
                      "'https://earthquake.usgs.gov/fdsnws/event/1/query?" + \
                      "format=%s&"%output_format + \
                      "starttime=%s&"%start_time + \
                      "endtime=%s&"%end_time + \
                      "minlatitude=%s&"%min_lat + \
                      "minlongitude=%s&"%min_lon + \
                      "maxlatitude=%s&"%max_lat + \
                      "maxlongitude=%s&"%max_lon + \
                      "minmagnitude=%s&"%min_mag + \
                      "maxmagnitude=%s'"%max_mag)
    os.system(wget_string)
#end func
    
if __name__ == '__main__':
    process()
#end func