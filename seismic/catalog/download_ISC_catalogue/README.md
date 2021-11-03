# Instructions
If executing from the command line, arguments aren't required, however you will be prompted to enter them at various points.


Execute from command line without specifying inputs
>
> python download_isc.py 
>Output path? ~
>Bounding box (min_lon, min_lat, max_lon, max_lat)? min_lon, min_lat, max_lon, max_lat
>Start date (YYYY-MM-DD)? YYYY-MM-DD
>Start time (hh:mm:ss)? hh:mm:ss
>End date (YYYY-MM-DD)? YYYY-MM-DD
>End time (hh:mm:ss)? hh:mm:ss
>Magnitude range (min_mag, max_mag)? min_mag, max_mag


Execute from command line, specifying inputs
>
> python download_isc.py --output_path ~ --start_date YYYY-MM-DD --start_time hh:mm:ss --end_date YYYY-MM-DD --end_time hh:mm:ss --bbox min_lon min_lat max_lon max_lat --mag_range min_mag max_mag
>