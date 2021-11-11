# Description
This script will download catalogues of events from the ISC web services portal, and produce .xml files in the QuakeML format.

# Instructions
Arguments:
> --bbox: Bounding box, specified as "min_lon min_lat max_lon max_lat".
> --start_date: Start date in format YYYY-MM-DD.
> --start_time: Start time in format hh:mm:ss.
> --end_date: End date in format YYYY-MM-DD.
> --end_time: End time in format hh:mm:ss.
> --mag_range: Magnitude range, specified as "min_mag max_mag".
> --output_path: Directory in which to save output xml files.
> --split_by_day: If true, one xml file will be produced for each day in the time span, or if else, each month in the time span.

Execute from command line, specifying inputs
>
> python download_isc.py --output_path ~ --start_date YYYY-MM-DD --start_time hh:mm:ss --end_date YYYY-MM-DD --end_time hh:mm:ss --bbox min_lon min_lat max_lon max_lat --mag_range min_mag max_mag --split_by_day True
>