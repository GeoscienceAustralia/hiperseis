# Description
These scripts will download catalogues of events from the USGS web services portal.


# download_usgs_hypocentres.py
This script will download hypocentres for events on the USGS database.


Arguments
---------
- bbox: Bounding box, specified as "min_lon min_lat max_lon max_lat".

- start_date: Start date in format YYYY-MM-DD.

- start_time: Start time in format hh:mm:ss.

- end_date: End date in format YYYY-MM-DD.

- end_time: End time in format hh:mm:ss.

- mag_range: Magnitude range, specified as "min_mag max_mag".

- output_path: Directory in which to save output files.

- output_format: Desired save file format. Possible values are 'csv', 'xml', 'text'.


Usage
-----
>
> python download_usgs_hypocentres.py --output_path ~ --start_date YYYY-MM-DD --start_time hh:mm:ss --end_date YYYY-MM-DD --end_time hh:mm:ss --bbox min_lon min_lat max_lon max_lat --mag_range min_mag max_mag --output_format format
>


# download_usgs_quakeml.py
This script will download quakeml files for events on the USGS database. Events on this database are stored in .zip files spanning a week each. There is no way to filter the returned events except by origin time. The script connects to ftp://hazards.cr.usgs.gov/ and extracts events from directory /NEICPDE/quakeml/.


Arguments
---------
- start_date: Start date in format YYYY-MM-DD.

- end_date: End date in format YYYY-MM-DD.

- output_path: Directory in which to save output xml files.


Usage
-----
>
> python download_usgs_catalogue.py --output_path ~ --start_date YYYY-MM-DD --end_date YYYY-MM-DD
>