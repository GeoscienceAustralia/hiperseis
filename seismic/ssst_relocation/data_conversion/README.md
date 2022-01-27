# Description
The scripts convert_for_relocation.py and convert_array_to_catalogue.py will convert input data into the format required by the Source Specific Station Term (SSST) relocation algorithm. Once relocation has been completed, the script convert_after_relocation.py will convert the output file into a format usable by the tomographic inversion software.

Workflow
--------
1) On a single processor, convert earthquake catalogue (output of /hiperseis/seismic/catalogue/merge_catalogues/python_new/merge_catalogues.py) and the automatic pick harvester output (output of /hiperseis/seismic/pick_harvester/pick.py) into numpy binary file format using convert_for_relocation.py.

2) If needed, on multiple processors, convert output of convert_for_relocation.py into 'SC3ML' .xml files for use by SeisComp3 using convert_array_to_catalogue.py.

3) Run the relocation algorithm.

4) Convert the output of the relocation algorithm (/hiperseis/seismic/ssst_relocation/relocation/main.py) into ASCII format using convert_after_relocation.py.


# convert_for_relocation.py

Script to convert earthquake catalogue and harvested picks from temporary networks into numpy binary format.


Arguments
---------
--event_file: name of output file from /hiperseis/seismic/catalogue/merge_catalogues/python_new/merge_catalogues.py.

--p_combined: name of P wave output file from /hiperseis/seismic/pick_harvester/pick.py.

--s_combined: name of S wave output file from /hiperseis/seismic/pick_harvester/pick.py.

--output_path: desired output directory.


Usage
-----
>
> python convert_for_relocation.py --event_file .../event_file.csv --p_combined .../p_combined.txt --s_combined .../s_combined.txt --output_path .../output_path/
>


# convert_array_to_catalogue.py
If using the iLoc earthquake relocation algorithm, a SeisComp3 database is required for data storage and input. This script will convert the output of convert_for_relocation.py into the correct format for upload to this database. Note that the output of conver_for_relocation.py is still used by the SSST relocation algorithm. 

Configuration
-------------
If a user wishes for iLoc to ignore certain networks when performing the relocation with iLoc, these must be specified in a configuration file. Two pieces of information are needed.
1) ignore_temp_networks: boolean. True if the user wishes to ignore certain networks.
2) temp_networks: list, with separator ', '.
The format for this configuraton file is the same as the example network_config.txt file.


Arguments
---------
--event_file: Name of file containing array of picks.

--config_file: Name of configuration file telling algorithm which networks to ignore.

--output_path: Output path.

--output_format: Format to use for the output .xml file. May be 'SC3ML', 'QuakeML', or many others.


Usage
-----
>
> mpirun -np $number_of_processors python convert_array_to_catalogue.py --event_file <event_file.npy> --config_file <network_config.txt> --output_path .../output_path/ --output_format <output_format>
>


# convert_after_relocation.py
After relocation, this script will convert the output numpy binary file into a format usable by the tomographic inversion software.


Arguments
---------
--event_file: Name of numpy binary file to extrack picks from.

--output_path: Output path.

--include_tcor: If a column is required for time corrections, set include_tcor = True. Note: this should only be used for plotting purposes. The tomographic inversion software cannot read the file if the time correction column is included.

    
Usage
-----
>
> python convert_after_relocation.py --event_file .../event_file.npy --output_path .../output_path/ --include_tcor True
>