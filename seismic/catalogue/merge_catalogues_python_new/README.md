# Description
These scripts will take catalogues of events and produce a single merged list of all events in the catalogue, and their picks. 

Workflow
--------
1) On multiple processors, read all input .xml files from Engdahl, ISC, USGS, and GA analyst reviewed datasets, and produce ‘pickle’ binary files to store the output. 

2) On a single processor, read in the binary files, and read in the hypocentres from Gary Gibson catalogues. Filter out duplicate events, and filter the catalogue to ensure the azimuthal gap, magnitude, and bounding box criteria are met. Save each event in a separate ‘pickle’ binary file.

3) On multiple processors, read in the binary files and convert the information into the .csv output format. Save the output from each processor to disk as a .csv file.

4) On a single processor, read in the .csv files, merge them, and produce a single .csv output file storing all events and their associated picks.


# Read_XML_Multiprocessing.py
This script is used to shorten the time required to read in catalogues from xml files. On multiple processors, catalogues are read in and merged, then output as a .pkl file. The output is to be read by Filter_Catalogues.py.


Arguments
---------
- paths: List of directories to read .xml files from. All subdirectories of these paths are also read recursively.

- already_read_path: Output directory from any previous run of the script. If --split_smaller was used in a previous run, the output file names (less file extension) are identical to the input file names, so the script will search the specified directory for these files. If it finds a match, the corresponding xml file is not required to be processed.

- output_path: Directory in which to save output .pkl files.

- split_smaller: If true, one .pkl file will be produced for each .xml file in the input directories, or else, one .pkl file is created per processor.


Usage
-----
>
> mpirun -np $nprocessors Read_XML_Multiprocessing.py --output_path ~ --already_read_path ~ --paths .../path1/ .../path2 ... --split_smaller True
>


# Read_USGS_XML_Multiprocessing.py
This script is almost identical to Read_USGS_XML_Multiprocessing, however since it is impossible to filter data from the USGS ftp server, a filtering routine is included to remove events which don't fit certain bounding box restrictions (see Filter_Catalogue.py description).


Arguments
---------
- paths: List of directories to read .xml files from. All subdirectories of these paths are also read recursively.

- already_read_path: Output directory from any previous run of the script. If --split_smaller was used in a previous run, the output file names (less file extension) are identical to the input file names, so the script will search the specified directory for these files. If it finds a match, the corresponding xml file is not required to be processed.

- output_path: Directory in which to save output .pkl files.

- split_smaller: If true, one .pkl file will be produced for each .xml file in the input directories, or else, one .pkl file is created per processor.


Usage
-----
>
> mpirun -np $nprocessors Read_USGS_XML_Multiprocessing.py --output_path ~ --already_read_path ~ --paths .../path1/ .../path2 ... --split_smaller True
>


# Filter_Catalogues.py
This script is used to read the output of Read_XML_Multiprocessing.py on a single processor, then merge and filter the catalogues.


Criteria
--------
- Events within the Australian bounding box [105, -45, 160, -13] are filtered to have azimuthal gap less than 180 degrees.

- Events in the Indian Ocean in bounding boxes [15, -60, 90, 30] and [90, -60, 180, 45] are filtered to have a magnitude larger than 5.

- Events in the rest of the world are filtered to have magnitude larger than 5.5.

- Events within +/- 15 seconds and within 50km of one another are counted as duplicates, and only one is kept.

- Hypocentres sourced from Gary Gibson are used to replace those in the GA catalogues if a duplicate is found.


Arguments
---------
- -GA, --GA_path: Directory in which GA event .xml files are stored.

- -EHB, --EHB_path: Directory in which Engdahl event .xml files are stored.

- -ISC, --ISC_path: Directory in which ISC event .xml files are stored.

- -PP, --preprocessed_path: The output directory for Read_XML_Multiprocessing.py.

- -GG, --GG_path: Directory in which Gary Gibson event .csv files are stored.

- -USGS, --USGS_path: Directory in which USGS event .csv files are stored.

- -O, --output_path: Directory in which to save output files.

- --read_from_preprocessed: True if xml files have been pre-processed using Read_XML_Multiprocessing.py, in which case .pkl files are read from preprocessed_path. If else, .xml files are read from GA_path and/or ISC_path and/or EHB_path.


Usage
-----
>
> python3 Filter_Catalogues.py --preprocessed_path .../preprocessed_path/ --GG_path .../GG_path/ --output_path .../output_path/ --read_from_preprocessed True/False
> 
> python3 Filter_Catalogues.py --GA_path .../GA_path/ --ISC-path .../ISC_path/ --EHB_path .../EHB_path/ --GG_path .../EHB_path/ --USGS_path .../USGS_path/ --output_path .../output_path/ 
>


# Convert_Catalogues_To_CSV.py
This script is used to read the output of Filter_Catalogues and convert the catalogues to .csv format on multiple processors. Using a single processor to do this takes days, while using multiple processors to filter and convert the catalogues wastes days worth of computation time while input files are being read in by Filter_Catalogues.py.


Arguments
---------
- -PP, --preprocessed_path: The output directory for Filter_Catalogues.py.

- -ST, --station_path: Directory where station xml file/s is/are stored.

- -O, --output_path: Directory in which to save output files.


Usage
-----
>
> mpirun -np $nprocessors Convert_Catalogues_CSV.py --preprocessed_path ~ --station_path ~ --output_path ~
>


# Merge_Catalogues.py
This script is used to merge the .csv output of Convert_Catalogues_To_CSV.py, and re-order the events by date.


Arguments
---------
- -I, --input_path: Directory with event .csv files to be merged.

- -O, --output_file: Name of required output file.


Usage
-----
>
> python3 Order_csv_output_by_date.py --input_file ~ --output_file ~
>
