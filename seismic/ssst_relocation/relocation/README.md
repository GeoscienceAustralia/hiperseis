# Description
These scripts contain the algorithm for performing Source Specific Station Term (SSST) relocation of a catalogue of earthquakes. The script "main.py" should be run to perform the relocation procedure. The script "tt_table_calculator.py" contains functions to generate pre-computed predicted travel time tables to speed up the relocation procedure. The script "Plots.py" may be run to create diagnostic plots before and after relocation of events. The files "Relocation.py, "Station_Corrections.py", "Travel_Times.py", and "relocation_helper.f" contain various modules used by the relocation procedure.

Workflow
--------
1) Create the travel time tables using "tt_table_calculator.py".

2) Generate a python module from "relocation_helper.f".

3) Perform relocation of earthquakes using "main.py".

4) Convert data to required output format using the script /hiperseis/seismic/ssst_relocation/data_conversion/convert_after_relocation.py and use "Plots.py" to generate diagnostic plots for events before and after relocation.


# relocation_helper.f
This file contains various Fortran subroutines which are converted into python functions for use by the SSST relocation algorithm. To do this, run the following command.

>
> python -m numpy.f2py -c relocation_helper.f -m relocation_helper
>

On linux, this will produce a file relocation_helper.so which can be read in as a python module. On Windows, this will produce a file relocation_helper.pyc.


# tt_table_calculator.py
This script is used to generate predicted travel time tables for a variety of types of P and S waves. The tables may be converted from the iLoc ak135 travel time tables, or generated using the TauP package. A configuration file "tt_table_config.txt" is used to instruct the algorithm how to construct the tables.


Arguments
---------
--config_file: name of configuration file.

--output_path: desired output directory.

--from_iloc_tables: set to True if a user wishes to convert the iLoc ak135 tables into the required format. If so, specify --iloc_path also.

--iloc_path: directory containing iLoc travel time tables, if required.

--model: name of model to use if using the TauP package for table generation, e.g. 'iasp91'.


Configuration
-------------
In the configuration file, for each table, three pieces of information are needed.
1) Phase(s) to compute travel times for.
2) Epicentral distance samples to use.
3) Depth samples to use.

This should be in the format of the file tt_table_config.txt.


Usage
-----
>
> python tt_table_calculator.py --config_file <configuration file> --output_path .../output_path/ --from_iloc_tables True --iloc_path .../iloc_path/
>
> python tt_table_calculator.py --config_file <configuration file> --output_path .../output_path/ --model <model name>
>


# main.py
This script performs the relocation procedure. The procedure is as follows.

1) Read pick information from input file created using /hiperseis/seismic/ssst_relocation/data_conversion/convert_for_relocation.py. 

2) For each pick, calculate predicted travel times for all available phases.

3) Redefine phase for each pick as that corresponding to the smallest travel time residual.

4) Compute time corrections for each station. These may either be static stationterms or source specific station terms. Apply the corrections to the picks.

5) Compute new best fitting hypocentres using minimisation algorithm. 

6) Reject new hypocentres if they differ by more than a certain threshold distancefrom original hypocentre. Otherwise update hypocentres for events.

The use of static station terms (SSTs) when performing earthquake relocation can account for differences in seismic travel times caused by heterogeneity in the region about a receiver. This method applies the median travel time residual of all arrivals at a station as a time correction.
Using source specific station terms (SSSTs) can improve the results further by accounting for heterogeneity along the path travelled between a source and a receiver. In this method a spatially varying time correction is applied to the picks. 
For the method employed in this workflow, time corrections are calculated in the following way.

1) Find all arrivals at a station of a particular phase.

2) For each event corresponding to these arrivals, find all other events within a threshold distance.

3) Compute the median travel time residual for the event and its nearest neighbours(within threshold distance).

4) Apply this median residual as a correction to the arrival time corresponding to this event.


Arguments
---------
--config_file: Name of configuration file to use.
    
--tt_table_path: Directory containing pre-computed travel time tables.
  
--input_path: Directory containing pick array save file and, if iteration > 0, list of events with unstable hypocentres.
  
--elcordir: Directory containing ellipticity correction coefficient tables.
 
--iteration: Iteration number.
 
--output_path: Directory to output updated pick array and list of events with unstable hypocentres

--relocation_algorithm: If relocation_algorithm is 'iloc' then the iLoc relocation algorithm is used, or if not, the relocation algorithm in relocation_helper.f is used.


Configuration
-------------
The configuration file for main.py should be that contained in relocation_config.txt. This has the following variables.
hypo_thr_dist_deg: float, representing the allowed distance a hypocentre may move during relocation before being classed as unstable.

hypo_thr_time_sec: float, representing the allowed number of seconds an event origin time may change by during relocation before being classed as unstable.

corr_thr_dist_deg: float, radius of sphere (in degrees) used for computing SSST corrections.

thr_p: float, time threshold for P wave phase redefinition.

thr_s: float, time threshold for S wave phase redefinition.

no_relocation: boolean. Set to True if only phase redefinition is desired.

correction_method: string. Set to 'SSST' if SSST corrections are desired, 'SST' if static station term corrections are desired, or 'None' for no time corrections.

phases: list, describing which phases to try to match arrivals to in the first pass of the phase redefinition algorithm. Additionally, if the default relocation algorithm is used, only these phases are used for relocation. Example: phases = P, Pg, Pn, Pb, S, Sg, Sn, Sb.

If the iLoc relocation algorithm is to be used, also use the following variables:
iloc_redefine_phases: boolean. If True, iLoc will be allowed to redefine phases during its relocation procedure.

iloc_use_rstt: boolean. If True, iLoc will be allowed to use the RSTT model for travel time calculation for Pg, Sg, Pn, and Sn phases.

If the default relocation algorithm is to be used, also use the following variables:
reloc_dlat: float, describing the initial spacing between grid points for the grid search algorithm. This decreases by 'reloc_sfrac' amount each iteration.

reloc_ddep: float, decribing the initial depth spacing between grid points for the grid search algorithm. This also decreases by 'reloc_sfrac' amount each iteration.

reloc_nit: integer, number of iterations to use for default relocation algorithm.

reloc_norm: integer. Set to 1 if L1 norm is used for minimisation routine, or 2 if L2 norm (RMS) is to be used.

reloc_sfrac: float, describing the proportion that reloc_dlat and reloc_ddep will shrink by each iteration. E.g. if reloc_sfrac = 0.5, then on each iteration, reloc_dlat and reloc_ddep are halved.

temp_networks: list, describing temporary seismic station deployments which are not to be used for relocation. Picks on these networks still have their phases redefined and residuals calculated. For example: temp_networks = 7B, 7D, 7E, 7F, 7G, 7J, 7Q, 7T, 7U, 7W, 7X, OA, W1, AQ.


Usage
-----
Run multiple iterations of this algorithm using the command below, e.g. by executing "for i in `seq 1 $number_of_iterations`; do <command>; done". On the first iteration, the algorithm only finds hypocentres which are unstable, using no time corrections. On all subsequent iterations, the configuration specified is followed.

>
> mpirun -np $number_of_processors main.py --config_file <configuration file> --tt_table_path .../tt_table_path/ --input_path .../input_path/ --elcordir .../elcordir/ --iteration $i --output_path .../output_path/
>