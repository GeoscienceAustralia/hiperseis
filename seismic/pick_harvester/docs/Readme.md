# SSST Workflow: Developer Documentation

This readme provides a high-level overview of the structure and design of the various components of the SSST 
workflow.

## Module Structure

![module](module.svg?raw=true "Title")

The workflow comprises three scripts: i) `pick.py`, ii) `ssst_relocate.py` and iii) `cluster_arrivals.py` -- refer to 
user documentation for further details on how to run these scripts [here](../Readme.md).

## Class Structure

![module](classes.svg?raw=true "Title")

`SSSTRelocator` is a subclass of `ParametricData`, inheriting generic parallelization approaches implemented in the 
baseclass -- some of the key features are outlined in the following section. `SSSTRelocator` also includes a number 
of automatically generated member variables (not shown in the class-diagram) that categorizes input events and 
arrivals, enabling more efficient and customizable processing. For example, if the input catalog comprises events 
 and arrivals from two sources (GA, GG), including automatically generated picks, the following member variables are 
generated internally:

1. `is_GA_event` (boolean, `dim(nEvents))`
2. `is_GG_event` (boolean, `dim(nEvents))`
3. `is_GA_arrival` (boolean, `dim(nArrivals))`
4. `is_GG_arrival` (boolean, `dim(nArrivals))`
5. `is_AUTO_arrival` (boolean, `dim(nArrivals))`

`TTInterpolator` is responsible for interpolating travel-time and associated derivatives. `Ellipticity` is responsible 
for computing ellipticity corrections. Note that `SSST_Result` also has to use the `Ellipticity` class, since, in the 
interests of saving disk-space, the hdf5 output of `ssst_relocate.py` does not contain ellipticity corrections at each 
iteration of the relocation procedure.

`NestedGrid` embodies a coarse-resolution, global outer `Grid` within which a finer inner `Grid` defines a region of 
interest for the clustering procedure in `cluster_arrivals.py`.

## Parallelization
MPI-collectives (e.g. AllGather) are not suitable for syncing large amounts of data across processes -- a number of 
 these failure modes encountered while developing this workflow are documented 
[here](https://ieeexplore.ieee.org/document/8425443). 
Disk-based `_sync` and `_gather` functions have been implemented in `ParametricData` and `SSSTRelocator` classes, 
respectively. `_sync` collects a fixed number of local elements of a global array on each rank and propagates them 
to every other rank  -- e.g. travel-time residuals for subsets of arrivals are computed on each rank in parallel and 
then synced across all ranks. `_gather` on the other hand collects a variable number of local elements and their 
respective indices within a global array on each rank and propagates them to every other rank -- e.g. when SSST 
corrections are computed, only a subset of local arrivals may participate in the process, e.g. when arrivals from 
certain high-quality sources are set not to have their phases redefined. 


