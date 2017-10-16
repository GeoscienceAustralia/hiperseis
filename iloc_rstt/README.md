From Seiscomp3 events to location using iLoc
============================================

Our workflow contains the following steps:

1. Export event xml from the seiscomp3 db. These events initially only have a 
few primary stations associations with mostly P arrivals. We will use the `ISF`
 functionality of [iLoc](http://www.seismology.hu/index.php/en/home/iloc). 
 Seiscomp3 provides executable (`sccnv`) that allows us to export events into 
 ISF format. This is how it works:
    ```bash
    sccnv -i event.xml -o ims10:isf.txt
    ```
2. Record the location (lat/lon) of the event and export primary stations 
waveforms (all channels) from all stations within a certain radius of the 
event. This radius will be, say, 20% more than the farthest (from event 
location) primary station already associated.

3. If any temporary station falls within this radius, we need to add one or max 
two temporary stations waveform (from each temporary stations network).

4. For a primary station for which we already have association, this station 
will have most likely have an associated P arrival, or may be an S arrival 
(confirm). Make sure the vertical channel is associated with the P wave. 
(Some rules required here if the P is associated with some other 
channels). If we have the S arrival, we are done with that station after 
confirming that the S waves came from either the `N` or `E` channels, where 
`E` and `N` are the North-South, and East-West components at the station.

5. Otherwise, we need to detect S wave arrival for this station. The 
procedure is as follows. We need to take a 50 seconds window (25 seconds 
either side) of the P arrival time and run the (S wave) picking algorithm. Note 
the time of the S arrival and write it (i.e. associate the arrival) in the 
`ISF` file. We will use `(sqrt(E^2 + N^2))` composite waveform for S waves.

6. For primary stations for which we don't already have associations (possible
 due to step 2 above), we need to run both the P and S picking algorithms and 
 add the picks in the `ISF` file. 

7. We repeat step 6 for temporary station identified in step 3.

8. We run `iloc` on the updated `ISF` file as the following:
    ```bash
    echo "isf_infile=isf.txt isf_outfile=isf.out" | iloc isf > isf.log
    ```
    
When an arrival is added in the ISF file, and used in `iloc`, `iloc` 
assumes that the arrival is associated with the event, and runs it's location
 algorithm on all the associated arrivals in the `ISF` file.  