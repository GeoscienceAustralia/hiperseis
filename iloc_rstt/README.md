From Seiscomp3 events to location using iLoc
============================================

## Using ISF file approach

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
event. This radius will be, say, 20% more than the farthest primary station 
(from the event location) already associated.

3. If any temporary station falls within this radius, we need to add one or max 
two temporary stations waveform (from each temporary stations network).

4. For a primary station for which we already have association, this station 
will most likely have an associated P arrival, or may be an S arrival 
(confirm). Make sure the vertical channel is associated with the P wave. 
(Some rules required here if the P is associated with some other 
channels). If we have the S arrival, we are done with that station after 
confirming that the S waves came from either the `N` or `E` channels, where 
`N` and `E` are the North-South and East-West components at the station.

5. Otherwise, we need to detect S wave arrival for this station. The 
procedure is as follows. We take a 50 seconds window (25 seconds either side)
 of the P arrival time and run the (S wave) picking algorithm. Note the time 
 of the S arrival and write it (i.e. associate the arrival) in the `ISF` 
 event file. We will use `(sqrt(E^2 + N^2))` composite waveform for S wave 
 detection.

6. For primary stations for which we don't already have associations (possible
 due to step 2 above), we need to run both the P and S picking algorithms and 
 add the picks in the `ISF` file. This is tricky since we have no reference 
 arrivals for these station and theoretical 1d arrival times may not exactly 
 match that of from the observed waveforms at the stations. We have to run 
 our P wave detection algorithm first on a long time window. This process can
  produce multiple picks and we need a process to select/reject such picks. 
  We can possibly follow a [`PhasePApy`](https://github.com/GeoscienceAustralia/PhasePApy) style pick filtering approach. 
  Once the P wave is detected, we can follow the same approach as step 5 for 
  the S wave detection. 

7. We repeat step 6 for temporary station identified in step 3.

8. We run `iloc` on the updated `ISF` file as the following:
    ```bash
    echo "isf_infile=isf.txt isf_outfile=isf.out" | iloc isf > isf.log
    ```
    
When an arrival is added in the `ISF` file, and used in `iloc`, `iloc` 
assumes that the arrival is associated with the event, and runs it's location
 algorithm on all the associated arrivals in the `ISF` file.
 
## Alternative approach to ISF: Use iloc with seiscomp3 db

According to the `iloc` manual, we can also use events directly from seiscomp3
 db in the following manner:
 
     ```bash
     echo "event_id update_db=1 do_gridsearch=0 depth=5" | iloc sc3db
     ```
The rest of the workflow remains the same, except the following steps:
 
1. we need to update the event xml with the additional picks,
2. insert the updated event into the seiscomp3 db,
3. use the iloc command on this updated event xml, and
4. optionally, insert the updated event back into seiscomp3 db.
