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

1. we need to update the original event xml with the additional picks,
2. insert the updated event into the seiscomp3 db,
3. use the iloc command on this updated event xml, and
4. optionally, insert the updated event back into seiscomp3 db.


## How to run iLoc using SQL with different approaches

1. For one event:
```
echo "`mysql -u alexei -sN -e"select ep.publicID, 'update_db=0 do_gridsearch=0 use_RSTT_PnSn=1 use_RSTT_PgLg=1' from Event e,PublicObject ep where e._oid=ep._oid and ep.publicID = 'quakeml:ga.ga.gov.au/event/00306594'" | sed 's/\t/ /g'`" | iloc sc3db > my.log
```

2. For a subset of events:
```
echo "`mysql -u alexei -sN -e"select ep.publicID, 'update_db=0 do_gridsearch=0 use_RSTT_PnSn=1 use_RSTT_PgLg=1' from Event e,PublicObject ep where e._oid=ep._oid and ep.publicID like '%306%'" | sed 's/\t/ /g'`" | iloc sc3db > subset.log
```

3. For all events:
```
echo "`mysql -u alexei -sN -e"select ep.publicID, 'update_db=0 do_gridsearch=0 use_RSTT_PnSn=1 use_RSTT_PgLg=1' from Event e,PublicObject ep where e._oid=ep._oid" | sed 's/\t/ /g’`" | iloc sc3db > all.log
```

## The combination of ak135 and RSTT locations
In the configuration file use 1 or 0 , where 1 is ak135 and 0 is RSTT
use_RSTT_PnSn=1
use_RSTT_PgLg=1


edit config.txt and change out_agency

out_agency=GAak135

Then run the script:

iloc sc3db < instructions.ak135.txt > ak135.log

## Useful SQL queries to manage event database

Set the verbose level to 3 to get more detailed info about iloc's inner workings.
In particular it lists the sql queries iloc makes.

1. Get_sc3_data

```
echo "quakeml:ga.ga.gov.au/event/00377167 update_db=0 do_gridsearch=0 verbose=3" | iloc sc3db > z.log
```

2. Get_sc3_event, get evid, orid, depth, event type and preferred magnitude id from Event and Origin:
```
select e._oid, o._oid, o.depth_value,
coalesce(e.typeCertainty, 'k'), coalesce(e.type, 'e'),
coalesce(e.preferredMagnitudeID, 'X')
from Event e, Origin o, alexeidb.PublicObject ep, alexeidb.PublicObject op
where ep.publicID='quakeml:ga.ga.gov.au/event/00377167' and ep._oid=e._oid
and op.publicID=preferredOriginID and op._oid=o._oid;
```

3. Get number of origins belonging to an event

```
select count(*) from OriginReference r, alexeidb.PublicObject ep
where ep.publicID='quakeml:ga.ga.gov.au/event/00377167'
and r._parent_oid=ep._oid;
        3 hypocentres
```

4. Get number of associated phase, readings and stations belonging to the
preferred origin

```
select count(r._oid), count(distinct p.waveformID_stationCode, p.waveformID_networkCode)
from Arrival r, Pick p, alexeidb.PublicObject rp
where r._parent_oid=1256734 and rp.publicID=r.pickID and rp._oid=p._oid
        7 phases, 4 readings, 4 stations
```

5. Get number of agencies that reported the event

```
select distinct p.creationInfo_agencyID
from Arrival r, Pick p, alexeidb.PublicObject rp
where r._parent_oid=1256734 and rp.publicID=r.pickID and rp._oid=p._oid
        1 agencies

    evid=689051 prime=1256734 numhyps=3 numphas=7
```

6. Get preferred origin

```
select coalesce(creationInfo_agencyID, 'BUD'), _oid,
    time_value, coalesce(time_value_ms, 0),
    latitude_value, longitude_value,
    coalesce(depth_value, 0), coalesce(depth_uncertainty, 0),
    coalesce(epicenterFixed, 0), coalesce(timeFixed, 0),
    coalesce(quality_associatedStationCount, 0), coalesce(quality_usedStationCount, 0),
    coalesce(quality_associatedPhaseCount, 0), coalesce(quality_usedPhaseCount, 0),
    coalesce(quality_minimumDistance, 0), coalesce(quality_maximumDistance, 0),
    coalesce(quality_azimuthalGap, 360), coalesce(quality_secondaryAzimuthalGap, 360),
    coalesce(uncertainty_minHorizontalUncertainty, 0),
    coalesce(uncertainty_maxHorizontalUncertainty, 0),
    coalesce(uncertainty_azimuthMaxHorizontalUncertainty, 0),
    coalesce(time_uncertainty, 0), coalesce(quality_standardError, 0),
    coalesce(creationInfo_author, 'BUD')
from Origin where _oid=1256734

        i=0 evid=689051 hypid=1256734 agency=RSTT
```

7. Get the rest of origins
```
select coalesce(o.creationInfo_agencyID, 'BUD'), o._oid,
    o.time_value, coalesce(o.time_value_ms, 0),
    o.latitude_value, o.longitude_value,
    coalesce(o.depth_value, 0), coalesce(o.depth_uncertainty, 0),
    coalesce(o.epicenterFixed, 0), coalesce(o.timeFixed, 0),
    coalesce(o.quality_associatedStationCount, 0), coalesce(o.quality_usedStationCount, 0),
    coalesce(o.quality_associatedPhaseCount, 0), coalesce(o.quality_usedPhaseCount, 0),
    coalesce(o.quality_minimumDistance, 0), coalesce(o.quality_maximumDistance, 0),
    coalesce(o.quality_azimuthalGap, 360), coalesce(o.quality_secondaryAzimuthalGap, 360),
    coalesce(o.uncertainty_confidenceEllipsoid_semiMinorAxisLength, 0),
    coalesce(o.uncertainty_confidenceEllipsoid_semiMajorAxisLength, 0),
    coalesce(o.uncertainty_confidenceEllipsoid_majorAxisAzimuth, 0),
    coalesce(o.time_uncertainty, 0), coalesce(o.quality_standardError, 0),
    coalesce(o.creationInfo_author, 'BUD')
from Origin o, OriginReference r, alexeidb.PublicObject rp
where r._parent_oid=689051 and rp.publicID=r.OriginID and rp._oid=o._oid
   and o._oid <> 1256734
   and coalesce(o.depth_value,0) between 0 and 700
order by o.creationInfo_agencyID, o.creationInfo_creationTime desc
        i=1 evid=689051 hypid=689037 agency=GA
        i=2 evid=689051 hypid=1256618 agency=ak135
```

8. Get network magnitudes belonging to the various origins

```
select coalesce(magnitude_value, -9), coalesce(nullif(type, ''), 'Q'),
    coalesce(magnitude_uncertainty, -1), coalesce(stationCount, -1),
    creationInfo_agencyID
from Magnitude m, alexeidb.PublicObject p
where m._parent_oid=1256734 and m._oid=p._oid
  and p.publicID !='quakeml:ga.ga.gov.au/magnitude/837217'
  and m.creationInfo_agencyID !='RSTT'
order by type
```
or
```
select coalesce(magnitude_value, -9), coalesce(nullif(type, ''), 'Q'),
    coalesce(magnitude_uncertainty, -1), coalesce(stationCount, -1),
    creationInfo_agencyID
from Magnitude m, alexeidb.PublicObject p
where m._parent_oid=689037 and m._oid=p._oid
  and p.publicID !='quakeml:ga.ga.gov.au/magnitude/837217'
  and m.creationInfo_agencyID !='RSTT'
order by type
```

or

```
select coalesce(magnitude_value, -9), coalesce(nullif(type, ''), 'Q'),
    coalesce(magnitude_uncertainty, -1), coalesce(stationCount, -1),
    creationInfo_agencyID
from Magnitude m, alexeidb.PublicObject p
where m._parent_oid=1256618 and m._oid=p._oid
  and p.publicID !='quakeml:ga.ga.gov.au/magnitude/837217'
  and m.creationInfo_agencyID !='RSTT'
order by type
```

9. Get preferred network magnitude
```
select coalesce(magnitude_value, -9), coalesce(nullif(type, ''), 'Q'),
    coalesce(magnitude_uncertainty, -1), coalesce(stationCount, -1),
    creationInfo_agencyID
from Magnitude m, alexeidb.PublicObject p
where m._oid=p._oid
  and p.publicID='quakeml:ga.ga.gov.au/magnitude/837217'
  and m.creationInfo_agencyID !='RSTT'
order by type
```


10. Get arrival data. Note the check for the valid operation time for stations!
This is why you may get an error message from iloc complaining about different
number of phases what it got from query #3

```
select r._oid, coalesce(p.creationInfo_agencyID, 'BUD'),
    p.time_value, coalesce(p.time_value_ms, 0),
    coalesce(r.phase_code, ''), coalesce(p.phaseHint_code, r.phase_code, ''),
    p.waveformID_stationCode, coalesce(p.waveformID_networkCode, '--'),
    coalesce(nullif(p.waveformID_locationCode, ''), '--'),
    coalesce(p.waveformID_channelCode, '???'),
    coalesce(p.horizontalSlowness_value, 0), coalesce(p.backazimuth_value, 0),
    coalesce(r.distance, 0), coalesce(p.polarity, ' '), coalesce(p.onset, ' '),
    s.latitude, s.longitude, coalesce(s.elevation, 0),
    r.pickID, coalesce(p.creationInfo_author, 'BUD')
from Arrival r, Pick p, alexeidb.PublicObject rp, alexeidb.Network n, alexeidb.Station s
where r._parent_oid=1256734 and rp.publicID=r.pickID and rp._oid=p._oid
  and p.waveformID_networkCode=n.code and p.waveformID_stationCode=s.code
  and s._parent_oid=n._oid
  and p.time_value between s.start and coalesce(s.end, '3001-01-01')
order by r._oid

```

11. Get amplitude data

```
select r._oid, a._oid, a.type,
    coalesce(a.amplitude_value, 0), coalesce(a.period_value, 0), coalesce(a.snr, 0),
    a.pickID, coalesce(a.waveformID_channelCode, '???')
from Amplitude a, Arrival r
where r._parent_oid=1256734 and r.pickID=a.pickid
order by r._oid
```

12.  Joins Event table with OriginReference and Origin tables to get input (GA)
and output (ak135, RSTT) locations

```
select ep.publicID,o.creationInfo_agencyID,
    o.time_value,o.time_value_ms,o.time_uncertainty,
    o.latitude_value,o.longitude_value,
    o.uncertainty_minHorizontalUncertainty,o.uncertainty_maxHorizontalUncertainty,
    o.uncertainty_azimuthMaxHorizontalUncertainty,
    o.depth_value,o.depth_uncertainty,o.depthType,
    o.quality_associatedPhaseCount,o.quality_usedPhaseCount,o.quality_usedStationCount,
    o.quality_standardError,o.quality_azimuthalGap,o.quality_secondaryAzimuthalGap
from Event e,PublicObject ep,Origin o,OriginReference r,PublicObject op
where e._oid=ep._oid and r._parent_oid=e._oid
  and r.originID=op.publicID and o._oid=op._oid
  and ep.publicID in ('quakeml:ga.ga.gov.au/event/00305736',
                      'quakeml:ga.ga.gov.au/event/00306594',
                      'quakeml:ga.ga.gov.au/event/00309646',
                      'quakeml:ga.ga.gov.au/event/00345974',
                      'quakeml:ga.ga.gov.au/event/00348254',
                      'quakeml:ga.ga.gov.au/event/00358495',
                      'quakeml:ga.ga.gov.au/event/00363318',
                      'quakeml:ga.ga.gov.au/event/00370363',
                      'quakeml:ga.ga.gov.au/event/00374747',
                      'quakeml:ga.ga.gov.au/event/00377167');

```


13. The query below list the stations and their time residuals at very close, local distances.  The result is ordered by station and distance.

```
select p.waveformID_stationCode, r.distance, r.azimuth,
    coalesce(r.phase_code, ''), r.timeResidual, ep.publicID
from Arrival r, Pick p, alexeidb.PublicObject rp,
     Event e, Origin o, alexeidb.PublicObject ep, alexeidb.PublicObject op
where ep._oid=e._oid and o.creationInfo_agencyID='GA'
  and op.publicID=preferredOriginID and op._oid=o._oid
  and r._parent_oid=o._oid and rp.publicID=r.pickID and rp._oid=p._oid
  and r.distance between 0 and 0.5
order by p.waveformID_stationCode,r.distance;
```

13. Delete events with no arrivals

```

create table tmp as select _oid from Event where _oid in (select e._oid from Event e,Origin o,PublicObject op where e.preferredOriginID=op.publicID and op._oid=o._oid and o._oid not in (select distinct _parent_oid from Arrival)); Query OK, 262 rows affected (0.16 sec)
Records: 262  Duplicates: 0  Warnings: 0

delete from Event where _oid in (select _oid from tmp); Query OK, 262 rows affected (1.18 sec)

commit;

```

14. Preferred origins with no arrivals

```

select ep.publicID,e._oid,e.preferredOriginID from Event e,PublicObject ep where e._oid=ep._oid and e.preferredOriginID not in (select OriginID from OriginReference) order by e._oid;

select ep.publicID,e._oid,e.preferredOriginID,o._oid
from Event e,PublicObject ep,Origin o,PublicObject op where e._oid=ep._oid and e.preferredOriginID=op.publicID and op._oid=o._oid and o._oid not in (select distinct _parent_oid from Arrival) order by o._oid;

```

15. Empty agency detection

```
mysql> select distinct creationInfo_agencyID from Origin;

mysql> select distinct creationInfo_agencyID from Event;

mysql> select distinct creationInfo_agencyID from Pick;

mysql> select distinct creationInfo_agencyID from Arrival;

mysql> select distinct creationInfo_agencyID from Magnitude;

```

16. Remove spaces from author field

```

select distinct creationInfo_author from Origin;


update Origin set creationInfo_author = replace(creationInfo_author, ' ', '’);
commit;


```

17. Check for missing stations in the inventory when arrival contains their record

```

select distinct waveformID_networkCode,waveformID_stationCode from Pick
where waveformID_stationCode not in (select code from Station) order by waveformID_stationCode;

```


## Manipulation of  RSTT

1. Extract nodes from RSTT

```

java -cp "/home/istvan/SLBM_Root.3.0.5.Linux/RSTT_NOGS/RSTT_NOGS.jar:/home/istvan/SLBM_Root.3.0.5.Linux/lib/slbmjni.jar" RSTT.GetNodes /home/istvan/SLBM_Root.3.0.5.Linux/models/rstt201404um.geotess -46 -13 116 156

```

2. Run RSTT_NOGS

```

(in the cloud machine the java part of RSTT is not yet installed)

java -cp "/home/centos/SLBM_Root.3.0.5.Linux/RSTT_NOGS/RSTT_NOGS.jar:/home/centos/SLBM_Root.3.0.5.Linux/lib/slbmjni.jar" RSTT.GetNodes
Usage: GetNodes <full path to RSTT model> latmin latmax lonmin lonmax
         Writes information for nodes in specifed bounds to RSTTnodes.out
         RSTTnodes.out has a node line 'node nodeID activeNodeID lat lon'.  note: nodeID applies to the global model and activeNodeID is a numbering system for the subset of selected nodes
                             crustal layer lines 'layer# depth (km) Vp (km/s) Vs (km/s)'
                             mantle layer lines 'layer# depth (km) Vp (km/s) Vs (km/s) Vp_gradient (km/s/km) Vs_gradient (km/s/km)'

                   OR

Usage: GetNodes <full path to RSTT model> nodeID    // prints the information for an individual node to the screen

                   OR

Usage: GetNodes <full path to RSTT model> <nodeIDfile>
         nodeIDfile contains a list of nodeIDs.  The list of nodeIDs is taken from the first column of the file.
         Output is to a file named <nodeIDfile>.nodeInfo.

                   OR

Usage: GetNodes <full path to RSTT model> <nodeIDfileExclude> <anyString>
         Outputs all nodeIDs except the ones listed in nodeIDfileExclude (nodeIDs in first column). Output file is RSTT_nodesMinus<nodeIDfileExclude>.


java -cp "/home/centos/SLBM_Root.3.0.5.Linux/RSTT_NOGS/RSTT_NOGS.jar:/home/centos/SLBM_Root.3.0.5.Linux/lib/slbmjni.jar" RSTT.SetNodes
Usage: SetNodes <input model path> <node update file> <output model path>
         <input model path> is the full path to the file containing the input model.
         <output model path> is the full path to the file containing the output model. The output model file is created if it doesn't exist.
         <node update file> (Use GetNodes for an example) is an ascii file containing the node update information. Multiple nodes can be updated.
            The first line of a node entry is 'node nodeID activeNodeID lat lon'. nodeID must be taken from the input model. activeNodeID is optional and is not currently used.
            9 lines specifying layer depths and velocites follow each node line.
            Lines 1-8 have entries 'layer# depth (km) Vp (km/s) Vs (km/s)' layer# starts with 0.
            Line 9 is similar to lines 1-8 with  additional entries that specify mantle velocity gradient for Vp and Vs. Units for gradient are km/s/km or 1/s.

```























