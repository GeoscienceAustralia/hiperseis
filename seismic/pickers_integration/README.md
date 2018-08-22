## How to use the pick extraction code

1. Create an input folder
2. Copy over the input event catalog to the input folder
3. Create an output folder
4. Invoke the command below:

```
    $ python3 extract_picks.py --help
    Usage: extract_picks.py [OPTIONS] SC3_HOST YEAR
    
      This script accomplishes the task of picking the p-phase and s-phase
      arrival times for given set of events and waveform data sets that have
      been ingested as miniseed data into a Seiscomp3 (SC3 for short) server
      SDS archive.
    
    Arguments:
    
      SC3_HOST    : the SC3 FDSN server where the permanent station waveform
      data is hosted
    
      YEAR        : the year to which the input event catalog corresponds
    
    Options:
      --input_folder TEXT   The input folder where ISC event catalog xml input
                            files are kept.
      --output_folder TEXT  The output folder where generated xml files are
                            generated.
      --parallel            Needs to be run in conjunction with MPI.(DO NOT USE.
                            NOT TESTED.)
      --help                Show this message and exit.
```

## The algorithm for pick extraction is as below

```
For a given year or a given time period:
	a. Query the SC3 server to get the stations that participated in seismic monitoring ans store in a list
	b. Parse the events corresponding to this time period into a list of obspy event objects
	c. iterate through the list of events; for each event:
		i.   extract the origin time of the event
		ii.  get the theoretical travel time for p-phase and s-phase arrivals for this station event pair
		ii.  for each station in the station list pick the p-phase arrival:
			A. query SC3 server to retrieve the seismogram for this station event pair within a time window around origin time for p-phase picking
			B. resample the seismogram to 20Hz
			C. create a list of frequency bands for p-phase that we would like to filter the seismogram with
			D. pick the p-phase and add to a list of p-picks if successful. For each frequency band:
				1. trim the seismogram to the relevant time window
				2. Apply STA/LTA picker to the seismogram and store triggers into a list
					- detrend and taper the seismogram
					- filter the seismogram with this frequency band and zerophase as True
					- calculate characteristic function for recursive STA/LTA
					- return STA/LTA if a trigger is detected
				3. Choose the best trigger
				4. For the best trigger, calculate AIC, AIC first derivative and AIC global minimum
				5. (debugging) call some plotting routines for validation
				6. If AIC is successful, calculate SNR
			E. create obspy pick object for each station that was successfully picked
		iv.  for each station that was successfully picked for p-phase for this event, pick the s-phase and add to a list of s-picks if successful
			Perform steps A, B, C, D and E for s-phase picking. The time window is broader and the LTA time window is bigger for s-phase
		v.   create an output obspy event object with the successful pick objects added and write to an xml file
```
The output events generated from this will be in the output folder.
