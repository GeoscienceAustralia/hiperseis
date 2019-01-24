#!/usr/env python
"""
Description:
    Reads waveforms (within a given time-range) from a Seiscomp3 server and
    dumps out ASDF files, along with a json file containing associated metadata
References:

CreationDate:   13/09/18
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     22/08/18   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import click
import os
import pyasdf

from obspy.core import UTCDateTime

from obspy.clients.fdsn import Client

import json

def make_ASDF_tag(tr, tag):
    # def make_ASDF_tag(ri, tag):
    data_name = "{net}.{sta}.{loc}.{cha}__{start}__{end}__{tag}".format(
        net=tr.stats.network,
        sta=tr.stats.station,
        loc=tr.stats.location,
        cha=tr.stats.channel,
        start=tr.stats.starttime.strftime("%Y-%m-%dT%H:%M:%S"),
        end=tr.stats.endtime.strftime("%Y-%m-%dT%H:%M:%S"),
        tag=tag)
    return data_name
# end func

def hasOverlap(stime1, etime1, stime2, etime2):
    result = 0
    
    if (etime1 is None): etime1 = UTCDateTime.now()
    if (etime2 is None): etime2 = UTCDateTime.now()
        
    if   (stime2 >= stime1 and stime2 <= etime1): result = 1
    elif (etime2 >= stime1 and etime2 <= etime1): result = 1
    elif (stime2 <= stime1 and etime2 >= etime1): result = 1
    
    return result
# end func

class DictToStr(dict):
    def __str__(self):
        return json.dumps(self)

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('host',
                type=str)
@click.argument('port',
                type=int)
@click.argument('output-path', required=True,
                type=click.Path(exists=True))
@click.argument('start-date', required=True,
                type=str)
@click.argument('end-date', required=True,
                type=str)
@click.option('--min-length-sec', type=int, default=None, help="Minimum length in seconds")                
@click.option('--merge-threshold', type=int, default=None, help="Merge traces if the number of traces fetched for an "
                                                                "interval exceeds this threshold")                
def process(host, port, output_path, start_date, end_date, min_length_sec, merge_threshold):
    """
    Example: python sc3toasdf.py 13.211.124.69 8081 /mnt/asdf_dump2/test/ 2018-01-01T00:00:00 2018-06-01T00:00:00 --merge-threshold 200
    """
    client = None
    try:
        client = Client('http://'+host+':%d'%port)

        if not client: raise Exception('Connection failed..')
    except Exception as e:
        print e
        print 'Failed to conncet to client. Aborting..'
        exit(0)
    # end try

    # Get inventory
    inv = client.get_stations()

    # Define increment, start and end times
    day = 24*3600
    stime = etime = None
    try:
        stime = UTCDateTime(start_date)
        etime = UTCDateTime(end_date)
    except Exception as e:
        print e
        print 'Incorrect start-date or end-date format. Aborting..'
    # end try

    # Start processing
    fn = os.path.join(output_path, '%d-%d.h5'%(stime.year, (etime-1).year))
    jsonfn = os.path.join(output_path, '%d-%d.json'%(stime.year, (etime-1).year))

    if(os.path.exists(fn)): os.remove(fn)
    ds = pyasdf.ASDFDataSet(fn, compression='gzip-3')
    
    jsonf = open(jsonfn, 'w+')
    jsonf.write('{')
    dictEntryCount = 0

    for network in inv.networks:
        if(not hasOverlap(stime, etime, network.start_date, network.end_date)): continue
        
        for station in network.stations:
            if(not hasOverlap(stime, etime, station.start_date, station.end_date)): continue
            
            ctime = stime 
            dcount = 0   
            trCount = 0
            while ctime < etime:                                    
                    
                cinv = None
                cwaveforms = None
                
                try:
                    cinv = client.get_stations(starttime=ctime, endtime=ctime+day, \
                            network=network.code, sta=station.code, loc='*', channel='*')

                    cwaveforms = client.get_waveforms(network.code, station.code, \
                            '*', '*', ctime, ctime+day)                    
                except Exception as e:
                    pass
                # end try
                
                if(cinv and cwaveforms):
                    
                    if(merge_threshold):
                        ntraces = len(cwaveforms)
                        if(ntraces > merge_threshold):
                            try:
                                cwaveforms = cwaveforms.merge(method=1, fill_value='interpolate')
                            except:
                                print 'Failed to merge traces. Moving along..'
                                ctime += day
                                dcount += 1
                                continue
                            # end try
                            print 'Merging stream with %d traces'%(ntraces)
                        # end if
                    # end if

                    for tr in cwaveforms:

                        if(tr.stats.npts == 0): continue
                        if(min_length_sec):
                            if(tr.stats.npts*tr.stats.delta < min_length_sec): continue
                        # end if

                        asdfTag = make_ASDF_tag(tr, "raw_recording").encode('ascii')
                        
                        try:
                            ds.add_waveforms(tr, tag='raw_recording')
                        except Exception as e:
                            print e
                            print 'Failed to append trace:'
                            print tr
                        # end try
                        
                        trdict = {"tr_starttime": tr.stats.starttime.timestamp,
                                  "tr_endtime": tr.stats.endtime.timestamp,
                                  "new_network": str(network.code),
                                  "new_station": str(station.code),
                                  "new_channel": str(tr.stats.channel),
                                  "new_location": str(tr.stats.location)}
                        if(dictEntryCount>0): jsonf.write(', ')
                        jsonf.write('"%s":%s'%(asdfTag, DictToStr(trdict)))
                        dictEntryCount += 1
                    # end for
                    
                    try:
                        ds.add_stationxml(cinv)
                    except Exception as e:
                        print e
                        print 'Failed to append inventory:'
                        print inv
                    # end try
                    
                    trCount += len(cwaveforms)
                # end if
                ctime += day
                dcount += 1
                
                if not (dcount % 100) : 
                    print 'Network: %s '%network.code + 'Station: %s '%station.code + \
                            str(ctime - day) + ' : added %d traces..'%trCount  
                    trCount = 0
                # end if
            # wend            
        # end for
    # end for
    
    print 'Closing asdf file..'
    del ds
    
    print 'Closing json database..'
    jsonf.write('}')
    jsonf.close()

    print 'Done.'
# end func

if (__name__ == '__main__'):
    process()
# end if
