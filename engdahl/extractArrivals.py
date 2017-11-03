#import urllib2
from ArrivalParser import ArrivalParser
from HDFLineParser import HDFLineParser
from itertools import groupby
import linecache
import obspy
from obspy import UTCDateTime
from obspy.core.event.header import PickOnset, PickPolarity
from obspy.core.event import Amplitude, Event, Magnitude, Origin, Pick,\
   OriginQuality, CreationInfo, WaveformStreamID, ResourceIdentifier, Arrival

'''
# This code is for associating stations with networks
# Commenting it out for now.
from obspy.clients.fdsn import Client

      client = Client(base_url='https://service.iris.edu/fdsnws/', debug=True)
      if not client.services.has_key('station') or not client.services.has_key('station'):
         print('The fdsn client of obspy was not initiated successfully. The network code will not be populated.')
      else:
         inv = client.get_stations(station=arrParser.station, level='network')
         networkCodes = [cod for cod in set([net.code for net in inv.networks]) if cod != 'SY']
         if len(networkCodes) > 1:
            print('Multiple network codes returned for the station: ' + arrParser.station +
                  '. Continuing with this arrival without network code assignment. ' + arrParser.__dict__)
         elif len(networkCodes) == 0:
            print('No network found for this station: ' + arrParser.station +
                  '. Continuing with this arrival without network code assignment. ' + arrParser.__dict__)
         else:
            arrParser.network = networkCodes[0]
'''

global eventCount
eventCount = 1
global originCount
originCount = 1
global pickCount
pickCount = 1

def setEventData(eventParser, arrivals, count):
   global originCount
   global eventCount
   global pickCount
   creation_info = CreationInfo(author='niket_engdahl_parser',
                                creation_time=UTCDateTime(),
                                agency_uri=ResourceIdentifier(id='smi:engdahl.ga.gov.au/ga-engdahl'),
                                agency_id='ga-engdahl')

#   magnitudeSurface = Magnitude(resource_id=ResourceIdentifier(id='smi:engdahl.ga.gov.au/origin/'+str(originCount)+'#netMag.Ms'),
#                         mag=eventParser.ms,
#                         magnitude_type='Ms',
#                         origin_id=ResourceIdentifier(id='smi:engdahl.ga.gov.au/origin/'+str(originCount)),
#                         azimuthal_gap=eventParser.openaz2,
#                         creation_info=creation_info)
   origin = Origin(resource_id=ResourceIdentifier(id='smi:engdahl.ga.gov.au/origin/'+str(originCount)),
                   time=UTCDateTime(int(str(2000 + int(eventParser.iyr))), int(eventParser.mon), int(eventParser.iday),
                                    int(eventParser.ihr), int(eventParser.min), int(eventParser.sec.split('.')[0]),
                                    int(eventParser.sec.split('.')[1] + '0')),
                   longitude=eventParser.glon,
                   latitude=eventParser.glat,
                   depth=float(eventParser.depth) * 1000, # engdahl files report kms, obspy expects m
                   depth_errors=eventParser.sedep,
                   method_id=ResourceIdentifier(id='EHB'),
                   earth_model_id=ResourceIdentifier(id='ak135'),
                   quality=OriginQuality(associated_phase_count=len(arrivals),
                            used_phase_count=len(arrivals),
                            standard_error=eventParser.se,
                            azimuthal_gap=eventParser.openaz2),
                   evaluation_mode='automatic',
                   creation_info=creation_info)

   magnitude = Magnitude(resource_id=ResourceIdentifier(id='smi:engdahl.ga.gov.au/origin/'+str(originCount)+'#netMag.Mb'),
                         mag=eventParser.mb,
                         magnitude_type='Mb',
                         origin_id=ResourceIdentifier(id='smi:engdahl.ga.gov.au/origin/'+str(originCount)),
                         azimuthal_gap=eventParser.openaz1,
                         creation_info=creation_info)

   originCount += 1

   pickList = []
   arrivalList = []
   pPhaseArrival = None
   for arrParser in arrivals:
      pickOnset = None
      pol = None

      if arrParser.phase and arrParser.phase.lower().startswith('p') and arrParser.year and arrParser.month and arrParser.day:
         pPhaseArrival = arrParser
      else:
         arrParser.year = pPhaseArrival.year
         arrParser.day = pPhaseArrival.day
         arrParser.month = pPhaseArrival.month
         arrParser.station = pPhaseArrival.station
         arrParser.delta = pPhaseArrival.delta
         arrParser.dtdd = pPhaseArrival.dtdd
         arrParser.backaz = pPhaseArrival.backaz
         arrParser.focalDip = pPhaseArrival.focalDip
         arrParser.angleAzimuth = pPhaseArrival.angleAzimuth

      if arrParser.phase1 == 'LR' or arrParser.phase2 == 'LR':
         continue

      if arrParser.phase1.startswith('i'):
         pickOnset = PickOnset.impulsive
         if arrParser.fm == '+':
            pol = PickPolarity.positive
         elif arrParser.fm == '-':
            pol = PickPolarity.negative
      elif arrParser.phase1.startswith('e'):
         pickOnset = PickOnset.emergent

      pick = Pick(resource_id=ResourceIdentifier(id='smi:engdahl.ga.gov.au/pick/'+str(pickCount)),
                  time=UTCDateTime(int(str(2000 + int(arrParser.year))), int(arrParser.month), int(arrParser.day),
                                   int(arrParser.hour), int(arrParser.minute), int(arrParser.second.split('.')[0]),
                                   int(arrParser.second.split('.')[1] + '0')),
                  waveform_id=WaveformStreamID(network_code='', station_code=arrParser.station, channel_code='BHZ'),
                  methodID=ResourceIdentifier('STA/LTA'),
                  backazimuth=arrParser.backaz if arrParser.backaz else None,
                  onset=pickOnset,
                  phase_hint=arrParser.phase,
                  polarity=pol,
                  evaluation_mode='automatic',
                  # TO-DO
                  comment='populate all the remaining fields here as key value',
                  creation_info=creation_info)
      if not arrParser.backaz:
          print "arrParser.backaz is empty. printing the arrParser for debugging"
      pickCount += 1
      pickList.append(pick)

      arrival = Arrival(pick_id=ResourceIdentifier(id='smi:engdahl.ga.gov.au/pick/'+str(pickCount-1)),
                        phase=arrParser.phase if arrParser.phase else None,
                        azimuth=arrParser.backaz if arrParser.backaz else None,
                        distance=arrParser.delta if arrParser.delta else None,
                        # if the * has some significance, it should be accounted for. ignoring for now.
                        time_residual=arrParser.residual.rstrip('*'),
                        time_weight=arrParser.wgt if arrParser.wgt else None,
                        backazimuth_weight=arrParser.wgt if arrParser.wgt else None)
      arrivalList.append(arrival)
      if not arrParser.wgt:
          print "arrParser.wgt is empty. printing the arrParser for debugging"
#          pprint.pprint(arrParser)

   origin.arrivals = arrivalList

   event = Event(resource_id=ResourceIdentifier(id='smi:engdahl.ga.gov.au/event/'+str(eventCount)),
                 creation_info=creation_info, event_type='earthquake')

   eventCount += 1

   event.picks = pickList
   event.origins = [origin, ]
   event.magnitudes = [magnitude, ]
   event.preferred_origin_id = origin.resource_id
   event.preferred_magnitude_id = magnitude.resource_id
   return event

def getArrivalSections(hdfFile, outFile):
   with open(outFile, 'r') as of:
      with open(hdfFile, 'r') as hf:
         count = 0
         for keyvaluetuple in groupby(of, lambda x: "PACIFIC EVENTS" in x.strip()):
            key, value = keyvaluetuple
            if not key:
               count = count + 1
               eventParser = HDFLineParser(linecache.getline(hdfFile, count).rstrip('\n'))

               listOfArrivals = []
               for groupline in value:
                  if len(groupline) == 133:
                     arrivalParser = ArrivalParser(groupline.rstrip('\n'))
                     listOfArrivals.append(arrivalParser)
                  elif len(groupline) == 115 and groupline.endswith('km'):
                     eventParser.seTime = groupline[groupline.substring('se h =   ') + 9 : groupline.substring(' sec     se lat')]
                     eventParser.seLat = groupline[groupline.substring('se lat =  ') + 10 : groupline.substring(' km     se lon')]
                     eventParser.seLon = groupline[groupline.substring('se lon =  ') + 10 : groupline.substring(' km     se depth')]

               ev = setEventData(eventParser, listOfArrivals, count)
               ev.write(ev.resource_id.id.split('/')[-1]+".xml", format='SC3ML')
               #remove the below if condition after testing, its only for debugging
               #if count > 2:
               #    break


def main():
   getArrivalSections('EHB.GA.HDF', 'EHB.GA.OUT')

if __name__ == "__main__":
   main()
