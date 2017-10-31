from ArrivalParser import ArrivalParser
from HDFLineParser import HDFLineParser
from itertools import groupby
import linecache
import obspy
from obspy import UTCDateTime
from obspy.core.event.header import PickOnset, PickPolarity
from obspy.core.event import Amplitude, Event, Magnitude, Origin, Pick,\
   OriginQuality, CreationInfo, WaveformStreamID, ResourceIdentifier, Arrival

global eventCount
eventCount = 1
global originCount
originCount = 1
global pickCount
pickCount = 1

def setEventData(eventParser, arrivals, count):
   creation_info = CreationInfo(author='niket_engdahl_data',
                                creation_time=UTCDateTime(),
                                agency_uri='NA',
                                agency_id='engdahl')
   origin = Origin(resource_id=ResourceIdentifier('engdahl:ga.gov.au/event/'+originCount),
                   time=UTCDateTime(eventParser.iyr, eventParser.mon, eventParser.iday,
                                    eventParser.ihr, eventParser.min, eventParser.sec.split('.')[0],
                                    eventParser.sec.split('.')[1] + '0'),
                   longitude=eventParser.glon,
                   latitude=eventParser.glat,
                   depth=eventParser.depth,
                   depth_errors=eventParser.sedep,
                   method_id='EHB',
                   earth_model_id='ak135',
                   quality=OriginQuality(associated_phase_count=len(arrivals),
                            used_phase_count=len(arrivals),
                            standard_error=eventParser.se,
                            azimuthal_gap=eventParser.openaz2),
                   evaluation_mode='automatic',
                   creation_info=creation_info)
   originCount += 1

   pickList = []
   arrivalList = []
   pPhaseArrival = None
   for arrParser in arrivals:
      pickOnset = None
      pol = None
      if arrParser.phase1.startswith('i'):
         pickOnset = PickOnset.impulsive
         pol = arrParser.fm
      elif arrParser.phase1.startswith('e'):
         pickOnset = PickOnset.emergent
      pick = Pick(resource_id=ResourceIdentifier('engdahl:ga.gov.au/pick/'+pickCount),
                  time=UTCDateTime(arrParser.year, arrParser.month, arrParser.day, arrParser.hour,
                                   arrParser.minute, arrParser.second.split('.')[0],
                                   arrParser.second.split('.')[1] + '0'),
                  waveform_id=WaveformStreamID(station_code=arrParser.station),
                  methodID=ResourceIdentifier('STA/LTA'),
                  backazimuth=arrParser.backaz,
                  onset=pickOnset,
                  phase_hint=arrParser.phase,
                  polarity=pol,
                  evaluation_mode='automatic',
                  # TO-DO
                  comment='populate all the remaining fields here as key value',
                  creation_info=creation_info)
      pickCount += 1
      pickList.append(pick)

      arrival = Arrival(pick_id=(pickCount-1),
                        phase=arrParser.phase[0],
                        azimuth=arrParser.backaz,
                        distance=arrParser.delta,
                        time_residual=arrParser.residual,
                        time_weight=arrParser.wgt,
                        backazimuth_weight=arrParser.wgt)
      arrivalList.append(arrival)

   origin.arrivals = arrivalList

   event = Event(resource_id=ResourceIdentifier('engdahl:ga.gov.au/event/'+eventCount),
                 creation_info=creation_info,
                 event_type='earthquake')

   eventCount += 1

   event.picks = pickList
   event.origins = [origin, ]
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

#               print eventParser.__dict__
               listOfArrivals = []
               for groupline in value:
                  if len(groupline) == 133:
                     arrivalParser = ArrivalParser(groupline.rstrip('\n'))
                     listOfArrivals.append(arrivalParser)
#                     print arrivalParser.__dict__
                  elif len(groupline) == 115 and groupline.endswith('km'):
                     eventParser.seTime = groupline[groupline.substring('se h =   ') + 9 : groupline.substring(' sec     se lat')]
                     eventParser.seLat = groupline[groupline.substring('se lat =  ') + 10 : groupline.substring(' km     se lon')]
                     eventParser.seLon = groupline[groupline.substring('se lon =  ') + 10 : groupline.substring(' km     se depth')]

               ev = setEventData(eventParser, listOfArrivals, count)
               ev.write(ev.resource_id.split('/')[-1]+".xml", format='SC3ML')
               #remove the below if condition after testing, its only for debugging
               if count > 5:
                   break


def main():
   getArrivalSections('EHB.GA.HDF', 'EHB.GA.OUT')

if __name__ == "__main__":
   main()
