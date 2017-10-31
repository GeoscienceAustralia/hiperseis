from ArrivalParser import ArrivalParser
from HDFLineParser import HDFLineParser
from itertools import groupby
import linecache
import obspy
from obspy import UTCDateTime
from obspy.core.event.header import PickOnset, PickPolarity
from obspy.core.event import Amplitude, Event, Magnitude, Origin, Pick,\
   OriginQuality, CreationInfo, WaveformStreamID, ResourceIdentifier, Arrival

def setEventData(eventParser, arrivals):
   event = Event()
   creation_info = CreationInfo(author='niket_engdahl_data',
                                creation_time=UTCDateTime(),
                                agency_uri='NA',
                                agency_id='engdahl')
   event.creation_info = creation_info
   event.event_type = 'earthquake'

   origin = Origin(time=UTCDateTime(eventParser.iyr, eventParser.mon, eventParser.iday,
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

   pickList = []
   pPhaseArrival = None
   for arrParser in arrivals:
      pickOnset = None
      pol = None
      if arrParser.phase1.startswith('i'):
         pickOnset = PickOnset.impulsive
         pol = arrParser.fm
      elif arrParser.phase1.startswith('e'):
         pickOnset = PickOnset.emergent
      pick = Pick(time=UTCDateTime(arrParser.year, arrParser.month, arrParser.day, arrParser.hour,
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
                  creation_info = creation_info)
      pickList.append(pick)
      arrival = Arrival()
   event.picks = pickList

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
#                     ev = setPickData(arrivalParser, ev)
#                     print arrivalParser.__dict__
                  elif len(groupline) == 115 and groupline.endswith('km'):
                     eventParser.seTime = groupline[groupline.substring('se h =   ') + 9 : groupline.substring(' sec     se lat')]
                     eventParser.seLat = groupline[groupline.substring('se lat =  ') + 10 : groupline.substring(' km     se lon')]
                     eventParser.seLon = groupline[groupline.substring('se lon =  ') + 10 : groupline.substring(' km     se depth')]

               ev = setEventData(eventParser, listOfArrivals)
               #ev.write(filename, format='SC3ML')


def main():
   getArrivalSections('EHB.GA.HDF', 'EHB.GA.OUT')

if __name__ == "__main__":
   main()
