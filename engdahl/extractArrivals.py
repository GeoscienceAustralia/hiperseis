from ArrivalParser import ArrivalParser
from HDFLineParser import HDFLineParser
from itertools import groupby
import linecache
import obspy
from obspy import UTCDateTime
from obspy.core.event import Amplitude, Event, Magnitude, Origin, Pick,\
   OriginQuality, CreationInfo, WaveformStreamID

#def setPickData()

def setEventData(eventParser, arrivals):
   event = Event()
   creation_info = CreationInfo(author='niket_engdahl_data',
                                creation_time=UTCDateTime(),
                                agency_uri='NA',
                                agency_id='engdahl')
   event.creation_info = creation_info
   event.event_type = 'earthquake'

   origin = Origin()
   seconds = eventParser.sec.split('.')[0]
   mseconds = eventParser.sec.split('.')[1] + '0'
   origin.time = UTCDateTime(eventParser.iyr, eventParser.mon, eventParser.iday,
                             eventParser.ihr, eventParser.min, seconds, mseconds)
   origin.longitude = eventParser.glon
   origin.latitude = eventParser.glat

   #a field called ser in hdf file represents standard error in position
   #there could be a way of deriving longitude and latitude error values
   #from this ser value. but needs to be clarified
   #origin.longitude_errors = eventParser.
   #origin.latitude_errors = eventParser.

   origin.depth = eventParser.depth
   origin.depth_errors = eventParser.sedep

   origin.method_id = 'EHB'
   origin.earth_model_id = 'ak135'

   origQual = OriginQuality(associated_phase_count=len(arrivals),
                            used_phase_count=len(arrivals),
                            standard_error=eventParser.se,
                            azimuthal_gap=eventParser.openaz2)
   origin.quality = origQual

   origin.evaluation_mode = 'automatic'
   creation_info = CreationInfo(author='niket_engdahl_data',
                                creation_time=UTCDateTime(),
                                agency_uri='NA',
                                agency_id='engdahl')
   origin.creation_info = creation_info

   pickList = []
   for arrParser in arrivals:
      pick = Pick(time=UTCDateTime(arrParser.year, arrParser.month, arrParser.day, arrParser.hour,
                                   arrParser.minute, arrParser.second.split('.')[0], arrParser.second.split('.')[1]),
                  waveform_id=WaveformStreamID(station_code=arrParser.station))
      pickList.append(pick)
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
