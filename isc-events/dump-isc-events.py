import urllib2
import httplib
import datetime
import obspy
from obspy import read_events
from StringIO import StringIO
import os

baseurl = 'http://www.isc.ac.uk/cgi-bin/web-db-v4?request=REVIEWED&out_format=QuakeML&bot_lat=&top_lat=&left_lon=&right_lon=&searchshape=CIRC&ctr_lat=-22.5&ctr_lon=127.5&radius=140&max_dist_units=deg&srn=&grn=&start_year=%s&start_month=%s&start_day=%s&start_time=%s&end_year=%s&end_month=%s&end_day=%s&end_time=%s&min_dep=&max_dep=&min_mag=3&max_mag=&req_mag_type=&req_mag_agcy=&min_def=&max_def=&prime_only=on&include_phases=on&include_magnitudes=on'

initialstart = datetime.datetime(2009, 4, 1, 0, 0, 0)
finalend = datetime.datetime(2014, 10, 1, 0, 0, 0)
tdelta = datetime.timedelta(hours=1, minutes=0, seconds=0)

def start_dump(start, end, runFailed=False):
   while (end < finalend):
      try:
         targetDir = os.getcwd()+'/'+str(start.year)+'/'+str(start.month)+'/'+str(start.day)
         if not os.path.exists(targetDir):
            os.makedirs(targetDir)

         if runFailed and os.path.isfile(targetDir+'/isc-event-'+start.strftime("%Y-%m-%d-%H:%M:%S")+'.xml'):
            start = end
            end = end + tdelta
            continue

         response = urllib2.urlopen(baseurl % (str(start.year), str(start.month), str(start.day),
                                               str(start.hour)+'%3A'+str(start.minute)+'%3A'+str(start.second),
                                               str(end.year), str(end.month), str(end.day),
                                               str(end.hour)+'%3A'+str(end.minute)+'%3A'+str(end.second)))
         qmlResponse = response.read()
         qmlResponse = qmlResponse.replace('<evaluationMode>a</evaluationMode>', '<evaluationMode>automatic</evaluationMode>')
         qmlResponse = qmlResponse.replace('<evaluationMode>m</evaluationMode>', '<evaluationMode>manual</evaluationMode>')

         # dump the quakeml for debugging purposes
         outf = open('/tmp/1.xml', 'w')
         outf.write(qmlResponse)
         outf.close()

         # handle empty and malformed quakeml
         if not '<event ' in qmlResponse:
            print('The returned quakeml for start='+str(start)+' and end='+str(end)+' returned no events')
            start = end
            end = end + tdelta
            continue

         cat = read_events(StringIO(qmlResponse))

         cat.write(targetDir+'/isc-event-'+start.strftime("%Y-%m-%d-%H:%M:%S")+'.xml', format='SC3ML')

         start = end
         end = end + tdelta
      except urllib2.HTTPError, e:
         print('Gettting events for start='+str(start)+' and end='+str(end)+' did not succeed.')
         print e
      except urllib2.URLError, e:
         print('Gettting events for start='+str(start)+' and end='+str(end)+' did not succeed.')
         print e

start_dump(initialstart, initialstart+tdelta)

# run the failed attempts once
start_dump(initialstart, initialstart+tdelta, runFailed=True)

