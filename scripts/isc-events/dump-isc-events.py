import urllib2
import httplib
import datetime
import obspy
from obspy import read_events
from StringIO import StringIO
import os

baseurl = 'http://www.isc.ac.uk/cgi-bin/web-db-v4?request=REVIEWED&out_format=CATQuakeML&bot_lat=&top_lat=&left_lon=&right_lon=&searchshape=CIRC&ctr_lat=-22.5&ctr_lon=127.5&radius=140&max_dist_units=deg&srn=&grn=&start_year=%s&start_month=%s&start_day=%s&start_time=%s&end_year=%s&end_month=%s&end_day=%s&end_time=%s&min_dep=&max_dep=&min_mag=4&max_mag=&req_mag_type=&req_mag_agcy=&prime_only=on&include_links=on'

initialstart = datetime.datetime(1993, 5, 1, 0, 0, 0)
finalend = datetime.datetime(1994, 1, 1, 0, 0, 0)
tdelta = datetime.timedelta(hours=1, minutes=0, seconds=0)

def start_dump(start, end, runFailed=False):
   temp_start = start
   temp_end = start + tdelta
   while (temp_end < end):
      try:
         targetDir = os.getcwd()+'/'+str(temp_start.year)+'/'+str(temp_start.month)+'/'+str(temp_start.day)
         if not os.path.exists(targetDir):
            os.makedirs(targetDir)

         if runFailed and os.path.isfile(targetDir+'/isc-event-'+temp_start.strftime("%Y-%m-%d-%H:%M:%S")+'.xml'):
            temp_start = temp_end
            temp_end = temp_end + tdelta
            continue

         response = urllib2.urlopen(baseurl % (str(temp_start.year), str(temp_start.month), str(temp_start.day),
                                               str(temp_start.hour)+'%3A'+str(temp_start.minute)+'%3A'+str(temp_start.second),
                                               str(temp_end.year), str(temp_end.month), str(temp_end.day),
                                               str(temp_end.hour)+'%3A'+str(temp_end.minute)+'%3A'+str(temp_end.second)))
         qmlResponse = response.read()
         qmlResponse = qmlResponse.replace('<evaluationMode>a</evaluationMode>', '<evaluationMode>automatic</evaluationMode>')
         qmlResponse = qmlResponse.replace('<evaluationMode>m</evaluationMode>', '<evaluationMode>manual</evaluationMode>')

         # dump the quakeml for debugging purposes
         outf = open('/tmp/1.xml', 'w')
         outf.write(qmlResponse)
         outf.close()

         # handle empty and malformed quakeml
         if not '<event ' in qmlResponse:
            if 'No events were found' in qmlResponse:
               print('There are no events recoreded at ISC for start='+str(temp_start)+' and end='+str(temp_end))
            else:
               print('The returned quakeml for start='+str(temp_start)+' and end='+str(temp_end)+' returned no events')
            temp_start = temp_end
            temp_end = temp_end + tdelta
            continue

         cat = read_events(StringIO(qmlResponse))

         cat.write(targetDir+'/isc-event-'+temp_start.strftime("%Y-%m-%d-%H:%M:%S")+'.xml', format='SC3ML')

         temp_start = temp_end
         temp_end = temp_end + tdelta
      except urllib2.HTTPError, e:
         print('Gettting events for start='+str(temp_start)+' and end='+str(temp_end)+' did not succeed.')
         print e
      except urllib2.URLError, e:
         print('Gettting events for start='+str(temp_start)+' and end='+str(temp_end)+' did not succeed.')
         print e

start_dump(initialstart, finalend)

# run the failed attempts once
# start_dump(initialstart, initialstart+tdelta, runFailed=True)

