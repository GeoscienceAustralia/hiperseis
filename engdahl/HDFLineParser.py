import struct
from obspy.core.event import Event, Pick, Amplitude, Magnitude, Origin

'''
Class that parses ascii data from an hdf format file
and the .out file to generate an object. The hdf format
is described below:
                            *HDF FILE FORMAT
*HDF file read statement
       read(1,100) 
     1 ahyp,isol,iseq,iyr,mon,iday,ihr,min,sec,ad,glat,glon,
     2 depth,iscdep,mb,ms,mw,ntot,ntel,ndep,igreg,se,ser,sedep,
     3 rstadel,openaz1,openaz2,az1,flen1,az2,flen2,avh,ievt
 100  format(a1,a3,a2,i2,2i3,1x,2i3,f6.2,a1,2f8.3,2f6.1,3f4.1,
     1 4i4,3f8.2,3f6.1,i4,f4.1,i4,f4.1,f5.1,i10)

Variable definitions:

ahyp         secondary teleseismic azimuth gap       a1
             blank = <= 180 deg
             Z     = >  180 deg
                    or
             primary teleseismic azimuth gap
             A     = <  180 deg
             B     = <  210 deg and  > 180 deg
             C     = <  240 deg and  > 210 deg
             D     = <  270 deg and  > 240 deg
             F     = >  270 deg
isol         solution type                           a3
             HEQ = origin time & hypocenter fixed
             DEQ = depth free (sedep <= 15 km)
             WEQ = depth fixed at waveform depth
             BEQ = depth fixed at BB depth
             FEQ = depth fixed by Engdahl
             LEQ = depth fixed by program (sedep > 15 km)
             XEQ = poor solution (ser > 35 km)
iseq         other info (b=blank)                    a2
             bx = earthquake location known to x km
             bc = only regional event data available
             Xb = explosion/cavity collapse
             Xx = explosion location known to x km
             Mb = CMT solution available
             Mx = CMT location known to x km
             Md = depth reviewed and accepted
             Mn = depth unreviewed but provisionally accepted 
                  based on => 5 depth phases and/or station(s) 
                  at distance(s) less than the focal depth
             Mr = depth under review
             Mf = depth set to regional depth estimate
             Mh = depth set to Harvard CMT depth
             Mx = depth reviewed but not accepted
             bd = depth reviewed and accepted
             bn = depth unreviewed, but provisionally accepted 
                  based on => 5 depth phases and/or station(s) 
                  at distance(s) less than the focal depth
             br = depth under review
             bf = depth set to regional depth estimate
             bx = depth reviewed but not accepted
iyr          year                                    i2
mon          month                                   i3
iday         day                                     i3
ihr          origin hour                             i4
min          origin minute                           i3
sec          origin second                           f6.2
ad           source agency                           a1
glat         geographic latitude                     f8.3
glon         geographic longitude                    f8.3
depth        depth                                   f6.1
iscdep       ISC depth                               f6.1
mb           mb magnitude                            f4.1
ms           Ms magnitude                            f4.1
mw           Mw magnitude                            f4.1
ntot         total number of observations used       i4
ntel         number of teleseismic observations      i4
             (delta > 28 deg) used                   
ndep         number of depth phase observations      i4
             (delta > 28 deg) used                   
igreg        Flinn-Engdahl region number             i4
se           standard error of observations used     f8.2
ser          standard error in position (km)         f8.2
sedep        standard error in depth (km)            f8.2
rstadel      distance to closest station             f6.1
openaz1      teleseismic azimuth gap                 f6.1
openaz2      secondary teleseismic azimuth gap       f6.1
az1          semi-axis azimuth                       f4.0
flen1        semi-axis length                        f4.1
az2          semi-axis azimuth                       f4.0
flen2        semi-axis length                        f4.1
             multiply product of len1 and len2 by pi 
             to get area of 90% confidence ellipse
avh          statistical geometric mean of axes      f5.1
iext         ISC event identifier                    i10
'''

class HDFParser:
   # the struct unpack format is as per the hdf format described above
   hdfStructFormat = '1s3s2s2s3s3s4s3s6s1s8s8s6s6s4s4s4s4s4s4s4s8s8s8s6s6s6s4s4s4s4s5s10s'

   def __init__(self, inputLine):
      self.inp = inputLine
      self._unpackEpicentre(self.inp)

   def _unpackEpicentre(self, epicentre):
      print EventParser.hdfStructFormat
      retTuple = struct.unpack(EventParser.hdfStructFormat, epicentre)
      self.ahyp = retTuple[0].strip()
      self.isol = retTuple[1].strip()
      self.iseq = retTuple[2].strip()
      self.iyr = retTuple[3].strip()
      self.mon = retTuple[4].strip()
      self.iday = retTuple[5].strip()
      self.ihr = retTuple[6].strip()
      self.min = retTuple[7].strip()
      self.sec = retTuple[8].strip()
      self.ad = retTuple[9].strip()
      self.glat = retTuple[10].strip()
      self.glon = retTuple[11].strip()
      self.depth = retTuple[12].strip()
      self.iscdep = retTuple[13].strip()
      self.mb = retTuple[14].strip()
      self.ms = retTuple[15].strip()
      self.mw = retTuple[16].strip()
      self.ntot = retTuple[17].strip()
      self.ntel = retTuple[18].strip()
      self.ndep = retTuple[19].strip()
      self.igreg = retTuple[20].strip()
      self.se = retTuple[21].strip()
      self.ser = retTuple[22].strip()
      self.sedep = retTuple[23].strip()
      self.rstadel = retTuple[24].strip()
      self.openaz1 = retTuple[25].strip()
      self.openaz2 = retTuple[26].strip()
      self.az1 = retTuple[27].strip()
      self.flen1 = retTuple[28].strip()
      self.az2 = retTuple[29].strip()
      self.flen2 = retTuple[30].strip()
      self.avh = retTuple[31].strip()
      self.iext = retTuple[32].strip()

def main():
   parser = HDFParser('ZFEQ   9  3 28  21 37 20.20   -0.770 120.269  15.0  15.0 0.0 0.0 0.0  11   0   0 265    0.00    0.00    0.00   0.5 360.0 360.0   0   0   0   0  0.0  99999999')
   print parser.__dict__

if __name__ == "__main__":
   main()
import struct
from obspy.core.event import Event, Pick, Amplitude, Magnitude, Origin

'''
Class that parses ascii data from an hdf format file
and the .out file to generate an object. The hdf format
is described below:
                            *HDF FILE FORMAT
*HDF file read statement
       read(1,100) 
     1 ahyp,isol,iseq,iyr,mon,iday,ihr,min,sec,ad,glat,glon,
     2 depth,iscdep,mb,ms,mw,ntot,ntel,ndep,igreg,se,ser,sedep,
     3 rstadel,openaz1,openaz2,az1,flen1,az2,flen2,avh,ievt
 100  format(a1,a3,a2,i2,2i3,1x,2i3,f6.2,a1,2f8.3,2f6.1,3f4.1,
     1 4i4,3f8.2,3f6.1,i4,f4.1,i4,f4.1,f5.1,i10)

Variable definitions:

ahyp         secondary teleseismic azimuth gap       a1
             blank = <= 180 deg
             Z     = >  180 deg
                    or
             primary teleseismic azimuth gap
             A     = <  180 deg
             B     = <  210 deg and  > 180 deg
             C     = <  240 deg and  > 210 deg
             D     = <  270 deg and  > 240 deg
             F     = >  270 deg
isol         solution type                           a3
             HEQ = origin time & hypocenter fixed
             DEQ = depth free (sedep <= 15 km)
             WEQ = depth fixed at waveform depth
             BEQ = depth fixed at BB depth
             FEQ = depth fixed by Engdahl
             LEQ = depth fixed by program (sedep > 15 km)
             XEQ = poor solution (ser > 35 km)
iseq         other info (b=blank)                    a2
             bx = earthquake location known to x km
             bc = only regional event data available
             Xb = explosion/cavity collapse
             Xx = explosion location known to x km
             Mb = CMT solution available
             Mx = CMT location known to x km
             Md = depth reviewed and accepted
             Mn = depth unreviewed but provisionally accepted 
                  based on => 5 depth phases and/or station(s) 
                  at distance(s) less than the focal depth
             Mr = depth under review
             Mf = depth set to regional depth estimate
             Mh = depth set to Harvard CMT depth
             Mx = depth reviewed but not accepted
             bd = depth reviewed and accepted
             bn = depth unreviewed, but provisionally accepted 
                  based on => 5 depth phases and/or station(s) 
                  at distance(s) less than the focal depth
             br = depth under review
             bf = depth set to regional depth estimate
             bx = depth reviewed but not accepted
iyr          year                                    i2
mon          month                                   i3
iday         day                                     i3
ihr          origin hour                             i4
min          origin minute                           i3
sec          origin second                           f6.2
ad           source agency                           a1
glat         geographic latitude                     f8.3
glon         geographic longitude                    f8.3
depth        depth                                   f6.1
iscdep       ISC depth                               f6.1
mb           mb magnitude                            f4.1
ms           Ms magnitude                            f4.1
mw           Mw magnitude                            f4.1
ntot         total number of observations used       i4
ntel         number of teleseismic observations      i4
             (delta > 28 deg) used                   
ndep         number of depth phase observations      i4
             (delta > 28 deg) used                   
igreg        Flinn-Engdahl region number             i4
se           standard error of observations used     f8.2
ser          standard error in position (km)         f8.2
sedep        standard error in depth (km)            f8.2
rstadel      distance to closest station             f6.1
openaz1      teleseismic azimuth gap                 f6.1
openaz2      secondary teleseismic azimuth gap       f6.1
az1          semi-axis azimuth                       f4.0
flen1        semi-axis length                        f4.1
az2          semi-axis azimuth                       f4.0
flen2        semi-axis length                        f4.1
             multiply product of len1 and len2 by pi 
             to get area of 90% confidence ellipse
avh          statistical geometric mean of axes      f5.1
iext         ISC event identifier                    i10
'''

class HDFLineParser:
   # the struct unpack format is as per the hdf format described above
   hdfStructFormat = '1s3s2s2s3s3s4s3s6s1s8s8s6s6s4s4s4s4s4s4s4s8s8s8s6s6s6s4s4s4s4s5s10s'

   def __init__(self, inputLine):
      self._unpackEpicentre(inputLine)

   def _unpackEpicentre(self, epicentre):
      retTuple = struct.unpack(HDFLineParser.hdfStructFormat, epicentre)
      self.ahyp = retTuple[0].strip()
      self.isol = retTuple[1].strip()
      self.iseq = retTuple[2].strip()
      self.iyr = retTuple[3].strip()
      self.mon = retTuple[4].strip()
      self.iday = retTuple[5].strip()
      self.ihr = retTuple[6].strip()
      self.min = retTuple[7].strip()
      self.sec = retTuple[8].strip()
      self.ad = retTuple[9].strip()
      self.glat = retTuple[10].strip()
      self.glon = retTuple[11].strip()
      self.depth = retTuple[12].strip()
      self.iscdep = retTuple[13].strip()
      self.mb = retTuple[14].strip()
      self.ms = retTuple[15].strip()
      self.mw = retTuple[16].strip()
      self.ntot = retTuple[17].strip()
      self.ntel = retTuple[18].strip()
      self.ndep = retTuple[19].strip()
      self.igreg = retTuple[20].strip()
      self.se = retTuple[21].strip()
      self.ser = retTuple[22].strip()
      self.sedep = retTuple[23].strip()
      self.rstadel = retTuple[24].strip()
      self.openaz1 = retTuple[25].strip()
      self.openaz2 = retTuple[26].strip()
      self.az1 = retTuple[27].strip()
      self.flen1 = retTuple[28].strip()
      self.az2 = retTuple[29].strip()
      self.flen2 = retTuple[30].strip()
      self.avh = retTuple[31].strip()
      self.iext = retTuple[32].strip()

