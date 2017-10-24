import struct
from obspy.core.event import Event, Pick, Amplitude, Magnitude, Origin

'''
Class that parses ascii data from one section of the
.OUT files. Each output file section looks like below:

  0110002900NORTHWEST PACIFIC EVENTS


  FEQ    iter = 0    5  9 16   h =  7 28 39.68   lat =  4.187  lon = 126.933   depth =  2.9   3.4  smdel =  0.000
  FEQ    iter = 1    5  9 16   h =  7 28 39.68   lat =  4.187  lon = 126.933   depth =  2.9   3.4  smdel =  0.000
  FEQ    iter = 2    5  9 16   h =  7 28 39.68   lat =  4.187  lon = 126.933   depth =  2.9   3.4  smdel =  0.000

  date   sta  delta  dtdd     backaz      focal angle    fm  phase               hr min sec  residual  pP   pwP  sP  PcP scor  wgt
                                           dip      az

  5 916 FITZ  22.18  10.67      3.42     33.85    183.26     P        P     P     7 33 37.06  -0.13                       0.0  0.5 G
  5 916 WR0   25.13   9.09    341.89     28.33    162.95     P        P     P     7 34  6.70   0.51                       0.0  0.6 G
  5 916 KUM   26.22  10.09     91.26     31.78    273.45     P        P     P     7 35  2.56  43.31*                      0.0  0.0 G
  5 916 MEEK  31.68   8.80     15.95     27.33    194.27     P        P     P     7 35  4.74   0.63                      -0.6  0.7 G
  5 916 FORT  34.78   8.64    358.03     26.80    178.30     P        P     P     7 35 32.05   0.65                      -0.3  0.8 G
  5 916 STKA  38.48   8.41    336.07     26.03    159.75     P        P     P     7 36  4.15   0.98                      -0.1  0.8 G
  5 916 KSM   16.82  12.63     80.43     41.23    261.25     Pn      eP    eP     7 32 39.66   3.22                       0.0  0.4
  5 916 KAKA  17.66  11.05    341.61     35.20    162.02  +  P       iP    iP     7 33  4.10  16.22*                      0.0  0.0
  5 916 GUMO  20.01   0.00    999.00    180.00     61.03              LR    LR   24  0  0.00   0.00*                      0.0  0.0
  5 916 FITZ  22.18  10.67      3.42     33.85    183.26     P       eP    eP     7 33 36.50  -0.69                       0.0  0.5
  5 916 FITZ  22.18  10.67      3.42     33.85    183.26     P        P     P     7 33 35.64  -1.55                       0.0  0.5
  5 916 WRAB  25.05   9.10    342.30     28.34    163.33     P       eP    eP     7 34  4.86  -0.60                       0.0  0.6
  5 916 WRA   25.06   9.10    342.33     28.34    163.36     P        P     P     7 34  5.09  -0.43                       0.0  0.6
                                                             S        S     S     7 38 27.07  -3.65                       0.0  0.1
  5 916 WB2   25.06   9.10    342.30     28.34    163.34     P       eP    eP     7 34  5.50  -0.05                       0.0  0.6
  5 916 JHJ   31.19   0.00    999.00    180.00     21.13              LR    LR   24  0  0.00   0.00*                      0.0  0.0
  5 916 MJAR  33.83   0.00    999.00    180.00     16.43              LR    LR   24  0  0.00   0.00*                      0.0  0.0
  5 916 FORT  34.78   8.64    358.03     26.80    178.30     P       eP    eP     7 35 31.70   0.30                      -0.3  0.8
  5 916 KLBR  36.64   8.53     15.46     26.43    193.18     P       eP    eP     7 35 46.50  -0.63                      -0.4  0.8
  5 916 MUN   37.38   8.48     17.80     26.27    195.10     P       eP    eP     7 35 53.80   0.33                      -0.4  0.8
  5 916 STKA  38.48   8.41    336.07     26.03    159.75  -  P       iP    iP     7 36  4.00   0.83                      -0.1  0.8
  5 916 STKA  38.48   8.41    336.07     26.03    159.75     P        P     P     7 36  3.10  -0.07                      -0.1  0.8
  5 916 SONM  46.93   0.00    999.00    180.00    341.12              LR    LR   24  0  0.00   0.00*                      0.0  0.0
  5 916 MKAR  57.30   7.07    123.61     21.64    324.99     P        P     P     7 38 28.53  -0.23                      -0.2  1.1
  5 916 MKAR  57.30   7.07    123.61     21.64    324.99     P        P     P     7 38 28.52  -0.24                      -0.2  1.1
  5 916 KURK  61.47   6.76    122.03     20.66    327.29     P       eP    eP     7 38 57.44   0.09                      -0.3  1.1
  5 916 BVAR  67.05   6.36    115.36     19.38    326.82     P        P     P     7 39 33.32  -0.86                      -0.1  1.2
  5 916 BRVK  67.12   6.35    115.28     19.36    326.82     P       eP    eP     7 39 34.34  -0.28                      -0.1  1.2
  5 916 VNDA  83.85   5.12    324.95     15.48    172.80     P        P     P     7 41  9.96  -1.21                       0.9  1.1
  5 916 MAW   84.25   5.08     64.35     15.38    200.26     P        P     P     7 41 13.39   0.67                       0.3  1.1
  5 916 ILAR  84.61   5.05    268.33     15.29     25.44     P        P     P     7 41 15.80   1.59                      -0.2  1.1
  5 916 KMBO  89.76   4.68     85.84     14.12    268.86     P        P     P     7 41 40.90  -1.48                       1.6  0.9

  se obs =   0.86 sec     se h =   0.45 sec     se lat =   7.02 km     se lon =  14.50 km     se depth =   0.00 km

'''

class ArrivalParser:
   # the struct unpack format of each arrival is as per the .OUT format described above
   arrivalStructFormat = '3s2s2s1x5s6s1x6s1x9s1x9s1s9s1x3s1x7s1x5s1x5s1x2s1x2s1x5s1x7s1x4s1x4s1x4s1x4s1x4s1x4s1x1s'

   def __init__(self, inputLine):
      self.inp = inputLine
      self._unpackArrival(self.inp)

   def _unpackArrival(self, arrival):
      print ArrivalParser.arrivalStructFormat
      retTuple = struct.unpack(ArrivalParser.arrivalStructFormat, arrival)
      self.year = retTuple[0].strip()
      self.month = retTuple[1].strip()
      self.day = retTuple[2].strip()
      self.station = retTuple[3].strip()
      self.delta = retTuple[4].strip()
      self.dtdd = retTuple[5].strip()
      self.backaz = retTuple[6].strip()
      self.focalDip = retTuple[7].strip()
      self.angleAzimuth = retTuple[8].strip()
      self.fm = retTuple[9].strip()
      self.phase = retTuple[10].strip()
      self.phase1 = retTuple[11].strip()
      self.phase2 = retTuple[12].strip()
      self.hour = retTuple[13].strip()
      self.minute = retTuple[14].strip()
      self.second = retTuple[15].strip()
      self.residual = retTuple[16].strip()
      self.pP = retTuple[17].strip()
      self.pwP = retTuple[18].strip()
      self.sP = retTuple[19].strip()
      self.PcP = retTuple[20].strip()
      self.scor = retTuple[21].strip()
      self.wgt = retTuple[22].strip()
      self.isGA = retTuple[23].strip()

def main():
   parser = ArrivalParser('  51025 KIPM  94.90   4.59    264.85    159.58     50.63     pP      epP   epP   19 54 31.12  -0.16* -0.2                 0.9  0.0 G')
   print parser.__dict__
   parser = ArrivalParser('                                                             PcP      PcP   PcP  11 28 54.70   1.86*           -9.8  1.9  0.6  0.0  ', )
   print parser.__dict__

if __name__ == "__main__":
   main()
