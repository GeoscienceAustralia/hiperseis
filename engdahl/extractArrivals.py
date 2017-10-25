from ArrivalParser import ArrivalParser
from HDFLineParser import HDFLineParser
from itertools import groupby
import linecache

def getArrivalSections(hdfFile, outFile):
   with open(outFile, 'r') as of:
      with open(hdfFile, 'r') as hf:
         count = 0
         for keyvaluetuple in groupby(of, lambda x: "PACIFIC EVENTS" in x.strip()):
            key, value = keyvaluetuple
            if not key:
               count = count + 1
               eventParser = HDFLineParser(linecache.getline(hdfFile, count).rstrip('\n'))
               print eventParser.__dict__
               for groupline in value:
                  if len(groupline) == 133:
                     arrivalParser = ArrivalParser(groupline.rstrip('\n'))
                     print arrivalParser.__dict__


def main():
   getArrivalSections('EHB.GA.HDF', 'EHB.GA.OUT')

if __name__ == "__main__":
   main()
