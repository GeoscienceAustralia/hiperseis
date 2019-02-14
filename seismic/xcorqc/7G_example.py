from obspy import UTCDateTime
import pyasdf

from seismic.ASDFdatabase.seisdb import SeisDB
from xcorqc import IntervalStackXCorr

# =========================== User Input Required =========================== #

# == Path to Temporary Station Data ==
tempDataFile        = '/g/data/ha3/Passive/_ANU/7G(2013-2015)/ASDF/7G(2013-2015).h5'
# Associated json database file
tempDataDBFile      = '/g/data/ha3/Passive/_ANU/7G(2013-2015)/ASDF/7G(2013-2015)_raw_dataDB.json'
# Station codes to analyse (set to ['*']) for all stations
tempStationCodes    =  ['CP43', 'CQ43']
# Decimation factor applied to traces
tempDecFactor       = 5

# == Time-Range, Windowing and Decimation ==
startTime           = "2015-01-01T00:00:00"
endTime             = "2015-06-01T00:00:00"
# Size of data-buffer (we don't want to fetch all data at once)
bufferSeconds       = 3600*24*10 # 10 days (should be a multiple of interval Seconds)
# Interval over which windows are stacked
intervalSeconds     = 3600*24*10 # 10 day
# Size of data-window to be stacked
windowSeconds       = 3600

# == Path to Reference Station Data ==
refDataFile         = '/g/data/ha3/rakib/_ANU/7G(2013-2015)/refData/stka.6m.h5'
# Station codes to analyse (set to ['*']) for all stations
refStationCodes     = ['*'] #['INKA']
# Decimation factor applied to traces
refDecFactor        = 5

# == Component to analyse ==
comp                = 'BHZ'

# == Output Path ==
outputDir           = '/g/data/ha3/Passive/newXcorr/'

# =========================== End of User Input ============================= #





# =========================== Begin Processing ============================= #

# Open data sets
refDs       = None
tempDs      = None
tempDsDb    = None

try:
    refDs       = pyasdf.ASDFDataSet(refDataFile, mode='r')
    tempDs      = pyasdf.ASDFDataSet(tempDataFile, mode='r')
    tempDsDb    = SeisDB(tempDataDBFile)
except:
    msg = 'Failed to open data sets. Check input file names..'
    raise Exception(msg)
#end try


startTime = UTCDateTime(startTime)
endTime   = UTCDateTime(endTime)

x, xCorrResDict, wcResDict = IntervalStackXCorr(refDs, tempDs, tempDsDb, startTime, endTime,
                                                refStationCodes, tempStationCodes, comp, refDecFactor,
                                                tempDecFactor, bufferSeconds, intervalSeconds,
                                                windowSeconds, flo=0.5, fhi=3, outputPath=outputDir,
                                                verbose=2)
