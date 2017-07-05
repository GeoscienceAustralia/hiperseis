import json
from obspy import read_inventory, read_events, UTCDateTime, Stream, read
from obspy.clients.fdsn.client import Client
from obspy.clients.fdsn.header import FDSNException
import pyasdf
import itertools
import sys
import os

class Client2ASDF(object):
	def __init__(self,refdir,filename,client="IRIS",network="AU"):
		"""Object to query IRIS or other client with a bounding box and time interval.
		   Returns ASDF containing station information and data for interval."""
		self._client = client
		self._network = network
		#self._query_ds = None
		self._refdir = refdir
		self._asdffilename = filename
		self.ref_stations = []
	def queryByBBoxInterval(self,bbox,timeinterval,bbpadding=2,verbose=False,event_id=None):
		""" Time interval is a tuple (starttime,endtime)
		"""
		assert len(timeinterval)==2,"timeinterval must be a tuple of ascending timestamps. len="+str(len(timeinterval))+" "+str(timeinterval)

		query_ds = pyasdf.ASDFDataSet(self._asdffilename)

		client = Client(self._client)
		ref_inv = client.get_stations(network=self._network,
						   starttime=UTCDateTime(timeinterval[0]),
						   endtime=UTCDateTime(timeinterval[1]),
						   minlongitude=bbox[0]-bbpadding,
						   maxlongitude=bbox[1]+bbpadding,
						   minlatitude=bbox[2]-bbpadding,
						   maxlatitude=bbox[3]+bbpadding,
						   level='channel')

		if verbose:
			print(ref_inv)

		ref_st = Stream()
		inv = []

		# go through inventory and request timeseries data
		for net in ref_inv:
			for stn in net:
				try:
					ref_st += client.get_waveforms(network=net.code, station=stn.code, channel='*', location='*',
							   starttime=UTCDateTime(timeinterval[0]),
							   endtime=UTCDateTime(timeinterval[1]))
				except FDSNException:
					print('No Data for Earthquake from Reference Station: ' + stn.code)

				else:
					self.ref_stations.append(net.code + '.' + stn.code)
					# append the inv
					inv.append(ref_inv.select(station=stn.code))

		# write ref traces into ASDF
		for st_inv in inv:
			query_ds.add_stationxml(st_inv)
		for tr in ref_st:
			query_ds.add_waveforms(tr, "reference_station", event_id=event_id)
			#tr.write(os.path.join(self._refdir, tr.id + ".MSEED"), format="MSEED") # Don't write miniseed
		if verbose:
			print("Wrote Reference Waveforms to ASDF file: " + temp_ASDF_out)
			print('\nEarthquake data query completed.')

		ref_inv.write(os.path.join(self._refdir, "ref_metadata.xml"), format="STATIONXML")

if __name__ == "__main__":
	dr = './'
	fn = 'tempasdf'
	pth = os.path.join(dr,fn)
	c = Client2ASDF(dr,fn)
	inka = (-27.741 , 140.746)
	p = 0.0001
	c.queryByBBoxInterval([inka[1]-p,inka[1]+p,inka[0]-p,inka[0]+p],("2014-06-01T00:00:00","2014-06-02T00:00:00"),1)
	# test that the file exists and was written
	assert os.path.isfile(pth), "ASDF file not written"
	os.remove(pth)
