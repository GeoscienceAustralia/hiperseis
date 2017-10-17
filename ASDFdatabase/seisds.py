import numpy as np
import time
import json
import pyasdf
from obspy import UTCDateTime, Stream, Trace
import random


class SeisDB(object):
    def __init__(self, json_file=False, generate_numpy_index=True):
        self._json_loaded = False
        self._valid_index = False
        self._use_numpy_index = False
        self._generate_numpy_index = generate_numpy_index
        if json_file:
            try:
                f = open(json_file, 'r')
                self._json_dict = json.load(f)
                self._json_file = json_file
                self._json_loaded = True
                self.generateIndex()
            except IOError as e:
                print "I/O error({0}): {1}".format(e.errorno, e.strerror)
            except ValueError as e:
                print "JSON Decoding has failed with a value error({0}): {1}".format(e.errorno, e.strerror)

    def generateIndex(self):
        assert self._json_loaded, "Invalid SeisDB object. Try loading a valid JSON file first."

        if self._generate_numpy_index:

            try:
                # dictionary with keys as integer index (corresponding to numpy array indexes) and value which is ASDF tag
                self._indexed_dict = {}
                # List for numpy array building
                self._index_dict_list = []
                # container for dtypes
                type_list = []
                # check if the dtype has been populated
                dtype_pop = False
                for _i, (key, value) in enumerate(self._json_dict.iteritems()):
                    self._indexed_dict[_i] = key
                    temp_list = []
                    for _j, (sub_key, sub_value) in enumerate(value.iteritems()):
                        # only add some of the attributes to the numpy array to speed up lookup

                        if sub_key == "tr_starttime":
                            temp_list.append(float(sub_value))
                            if not dtype_pop:
                                type_list.append(('st', float))
                        elif sub_key == "tr_endtime":
                            temp_list.append(float(sub_value))
                            if not dtype_pop:
                                type_list.append(('et', float))
                        elif sub_key == "new_network":
                            temp_list.append(str(sub_value))
                            if not dtype_pop:
                                type_list.append(('net', 'S2'))
                        elif sub_key == "new_station":
                            temp_list.append(str(sub_value))
                            if not dtype_pop:
                                type_list.append(('sta', 'S5'))
                        elif sub_key == "new_channel":
                            temp_list.append(str(sub_value))
                            if not dtype_pop:
                                type_list.append(('cha', 'S3'))
                        elif sub_key == "new_location":
                            temp_list.append(str(sub_value))
                            if not dtype_pop:
                                type_list.append(('loc', 'S2'))

                    dtype_pop = True

                    self._index_dict_list.append(tuple(temp_list))
                dt = np.dtype(type_list)
                self._indexed_np_array = np.array(self._index_dict_list, dtype=dt)
                self._use_numpy_index = True
                self._valid_index = True

            except KeyError as e:
                print "Indexing JSON dictionary has failed with a key error({0}): {1}".format(e.errorno, e.strerror)

        else:
            try:
                # dictionary with keys as integer index (corresponding to numpy array indexes) and value which is ASDF tag
                self._indexed_dict = {}
                # new dictionary to be sure that starttime and endtime fields are float
                self._formatted_dict = {}
                for _i, (key, value) in enumerate(self._json_dict.iteritems()):
                    self._indexed_dict[_i] = key
                    temp_dict = {}
                    for _j, (sub_key, sub_value) in enumerate(value.iteritems()):
                        if sub_key == "tr_starttime":
                            temp_dict[sub_key] = float(sub_value)
                        elif sub_key == "tr_endtime":
                            temp_dict[sub_key] = float(sub_value)
                        else:
                            temp_dict[sub_key] = sub_value

                    self._formatted_dict[_i] = temp_dict
                self._valid_index = True
            except KeyError as e:
                print "Indexing JSON dictionary has failed with a key error({0}): {1}".format(e.errorno, e.strerror)

    def queryByTime(self, sta, chan, query_starttime, query_endtime):
        qs = query_starttime
        qe = query_endtime
        assert self._json_loaded, "Invalid SeisDB object. Try loading a valid JSON file first."
        assert self._valid_index, "Invalid SeisDB object. Index has not been generated."
        if not self._use_numpy_index:
            indices = []
            for _i, key in enumerate(self._formatted_dict.keys()):
                matched_entry = self._formatted_dict[key]
                if ((matched_entry['tr_starttime'] <= qs < matched_entry['tr_endtime'])
                    or (qs <= matched_entry['tr_starttime'] and matched_entry['tr_starttime'] < qe)) \
                        and ((matched_entry['new_station'] in sta) and (matched_entry['new_channel'] in chan)):
                    indices.append(_i)
            # indices_array = np.array(indices)
            # Print output
            # print(indices_array)
            # for index in indices_array:
            #    print(self._indexed_dict[index]['ASDF_tag'])
            # print(len(indices_array))
            # return {k:self._indexed_dict[self._index_dict_list[k]] for k in indices if k in self._index_dict_list}
            return {k: {"ASDF_tag": self._indexed_dict[k],
                        "new_station": self._formatted_dict[k]["new_station"],
                        "new_network": self._formatted_dict[k]["new_network"]}
                    for k in indices}
        else:
            _indexed_np_array_masked = np.where((np.in1d(self._indexed_np_array['sta'], sta))
                           & (np.in1d(self._indexed_np_array['cha'], chan))
                           & np.logical_or(np.logical_and(self._indexed_np_array['st'] <= qs,  qs < self._indexed_np_array['et']),
                                           (np.logical_and(qs <= self._indexed_np_array['st'],
                                                           self._indexed_np_array['st'] < qe))))
            # print(_indexed_np_array_masked[0])
            # for index in _indexed_np_array_masked[0]:
            #    print(self._indexed_np_array[index, 6])
            # print(len(_indexed_np_array_masked[0]))
            # print(self._index_dict_list[0])
            # return {k:self._indexed_dict[self._indexed_np_array[k]] for k in _indexed_np_array_masked[0] if k in self._indexed_np_array}
            return {k: {"ASDF_tag": self._indexed_dict[k],
                        "new_station": self._indexed_np_array['sta'][k],
                        "new_network": self._indexed_np_array['net'][k]}
                    for k in _indexed_np_array_masked[0]}

    def queryToStream(self, ds, query, query_starttime, query_endtime,
                      decimation_factor=None):
        """
        Method to use output from seisds query to return waveform streams

        :type ds: ASDFDataSet
        :param ds: ASDFDataSet to fetch waveforms from
        :type query: dict
        :param query: Dictionary returned by function 'queryByTime'
        :type query_starttime: UTCDateTime or str
        :param query_starttime: Start time of waveform
        :type query_endtime: UTCDateTime or str
        :param query_endtime: End time of waveform
        :type decimation_factor: int
        :param decimation_factor: Decimation factor applied to returned waveforms

        :return : Stream object with merged waveform data, in which gaps are masked
        """

        # Open a new st object
        st = Stream()

        for matched_info in query.values():

            # read the data from the ASDF into stream
            temp_tr = ds.waveforms[matched_info["new_network"] + '.' + matched_info["new_station"]][
                matched_info["ASDF_tag"]][0]

            # trim trace to start and endtime
            temp_tr.trim(starttime=UTCDateTime(query_starttime),
                         endtime=UTCDateTime(query_endtime))

            # append the asdf id tag into the trace stats so that the original data is accesbale
            temp_tr.stats.asdf.orig_id = matched_info["ASDF_tag"]

            # Decimate trace
            if(decimation_factor is not None): temp_tr.decimate(decimation_factor)

            # append trace to stream
            st += temp_tr

            # free memory
            temp_tr = None
        #end for

        if st.__nonzero__():
            # Attempt to merge all traces with matching ID'S in place
            #print('')
            #print('Merging %d Traces ....' % len(st))
            # filling no data with 0
            #st.print_gaps()
            st.merge()
            return st
        else:
            return None
    #end func

    def fetchDataByTime(self, ds, sta, chan, query_starttime, query_endtime,
                        decimation_factor=None):
        """
        Wrapper to use output from seisds query to return waveform streams

        :type ds: ASDFDataSet
        :param ds: ASDFDataSet to fetch waveforms from
        :type sta: list
        :param sta: List of station names
        :type chan: str
        :param chan: Wildcard to be used to filter channels
        :type query_starttime: UTCDateTime or str
        :param query_starttime: Start time of waveform
        :type query_endtime: UTCDateTime or str
        :param query_endtime: End time of waveform
        :type decimation_factor: int
        :param decimation_factor: Decimation factor applied to returned waveforms

        :return : Stream object with merged waveform data, in which gaps are masked
        """

        qr = self.queryByTime(sta, chan, query_starttime, query_endtime)

        st = self.queryToStream(ds, qr, query_starttime, query_endtime,
                                decimation_factor=decimation_factor)

        return st
    #end func
#end class

if __name__ == "__main__":
    print "Testing db access"
    #sta = ["SQ2A1", "SQ2F6"]
    sta = ["CP43", "CQ43"]
    chan = ["BHZ"]
    query_stime = UTCDateTime(2015, 1, 1, 00).timestamp
    query_etime = UTCDateTime(2015, 1, 2, 00).timestamp

    db = '/g/data1/ha3/Passive/_ANU/7F(2013-2014)/raw_DATA/7F(2013-2014)_raw_dataDB.json'

    #db = '/g/data/ha3/rakib/_ANU/7G(2013-2015)/ASDF/7G(2013-2015)_raw_dataDB.json'

    db = '/g/data1/ha3/Passive/_ANU/7G(2013-2015)/ASDF/7G(2013-2015)_raw_dataDB.json'
    ds = '/g/data1/ha3/Passive/_ANU/7G(2013-2015)/ASDF/7G(2013-2015).h5'

    print "Sta = " + str(sta)
    print "Chan = " + str(chan)
    print "Query start time = " + str(query_stime)
    print "Query end time = " + str(query_etime)
    print "JSON file = " + str(db)
    sdb = SeisDB(db)

    '''
    res = sdb.queryByTime(sta, chan, query_stime, query_etime)
    # print res
    for index, element in res.iteritems():
        print "Index = " + str(index)
        print "Element = " + str(element)
    '''

    h5ds = pyasdf.ASDFDataSet(ds)
    st = sdb.fetchDataByTime(h5ds, sta, chan, query_stime, query_etime)
    del h5ds