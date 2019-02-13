import numpy as np
import time
import json
from obspy import UTCDateTime
import random
import logging

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
                logging.error("I/O error({0}): {1}".format(e.errorno, e.strerror))
            except ValueError as e:
                logging.error("JSON Decoding has failed with a value error({0}): {1}".format(e.errorno, e.strerror))

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
                for _i, (key, value) in enumerate(self._json_dict.items()):
                    self._indexed_dict[_i] = key
                    temp_list = []

                    # get ASDF tag
                    tag = key.split('__')[3]

                    temp_list.append(tag)
                    if not dtype_pop:
                        type_list.append(('tag', 'S100'))

                    logging.info(value)
                    logging.info(type(value))

                    for sub_key, sub_value in value.items():
                        logging.info("{}, {}".format(sub_key, sub_value))

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

                self.dt = np.dtype(type_list)
                self._indexed_np_array = np.array(self._index_dict_list, dtype=self.dt)
                logging.info(self._indexed_np_array)
                self._use_numpy_index = True
                self._valid_index = True

            except KeyError as e:
                logging.error("Indexing JSON dictionary has failed with a key error({0}): {1}".format(e.errorno, e.strerror))

        else:
            try:
                # dictionary with keys as integer index (corresponding to numpy array indexes) and value which is ASDF tag
                self._indexed_dict = {}
                # new dictionary to be sure that starttime and endtime fields are float
                self._formatted_dict = {}
                for _i, (key, value) in enumerate(self._json_dict.items()):
                    self._indexed_dict[_i] = key
                    temp_dict = {}
                    for _j, (sub_key, sub_value) in enumerate(value.items()):
                        if sub_key == "tr_starttime":
                            temp_dict[sub_key] = float(sub_value)
                        elif sub_key == "tr_endtime":
                            temp_dict[sub_key] = float(sub_value)
                        else:
                            temp_dict[sub_key] = sub_value

                    self._formatted_dict[_i] = temp_dict
                self._valid_index = True
            except KeyError as e:
                logging.info("Indexing JSON dictionary has failed with a key error({0}): {1}".format(e.errorno, e.strerror))

    def queryByTime(self, net, sta, chan, tags, query_starttime, query_endtime):
        qs = query_starttime
        qe = query_endtime
        assert self._json_loaded, "Invalid SeisDB object. Try loading a valid JSON file first."
        assert self._valid_index, "Invalid SeisDB object. Index has not been generated."
        if not self._use_numpy_index:
            indices = []
            for _i, key in enumerate(self._formatted_dict.keys()):
                # check if tag matches
                # get ASDF tag
                tag = key.split('__')[3]
                if not tag in tags:
                    continue
                matched_entry = self._formatted_dict[key]
                if ((matched_entry['tr_starttime'] <= qs < matched_entry['tr_endtime'])
                    or (qs <= matched_entry['tr_starttime'] and matched_entry['tr_starttime'] < qe)) \
                        and ((matched_entry['new_network'] in net) and (matched_entry['new_station'] in sta) and (
                                    matched_entry['new_channel'] in chan)):
                    indices.append(_i)
            # indices_array = np.array(indices)
            # Print output
            # logging.debug(indices_array)
            # for index in indices_array:
            #    logging.debug(self._indexed_dict[index]['ASDF_tag'])
            # logging.debug(len(indices_array))
            # return {k:self._indexed_dict[self._index_dict_list[k]] for k in indices if k in self._index_dict_list}
            return {k: {"ASDF_tag": self._indexed_dict[k],
                        "new_station": self._formatted_dict[k]["new_station"],
                        "new_network": self._formatted_dict[k]["new_network"]}
                    for k in indices}
        else:
            _indexed_np_array_masked = np.where(
                (np.in1d(self._indexed_np_array['net'], net)) & (np.in1d(self._indexed_np_array['sta'], sta))
                & (np.in1d(self._indexed_np_array['cha'], chan)) & (np.in1d(self._indexed_np_array['tag'], tags))
                & np.logical_or(np.logical_and(self._indexed_np_array['st'] <= qs, qs < self._indexed_np_array['et']),
                                (np.logical_and(qs <= self._indexed_np_array['st'],
                                                self._indexed_np_array['st'] < qe))))
            # logging.debug(_indexed_np_array_masked[0])
            # for index in _indexed_np_array_masked[0]:
            #    logging.debug(self._indexed_np_array[index, 6])
            # logging.debug(len(_indexed_np_array_masked[0]))
            # logging.debug(self._index_dict_list[0])
            # return {k:self._indexed_dict[self._indexed_np_array[k]] for k in _indexed_np_array_masked[0] if k in self._indexed_np_array}
            return {k: {"ASDF_tag": self._indexed_dict[k],
                        "new_station": self._indexed_np_array['sta'][k],
                        "new_network": self._indexed_np_array['net'][k]}
                    for k in _indexed_np_array_masked[0]}

    def get_recording_intervals(self, net, sta, chan, tags):
        # logging.debug(self._indexed_np_array)
        logging.info(net, sta, chan, tags)
        _indexed_np_array_masked = np.where(
            (np.in1d(self._indexed_np_array['net'], net)) & (np.in1d(self._indexed_np_array['sta'], sta)) & (
                np.in1d(self._indexed_np_array['cha'], chan)) & (np.in1d(self._indexed_np_array['tag'], tags)))

        _masked_np_array = self._indexed_np_array[_indexed_np_array_masked]

        _st_array = _masked_np_array["st"]
        _et_array = _masked_np_array["et"]

        _arg_sorted_st_array = np.argsort(_st_array)

        _sorted_st_array = _st_array[_arg_sorted_st_array]
        _sorted_et_array = _et_array[_arg_sorted_st_array]

        logging.info(_sorted_st_array)
        logging.info(_sorted_et_array)

        if _sorted_st_array.shape[0] == 1:
            #i.e. there are no gaps
            logging.info(_sorted_st_array[0])
            return (np.array([[_sorted_st_array[0]], [_sorted_et_array[0]]]))


        # offset the starttime array so that we start from the second recorded waveform
        _offset_st_array = _sorted_st_array[1:]

        logging.info(_offset_st_array)

        _diff_array = _offset_st_array - _sorted_et_array[:-1]

        logging.info(_diff_array)

        # now get indexes when gap (diff positive) or overlap (diff negative)
        # remember to add 1 to index because of offset
        _sorted_gaps_index = np.where(_diff_array > 1)

        logging.info(_sorted_gaps_index)

        # if it is empty then the gaps are smaller than 1 second ignore them
        if _sorted_gaps_index[0].shape[0]==0:
            return (np.array([[_sorted_st_array[0]], [_sorted_et_array[-1]]]))



        # recording intervals:
        # rec_start_list_after gaps = list(_sorted_st_array[_sorted_gaps_index])
        rec_end_list_after_gaps = list(_sorted_et_array[_sorted_gaps_index])

        logging.info(rec_end_list_after_gaps)


        rec_start_list =[]
        rec_end_list =[]

        int_id = 0
        for _i in range(len(rec_end_list_after_gaps)):
            rec_start_list.append(_sorted_st_array[int_id])
            rec_end_list.append(rec_end_list_after_gaps[_i])
            int_id = _sorted_gaps_index[0][_i] + 1
            if _i == len(rec_end_list_after_gaps) - 1:
                # last interval
                # logging.debug(_sorted_st_array[int_id], _sorted_et_array[-1])
                rec_start_list.append(_sorted_st_array[int_id])
                rec_end_list.append(_sorted_et_array[-1])

        return (np.array([rec_start_list, rec_end_list]))

    def get_data_percentage(self, net, sta, chan, tags):

        # call intervals calc
        intervals_array = self.get_recording_intervals(net, sta, chan, tags)

        gaps_array = self.get_gaps_intervals(net, sta, chan, tags)


        logging.info("Getting data percentage")

        logging.info(str([intervals_array[1][x] - intervals_array[0][x] for x in range(len(intervals_array[0]))]))

        total_rec_seconds = sum([intervals_array[1][x] - intervals_array[0][x] for x in range(len(intervals_array[0]))])
        total_gaps_seconds = sum([gaps_array[1][x] - gaps_array[0][x] for x in range(len(gaps_array[0]))])

        logging.info("{}, {}".format(total_rec_seconds, total_gaps_seconds))
        logging.info("Percentage of Gaps: {}".format(100*(total_gaps_seconds/total_rec_seconds)))
        logging.info("Percentage of Recording: {}".format(100 * (1-(total_gaps_seconds / total_rec_seconds))))

    def get_gaps_intervals(self, net, sta, chan, tags):

        # logging.debug(self._indexed_np_array)
        _indexed_np_array_masked = np.where(
            (np.in1d(self._indexed_np_array['net'], net)) & (np.in1d(self._indexed_np_array['sta'], sta)) & (
            np.in1d(self._indexed_np_array['cha'], chan)) & (np.in1d(self._indexed_np_array['tag'], tags)))
        # logging.debug(_indexed_np_array_masked)

        _masked_np_array = self._indexed_np_array[_indexed_np_array_masked]
        # logging.debug(_masked_np_array)

        # logging.debug(str(np.sort(_masked_np_array, order='st')))

        _st_array = _masked_np_array["st"]
        _et_array = _masked_np_array["et"]

        # logging.debug(_st_array)

        _arg_sorted_st_array = np.argsort(_st_array)

        _sorted_st_array = _st_array[_arg_sorted_st_array]
        _sorted_et_array = _et_array[_arg_sorted_st_array]

        # logging.debug(_arg_sorted_st_array)
        # logging.debug(_sorted_st_array)
        # logging.debug(_sorted_et_array)

        # offset the starttime array so that we start from the second recorded waveform
        _offset_st_array = _sorted_st_array[1:]

        # logging.debug(_offset_st_array)

        _diff_array = _offset_st_array - _sorted_et_array[:-1]

        # logging.debug(_diff_array)

        # now get indexes when gap (diff positive) or overlap (diff negative)
        # remember to add 1 to index because of offset
        _sorted_gaps_index = np.where(_diff_array > 1)
        # _sorted_ovlps_index = np.where(_diff_array < 1)

        # logging.debug(_sorted_gaps_index)

        # now get the gap_starttimes and gap endtimes
        gaps_start_list = list(_sorted_et_array[_sorted_gaps_index])
        gaps_end_list = list(_offset_st_array[_sorted_gaps_index])

        return (np.array([gaps_start_list, gaps_end_list]))

    def get_unique_information(self):
        """
        Method to retreive the unique channels and tags within an ASDF file from the JSON Database
        :return: (unique channels, unique tags)
        """
        assert self._json_loaded, "Invalid SeisDB object. Try loading a valid JSON file first."
        assert self._valid_index, "Invalid SeisDB object. Index has not been generated."
        if not self._use_numpy_index:
            logging.error("Must be using Numpy index")

        else:
            return (np.unique(self._indexed_np_array['cha']), np.unique(self._indexed_np_array['tag']))

    def retrieve_full_db_entry(self, json_key):
        """
        Method to take an output from queryByTime and get the full information from the original JSON db
        :return: full DB
        """

        # logging.debug(self._json_dict[json_key])
        return (self._json_dict[json_key])

    def is_chan_related(self, chan, net, sta, loc):
        """
        Method to test if a channel code is related to a given net/stn/tag/loc
        :param chan: Channel Code, String
        :param net: Network Code, String
        :param sta: Station Code, String
        :param sta: Location Code, String
        :return: True/False, Bool
        """
        _indexed_np_array_masked = np.where(
            (np.in1d(self._indexed_np_array['net'], net)) & (np.in1d(self._indexed_np_array['sta'], sta))
            & (np.in1d(self._indexed_np_array['cha'], chan)) & (np.in1d(self._indexed_np_array['loc'], loc)))

        _masked_np_array = self._indexed_np_array[_indexed_np_array_masked]

        # returns true if there are traces in ASDF with chanel relating, or false otherwise
        return(bool(_masked_np_array.size))


if __name__ == "__main__":
    print("Testing db access")
    sta = ["SQ2A1", "SQ2F6"]
    chan = ["BHZ"]
    query_time = UTCDateTime(2014, 3, 17, 00).timestamp
    db = '/g/data1/ha3/Passive/_ANU/7F(2013-2014)/raw_DATA/7F(2013-2014)_raw_dataDB.json'
    print("Sta = " + str(sta))
    print("Chan = " + str(chan))
    print("Query time = " + str(query_time))
    print("JSON file = " + str(db))
    sdb = SeisDB(db)
    res = sdb.queryByTime(sta, chan, query_time)
    for index, element in res.items():
        logging.info("Index = " + str(index))
        logging.info("Element = " + str(element))

