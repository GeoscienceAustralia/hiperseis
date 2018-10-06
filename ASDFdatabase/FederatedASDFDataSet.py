"""
Description:
    Class for providing fast access to data contained within a set of ASDF files

References:

CreationDate:   03/09/18
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     03/19/18   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

from mpi4py import MPI
import os
import glob
import atexit
import logging
import pickle
import numpy as np

from obspy.core import Stream, UTCDateTime
from obspy import read, Trace
import pyasdf
import ujson as json
from collections import defaultdict
from rtree import index

def setup_logger(name, log_file, level=logging.INFO):
    """
    Function to setup a logger; adapted from stackoverflow
    """
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    handler = logging.FileHandler(log_file, mode='w')
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)
    return logger
# end func

def tree():
    def the_tree():
        return defaultdict(the_tree)
    # end func
    return the_tree()
# end func

def split_list(lst, npartitions):
    result = []
    for i in np.arange(npartitions):
        result.append([])
    # end for
    count = 0
    for iproc in np.arange(npartitions):
        for i in np.arange(np.divide(len(lst), npartitions)):
            result[iproc].append(lst[count])
            count += 1
    # end for
    for iproc in np.arange(np.mod(len(lst), npartitions)):
        result[iproc].append(lst[count])
        count += 1
    # end for

    return result
# end func

class FederatedASDFDataSet():
    @staticmethod
    def hasOverlap(stime1, etime1, stime2, etime2):
        result = 0

        #print stime1, etime1, stime2, etime2

        if (etime1 is None): etime1 = UTCDateTime.now().timestamp
        if (etime2 is None): etime2 = UTCDateTime.now().timestamp

        if (stime2 >= stime1 and stime2 <= etime1):
            result = 1
        elif (etime2 >= stime1 and etime2 <= etime1):
            result = 1
        elif (stime2 <= stime1 and etime2 >= etime1):
            result = 1

        return result

    # end func

    def __init__(self, asdf_source, use_json_db=False, logger=None):
        """

        :param asdf_source: path to a text file containing three space-delimited columns:
               1. path to an ASDF file
               2. Start-time (in UTCDateTime format) of data in the ASDF file (used for speeding up data access)
               3. End-time (in UTCDateTime format) of data in the ASDF file (used for speeding up data access)
               Entries can be commented out with '#'
        """
        self.comm = MPI.COMM_WORLD
        self.nproc = self.comm.Get_size()
        self.rank = self.comm.Get_rank()

        self.use_json_db = use_json_db
        self.logger = logger
        self.asdf_source = None
        self.asdf_file_names = []
        self.json_db_file_names = []
        self.tree_list = []
        self.asdf_tags_list = []
        self.asdf_station_coordinates = []
        self.net_sta_start_time = tree()
        self.net_sta_end_time = tree()
        self.waveform_start_time_list = []
        self.waveform_end_time_list = []

        if(type(asdf_source)==str):
            self.asdf_source = asdf_source

            fileContents = filter(len, open(self.asdf_source).read().splitlines())

            # collate filename and time bounds
            for i in range(len(fileContents)):
                if(fileContents[i][0]=='#'): continue # filter commented lines

                fn = fileContents[i].strip(' \t\n\r\n')
                self.asdf_file_names.append(fn)

                # look for json db
                files = glob.glob(os.path.dirname(fn) + '/*.json')
                found_json_db = False
                for f in files:
                    if (os.path.splitext(os.path.basename(fn))[0] in f):
                        found_json_db = True
                        self.json_db_file_names.append(f)
                        break
                    # end if
                # end for
                if(not found_json_db): self.json_db_file_names.append(None)
            # end for
        else:
            raise NameError('Invalid value for asdf_source..')
        # end if

        self.asdf_datasets = []
        for ifn, fn in enumerate(self.asdf_file_names):
            if self.logger:self.logger.info('Opening ASDF file %s..'%(fn))

            if(os.path.exists(fn)):
                ds = pyasdf.ASDFDataSet(fn, mode='r')
                self.asdf_datasets.append(ds)
                self.asdf_station_coordinates.append(ds.get_all_coordinates())
            else:
                raise NameError('File not found: %s..'%fn)
            # end if

            self.asdf_tags_list.append([])
            self.tree_list.append(None)
            self.waveform_start_time_list.append(None)
            self.waveform_end_time_list.append(None)
        # end func

        # Create indices
        if(self.use_json_db): self.create_index()
        else: self.create_sparse_index()

        atexit.register(self.cleanup) # needed for closing asdf files at exit
    # end func

    def extract_time_ranges(self):
        print 'Extracting time ranges on rank %d'%(self.rank)
        for it, val in enumerate(self.tree_list):
            if(val):
                min = []
                max = []
                for (nk, nv) in val.iteritems():
                    for (sk, sv) in nv.iteritems():

                        if (type(self.net_sta_start_time[nk][sk]) == defaultdict):
                            self.net_sta_start_time[nk][sk] = 1e32
                        # end if
                        if (type(self.net_sta_end_time[nk][sk]) == defaultdict):
                            self.net_sta_end_time[nk][sk] = -1e32
                        # end if

                        for (lk, lv) in sv.iteritems():
                            for (ck, cv) in lv.iteritems():
                                if (type(cv) == index.Index):
                                    if (len(cv.leaves()[0][1])==0): continue

                                    bounds = cv.bounds

                                    if(self.net_sta_start_time[nk][sk] > bounds[0]):
                                        self.net_sta_start_time[nk][sk] = bounds[0]
                                    if(self.net_sta_end_time[nk][sk] < bounds[2]):
                                        self.net_sta_end_time[nk][sk] = bounds[2]

                                    min.append(bounds[0])
                                    max.append(bounds[2])
                                # end if
                            # end for
                        # end for
                    # end for
                # end for
                self.waveform_start_time_list[it] = np.min(np.array(min))
                self.waveform_end_time_list[it] = np.max(np.array(max))
            # end if
        # end for
    # end func

    def create_index(self):
        for ifn, fn in enumerate(self.json_db_file_names):
            if(fn):
                print 'Creating index for %s..'%(fn)
                d = json.load(open(fn, 'r'))

                t = tree()
                fieldMap = None
                indices = []
                for ki, (k, v) in enumerate(d.iteritems()):

                    if (fieldMap == None):
                        fieldMap = defaultdict()
                        for i, sk in enumerate(v.keys()):
                            fieldMap[sk] = i
                        # end for
                        indices = [fieldMap['new_network'],
                                   fieldMap['new_station'],
                                   fieldMap['new_channel'],
                                   fieldMap['new_location'],
                                   fieldMap['tr_starttime'],
                                   fieldMap['tr_endtime']]
                    # end if

                    values = v.values()
                    network = values[indices[0]]
                    station = values[indices[1]]
                    channel = values[indices[2]]
                    location = values[indices[3]]
                    tr_st = values[indices[4]]
                    tr_et = values[indices[5]]

                    if (type(t[network][station][location][channel]) == defaultdict):
                        t[network][station][location][channel] = index.Index()
                    # end if

                    t[network][station][location][channel].insert(ki, (tr_st, 1, tr_et, 1))
                # end for

                self.tree_list[ifn] = t
                self.asdf_tags_list[ifn] = d.keys()

                del d
                print 'Done creating index'
            # end if
        # end for
        self.extract_time_ranges()
    # end func

    def create_sparse_index(self, balance_by_tags=True):
        def decode_tag(tag, type='raw_recording'):
            if (type not in tag): return None
            try:
                tokens = tag.split('.')
                nc, sc, lc = tokens[0], tokens[1], tokens[2]

                tokens = tokens[3].split('__')
                cc = tokens[0]
                starttime = UTCDateTime(tokens[1]).timestamp
                endttime  = UTCDateTime(tokens[2]).timestamp

                return nc, sc, lc, cc, starttime, endttime
            except:
                return None
            # end func
        # end func

        tagsCount = 0
        for ids, ds in enumerate(self.asdf_datasets):
            print 'Creating index for %s..' % (self.asdf_file_names[ids])

            keys = ds.get_all_coordinates().keys()
            keys = split_list(keys, self.nproc)

            data = []
            #print 'Found %d keys'%(len(keys))
            for ikey, key in enumerate(keys[self.rank]):
                sta = ds.waveforms[key]
                #print 'Loading key number %d: %s'%(ikey, key)
                for tag in sta.list():

                    result = decode_tag(tag)
                    if (result):
                        network, station, location, channel, tr_st, tr_et = result
                        data.append([network, station, location, channel, tr_st, tr_et, tag])
                    # end if
                # end for
            # end for

            if(balance_by_tags):
                data = self.comm.allgather(data)
                data = [item for sublist in data for item in sublist]
                data = split_list(data, self.nproc)
                data = data[self.rank]
            # end if

            ki = 0
            tags_list = []
            t = tree()
            for row in data:
                network, station, location, channel, tr_st, tr_et, tag = row

                if (type(t[network][station][location][channel]) == defaultdict):
                    t[network][station][location][channel] = index.Index()
                # end if

                t[network][station][location][channel].insert(ki, (tr_st, 1, tr_et, 1))

                tags_list.append(tag)
                ki += 1
            # end for

            self.tree_list[ids] = t
            self.asdf_tags_list[ids] = tags_list

            tagsCount += len(data)
        # end for
        print 'Created sparse indices on rank %d for %d waveforms' % (self.rank, tagsCount)

        self.extract_time_ranges()
    # end func

    def get_stations(self, starttime, endtime, network=None, station=None, location=None, channel=None):
        starttime = UTCDateTime(starttime).timestamp
        endtime = UTCDateTime(endtime).timestamp
        # create a list of asdf datasets that may contain queried data
        dslistIndices = []
        for i, (stime, etime) in enumerate(zip(self.waveform_start_time_list,
                                               self.waveform_end_time_list)):
            if(stime==None or etime==None): continue
            if(FederatedASDFDataSet.hasOverlap(starttime, endtime, stime, etime)):
                dslistIndices.append(i)
            # end if
        # end for

        results = []
        for i in dslistIndices:
            #print 'Accessing file: %s'%(self.asdf_file_names[i])

            ds = self.asdf_datasets[i]

            val = self.tree_list[i]
            if(val): # has json db
                for (nk, nv) in val.iteritems():
                    if(network):
                        nk = network
                        nv = val[nk]
                    # end if
                    for (sk, sv) in nv.iteritems():
                        if (station):
                            sk = station
                            sv = nv[sk]
                        # end if
                        for (lk, lv) in sv.iteritems():
                            if (location):
                                lk = location
                                lv = sv[lk]
                            # end if
                            for (ck, cv) in lv.iteritems():
                                if (channel):
                                    ck = channel
                                    cv = lv[ck]
                                # end if

                                if (type(cv) == index.Index):
                                    tag_indices = list(cv.intersection((starttime, 1,
                                                                        endtime, 1)))

                                    if(len(tag_indices)):
                                        rv = [nk, sk, lk, ck,
                                              self.asdf_station_coordinates[i]['%s.%s'%(nk, sk)]['longitude'],
                                              self.asdf_station_coordinates[i]['%s.%s'%(nk, sk)]['latitude']]
                                        results.append(rv)
                                # end if
                                if(channel): break
                            # end for
                            if(location): break
                        # end for
                        if(station): break
                    # end for
                    if(network): break
                # end for
            else: # no json db
                '''
                TODO:
                Gather results based on ifilter interface of pyasdf -- it's quite slow btw.                
                '''
                pass
            # end if
        # end for

        return results
    # end func

    def get_waveforms(self, network, station, location, channel, starttime,
                      endtime, tag, automerge=False):

        starttime = UTCDateTime(starttime).timestamp
        endtime = UTCDateTime(endtime).timestamp
        # create a list of asdf datasets that may contain queried data
        dslistIndices = []

        for i, (stime, etime) in enumerate(zip(self.waveform_start_time_list,
                                               self.waveform_end_time_list)):
            if(stime==None or etime==None): continue
            if(FederatedASDFDataSet.hasOverlap(starttime, endtime, stime, etime)):
                dslistIndices.append(i)
            # end for
        # end for

        s = Stream()
        for i in dslistIndices:
            #print 'Accessing file: %s'%(self.asdf_file_names[i])

            ds = self.asdf_datasets[i]

            cs = Stream()
            if(self.tree_list[i]):
                val = self.tree_list[i][network][station][location][channel]
                if(type(val) == index.Index):
                    tag_indices = list(val.intersection((starttime, 1,
                                                         endtime, 1)))
                    station_data = ds.waveforms['%s.%s'%(network, station)]
                    for ti in tag_indices:
                        try:
                            cs += station_data[self.asdf_tags_list[i][ti]]
                        except:
                            pass
                        # end try
                    # end for
                # end if

                if(automerge):
                    try:
                        cs.merge(method=-1)
                    except:
                        pass
                    # end try
                # end if
            else:
                cs = ds.get_waveforms(network, station, location, channel, starttime,
                                      endtime, tag, automerge)
            # end if

            # Trim traces
            for t in cs:
                t.trim(starttime=UTCDateTime(starttime),
                       endtime=UTCDateTime(endtime))
                s += t
            # end for
        # end for

        return s
    # end func

    def local_net_sta_list(self):
        for it, val in enumerate(self.tree_list):
            if(val):
                for (nk, nv) in val.iteritems():
                    for (sk, sv) in nv.iteritems():
                        yield nk, sk, \
                              UTCDateTime(self.net_sta_start_time[nk][sk]), \
                              UTCDateTime(self.net_sta_end_time[nk][sk])
                    # end for
                # end for
            # end if
        # end for
    # end func

    def cleanup(self):
        for i, ds in enumerate(self.asdf_datasets):
            if(self.logger): self.logger.info('Closing ASDF file %s..'%(self.asdf_file_names[i]))
            del ds
        # end for
    # end func
# end class

if __name__=="__main__":
    fn = os.path.join('/tmp', 'test.log')
    logger = setup_logger('main', fn)

    fds = FederatedASDFDataSet('/g/data/ha3/rakib/tmp/a.txt', logger)
    #s = fds.get_waveforms('AU', 'QIS', '*', 'BHZ',
    #                      '2010-06-01T00:00:00', '2010-06-01T00:06:00',
    #                      'raw_recording', automerge=False)
    #print s
