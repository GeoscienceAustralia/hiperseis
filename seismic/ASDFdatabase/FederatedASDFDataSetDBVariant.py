"""
Description:
    Class for providing fast access to data contained within a set of ASDF files
    A reusable sqlite database is created for fast access to waveform data
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
import sqlite3
import psutil
import hashlib
from functools import partial

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

def split_list(lst, npartitions):
    k, m = divmod(len(lst), npartitions)
    return [lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(npartitions)]
# end func

def split_list_by_timespan(l, n):
    lmin = np.min(l[:, 1])
    lmax = np.max(l[:, 1])

    span = float(lmax - lmin) / n

    result = [[] for i in range(n)]
    s = lmin
    for i in range(n):
        e = s + span
        indices = np.where((l[:, 1] >= s) & (l[:, 1] <= e))
        result[i] = l[indices]
        s = e

        if(len(result[i]) == len(l)): break
    # end for

    return result
# end func

class FederatedASDFDataSetDBVariant():
    def __init__(self, asdf_source, logger=None):
        """
        :param asdf_source: path to a text file containing a list of ASDF files:
               Entries can be commented out with '#'
        :param logger: logger instance
        """

        self.comm = MPI.COMM_WORLD
        self.nproc = self.comm.Get_size()
        self.rank = self.comm.Get_rank()

        self.logger = logger
        self.asdf_source = None
        self.asdf_file_names = []
        self.asdf_station_coordinates = []

        if(type(asdf_source)==str):
            self.asdf_source = asdf_source
            self.source_sha1 = hashlib.sha1(open(self.asdf_source).read()).hexdigest()
            fileContents = filter(len, open(self.asdf_source).read().splitlines())

            # collate file names
            for i in range(len(fileContents)):
                if(fileContents[i][0]=='#'): continue # filter commented lines

                fn = fileContents[i].strip(' \t\n\r\n')
                self.asdf_file_names.append(fn)
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
                self.asdf_station_coordinates.append(defaultdict(list))
            else:
                raise NameError('File not found: %s..'%fn)
            # end if
        # end func

        # Create database
        self.conn = None
        self.db_fn = os.path.join(os.path.dirname(self.asdf_source),  self.source_sha1 + '.db')
        self.create_database()

        atexit.register(self.cleanup) # needed for closing asdf files at exit
    # end func

    def create_database(self):
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

        dbFound = os.path.exists(self.db_fn)
        self.comm.Barrier()

        if(dbFound):
            print 'Found database: %s'%(self.db_fn)
            self.conn = sqlite3.connect(self.db_fn)
        else:
            if(self.rank==0):
                self.conn = sqlite3.connect(self.db_fn)
                self.conn.execute('create table wdb(ds_id smallint, net varchar(6), sta varchar(6), loc varchar(6), '
                                  'cha varchar(6), st double, et double, tag text)')
                self.conn.execute('create table netsta(ds_id smallint, net varchar(6), sta varchar(6), lon double, '
                                  'lat double)')
                metadatalist = []
                for ids, ds in enumerate(self.asdf_datasets):
                    coords_dict = ds.get_all_coordinates()
                    for k in coords_dict.keys():
                        lon = coords_dict[k]['longitude']
                        lat = coords_dict[k]['latitude']
                        nc, sc = k.split('.')
                        metadatalist.append([ids, nc, sc, lon, lat])
                    # end for
                # end for
                self.conn.executemany('insert into netsta(ds_id, net, sta, lon, lat) values '
                                      '(?, ?, ?, ?, ?)', metadatalist)
            # end if

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
                            data.append([ids, network, station, location, channel, tr_st, tr_et, tag])
                        # end if
                    # end for
                # end for

                data = self.comm.gather(data, root=0)
                if(self.rank==0):
                    data = [item for sublist in data for item in sublist]
                    self.conn.executemany('insert into wdb(ds_id, net, sta, loc, cha, st, et, tag) values '
                                          '(?, ?, ?, ?, ?, ?, ?, ?)', data)
                    print 'Inserted %d entries on rank %d'%(len(data), self.rank)
                    tagsCount += len(data)
                    # end if
            # end for

            if(self.rank==0):
                self.conn.execute('create index allindex on wdb(net, sta, loc, cha, st, et)')
                self.conn.execute('create index netstaindex on netsta(ds_id, net, sta)')
                print 'Created database on rank %d for %d waveforms (%5.2f MB)' % \
                      (self.rank, tagsCount, round(psutil.Process().memory_info().rss / 1024. / 1024., 2))

                if (self.rank == 0): self.conn.close()
            # end if
            self.comm.Barrier()
            self.conn = sqlite3.connect(self.db_fn)
        # end if

        # Load metadata
        rows = self.conn.execute('select * from netsta').fetchall()
        for row in rows:
            ds_id, net, sta, lon, lat = row
            self.asdf_station_coordinates[ds_id]['%s.%s' % (net.strip(), sta.strip())] = [lon, lat]
        # end for
    # end func

    def get_global_time_range(self, network, station, location=None, channel=None):
        query = "select min(st), max(et) from wdb where net='%s' and sta='%s' "%(network, station)

        if (location):
            query += "and loc='%s' "%(location)
        if (channel):
            query += "and cha='%s' "%(channel)

        row = self.conn.execute(query).fetchall()[0]

        min = UTCDateTime(row[0]) if row[0] else UTCDateTime(4102444800.0)
        max = UTCDateTime(row[1]) if row[1] else UTCDateTime(-2208988800.0)
        return min, max
    # end func

    def get_stations(self, starttime, endtime, network=None, station=None, location=None, channel=None):
        starttime = UTCDateTime(starttime).timestamp
        endtime = UTCDateTime(endtime).timestamp

        query = 'select * from wdb where '
        if (network): query += " net='%s' "%(network)
        if (station):
            if(network): query += "and sta='%s' "%(station)
            else: query += "sta='%s' "%(station)
        if (location):
            if(network or station): query += "and loc='%s' "%(location)
            else: query += "loc='%s' "%(location)
        if (channel):
            if(network or station or location): query += "and cha='%s' "%(channel)
            else: query += "cha='%s' "%(channel)
        if (network or station or location or channel): query += ' and '
        query += ' et>=%f and st<=%f' \
                 % (starttime, endtime)
        query += ' group by net, sta, loc, cha'

        rows = self.conn.execute(query).fetchall()
        results = []
        for row in rows:
            ds_id, net, sta, loc, cha, st, et, tag = row

            rv = [net, sta, loc, cha,
                  self.asdf_station_coordinates[ds_id]['%s.%s' % (net, sta)][0],
                  self.asdf_station_coordinates[ds_id]['%s.%s' % (net, sta)][1]]
            results.append(rv)
        # end for

        return results
    # end func

    def get_waveforms(self, network, station, location, channel, starttime,
                      endtime, automerge=False, trace_count_threshold=200):

        starttime = UTCDateTime(starttime).timestamp
        endtime = UTCDateTime(endtime).timestamp

        query = "select * from wdb where net='%s' and sta='%s' and loc='%s' and cha='%s' " \
                %(network, station, location, channel) + \
                "and et>=%f and st<=%f" \
                 % (starttime, endtime)

        rows = self.conn.execute(query).fetchall()
        s = Stream()

        if(len(rows) > trace_count_threshold): return s

        #print 'rank: %d net: %s sta: %s loc:%s cha:%s ntraces: %d'%(self.rank,
        #                                                           network, station,
        #                                                           location, channel, len(rows))
        #return s

        for row in rows:
            ds_id, net, sta, loc, cha, st, et, tag = row
            station_data = self.asdf_datasets[ds_id].waveforms['%s.%s'%(net, sta)]
            try:
                s += station_data[tag]
            except:
                pass
            # end try
        # end for

        if(automerge):
            try:
                s.merge(method=-1)
            except:
                pass
            # end try
        # end if

        # Trim traces
        for t in s:
            t.trim(starttime=UTCDateTime(starttime),
                   endtime=UTCDateTime(endtime))
        # end for

        return s
    # end func

    def local_net_sta_list(self):
        workload = None
        if(self.rank==0):
            workload = []
            for i in np.arange(self.nproc):
                workload.append(defaultdict(partial(defaultdict, list)))
            # end for

            nets = self.conn.execute('select distinct net from wdb').fetchall()
            for net in nets:
                net = net[0]
                stas = self.conn.execute("select distinct sta from wdb where net='%s'"%(net)).fetchall()
                for sta in stas:
                    sta = sta[0]

                    tbounds = self.conn.execute("select st, et from wdb where net='%s' and sta='%s' order by et"
                                                %(net, sta)).fetchall()

                    if(len(tbounds)==0): continue
                    #tbounds = np.array(tbounds)
                    #tbounds = split_list_by_timespan(tbounds, self.nproc)
                    tbounds = split_list(tbounds, self.nproc)

                    for iproc in np.arange(self.nproc):
                        if (len(tbounds[iproc])):
                            arr = np.array(tbounds[iproc])
                            workload[iproc][net][sta] = np.array([np.min(arr[:,0]),
                                                                  np.max(arr[:,1])])
                        # end for
                    # end for
                # end for
            # end for
        # end if

        workload = self.comm.scatter(workload, root=0)
        for (nk, nv) in workload.iteritems():
            for (sk, sv) in nv.iteritems():
                start_time = None
                end_time = None
                try:
                    start_time = UTCDateTime(workload[nk][sk][0])
                    end_time = UTCDateTime(workload[nk][sk][1])
                except:
                    continue
                # end try

                yield nk, sk, start_time, end_time
            # end for
        # end for
    # end func

    def cleanup(self):
        for i, ds in enumerate(self.asdf_datasets):
            if(self.logger): self.logger.info('Closing ASDF file %s..'%(self.asdf_file_names[i]))
            del ds
        # end for

        self.conn.close()
    # end func
# end class

if __name__=="__main__":
    fn = os.path.join('/tmp', 'test.log')
    logger = setup_logger('main', fn)

    fds = FederatedASDFDataSetDBVariant('/g/data/ha3/rakib/tmp/a.txt', logger)
    #s = fds.get_waveforms('AU', 'QIS', '*', 'BHZ',
    #                      '2010-06-01T00:00:00', '2010-06-01T00:06:00',
    #                      'raw_recording', automerge=False)
    #print s
