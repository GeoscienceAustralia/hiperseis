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
from future.utils import iteritems
from builtins import range

from mpi4py import MPI
import os
import glob
import atexit
import logging
from ordered_set import OrderedSet as set
import numpy as np

from obspy.core import Stream, UTCDateTime
import pyasdf
from pyasdf.exceptions import ASDFValueError
from collections import defaultdict
import sqlite3
import hashlib
from functools import partial
from seismic.ASDFdatabase.utils import MIN_DATE, MAX_DATE
import pickle as cPickle
import pandas as pd
from rtree import index
import traceback

logging.basicConfig()

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
    return [lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(npartitions)]
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

class _FederatedASDFDataSetImpl():
    def __init__(self, asdf_source, logger=None, single_item_read_limit_in_mb=1024):
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

        if isinstance(asdf_source, str):
            self.asdf_source = asdf_source
            self.source_sha1 = hashlib.sha1(open(self.asdf_source).read().encode('utf-8')).hexdigest()
            fileContents = list(filter(len, open(self.asdf_source).read().splitlines()))

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
                ds.single_item_read_limit_in_mb = single_item_read_limit_in_mb
                self.asdf_datasets.append(ds)
                self.asdf_station_coordinates.append(defaultdict(list))
            else:
                raise NameError('File not found: %s..'%fn)
            # end if
        # end func

        # Create database
        self.conn = None
        self.db_fn = os.path.join(os.path.dirname(self.asdf_source),  self.source_sha1 + '.db')
        self.masterinv = None
        self.create_database()
        self._load_corrections()

        atexit.register(self.cleanup) # needed for closing asdf files at exit
    # end func

    def _load_corrections(self):
        self.correction_files = []
        self.correction_map_tree = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        self.correction_map_bounds = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        self.correction_map_values = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

        # check to see if corrections are to be applied
        self.corrections_enabled = False
        if('GPS_CLOCK_CORRECTION' in os.environ.keys()):
            try:
                self.corrections_enabled = np.bool_(np.int_(os.environ['GPS_CLOCK_CORRECTION']))
            except Exception as e:
                print(str(e))
                assert 0, 'Invalid value for GPS_CLOCK_CORRECTION: {}. Must be 1 or 0. Aborting..'.format(os.environ['GPS_CLOCK_CORRECTION'])
            # end try
        # end if
        
        if(not self.corrections_enabled): return

        pattern = os.path.join(os.path.dirname(self.asdf_source), '.corrections/*.csv')
        fnames = glob.glob(pattern)

        if(len(fnames)): print('Loading corrections..')

        dtypes = {'net':str, 'sta':str, 'loc':str, 'comp':str, 'date':str,
                  'clock_correction':str}
        for fname in fnames:
            df = pd.read_csv(fname, delimiter=',', header=0, dtype=dtypes, na_filter=False)

            try:
                corr_count = 0
                for i in np.arange(len(df)):
                    net = df['net'][i]
                    sta = df['sta'][i]
                    loc = df['loc'][i]
                    corr = df['clock_correction'][i]

                    if(corr == 'NOXCOR'): continue
                    else: corr = float(df['clock_correction'][i])

                    st = UTCDateTime(df['date'][i]).timestamp
                    et = st + 24*3600

                    if(type(self.correction_map_tree[net][sta][loc]) != index.Index):
                        self.correction_map_tree[net][sta][loc] = index.Index()
                        self.correction_map_bounds[net][sta][loc] = []
                        self.correction_map_values[net][sta][loc] = []
                    # end if

                    self.correction_map_tree[net][sta][loc].insert(corr_count, (st, 1, et, 1))
                    self.correction_map_bounds[net][sta][loc].append([st, et])
                    self.correction_map_values[net][sta][loc].append(corr)
                    corr_count += 1
                # end for
            except Exception as e:
                print ('Warning: failed to read corrections file {} with error({}). '
                       'Continuing along..'.format(fname, traceback.format_exc()))
            #end try
        # end for
    #end func

    def _get_correction(self, net, sta, loc, st, et):
        tindex = self.correction_map_tree[net][sta][loc]
        if(type(tindex) != index.Index):
            return None
        else:
            epsilon = 1e-5
            indices = list(tindex.intersection((st.timestamp+epsilon, 1, et.timestamp-epsilon, 1)))

            if(len(indices)):
                if(len(indices) > 1):
                    print('Warning: multivalued corrections found ({}.{}.{}: {} - {}). '
                          'Ignoring and moving along..'.format(net, sta, loc, st, et))
                    return None
                # end if

                cst, cet = self.correction_map_bounds[net][sta][loc][indices[0]]
                a = np.fmax(st.timestamp, cst)
                b = np.fmin(et.timestamp, cet)

                if(a > b): # sanity check
                    raise ValueError('Error encountered in _get_correction (({}.{}.{}: {} - {})). '
                                     'Aborting..'.format(net, sta, loc, st, et))
                # end if

                # return overlap and correction
                return [UTCDateTime(a), UTCDateTime(b)], self.correction_map_values[net][sta][loc][indices[0]]
            else:
                return None
            # end if
        #end if
    #end func

    def _apply_correction(self, stream):

        def day_split(trc):
            rstream = Stream()
            daySeconds = 24 * 3600  # seconds

            st = trc.stats.starttime
            et = trc.stats.endtime
            dayAlignedStartTime = UTCDateTime(year=st.year, month=st.month, day=st.day)
            dayAlignedEndTime = dayAlignedStartTime + daySeconds

            ct = st
            trcCopy = trc.copy()
            while (ct < et):
                step = daySeconds

                if (ct + step > dayAlignedEndTime):
                    step = dayAlignedEndTime - ct
                # end if

                rstream += trcCopy.slice(ct, ct + step, nearest_sample=False)

                ct += step
                dayAlignedEndTime += daySeconds
            # wend

            return rstream
        # end func

        resultStream = Stream()
        for mtr in stream:
            dayStream = day_split(mtr)

            for tr in dayStream:
                net = tr.stats.network
                sta = tr.stats.station
                loc = tr.stats.location
                st  = tr.stats.starttime
                et  = tr.stats.endtime
            
                if(st == et):
                    resultStream.append(tr)
                    continue
                # end if
                
                result = self._get_correction(net, sta, loc, st, et)

                if(result):
                    ost, oet = result[0]
                    corr = result[1]

                    trCorrected = tr.copy().slice(ost, oet)
                    trCorrected.stats.starttime -= corr
                    currStream = Stream([tr.slice(st, ost),
                                         trCorrected,
                                         tr.slice(oet, et)])

                    # overlap resulting from correction applied is discarded (method=0)
                    currStream.merge()
                    if(len(currStream) == 1):
                        resultStream.append(currStream[0])
                    else:
                        raise ValueError('Merge error in _apply_correction')
                    # end if
                else:
                    resultStream.append(tr)
                # end if
            # end for
        # end for

        return resultStream
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
                endtime  = UTCDateTime(tokens[2]).timestamp

                return nc, sc, lc, cc, starttime, endtime
            except Exception:
                if self.logger:
                    self.logger.error("Failed to decode tag {}".format(tag))
                return None
            # end try
        # end func

        dbFound = os.path.exists(self.db_fn)
        self.comm.Barrier()

        if(dbFound):
            print(('Found database: %s'%(self.db_fn)))
            self.conn = sqlite3.connect(self.db_fn)
        else:
            if(self.rank==0):
                self.conn = sqlite3.connect(self.db_fn)
                self.conn.execute('create table wdb(ds_id smallint, net varchar(6), sta varchar(6), loc varchar(6), '
                                  'cha varchar(6), st double, et double, tag text)')
                self.conn.execute('create table netsta(ds_id smallint, net varchar(6), sta varchar(6), lon double, '
                                  'lat double, elev_m double)')
                self.conn.execute('create table masterinv(inv blob)')

                metadatalist = []
                masterinv = None
                for ids, ds in enumerate(self.asdf_datasets):
                    coords_dict = ds.get_all_coordinates()

                    for k in coords_dict.keys():
                        if(not masterinv):
                            masterinv = ds.waveforms[k].StationXML
                        else:
                            try:
                                masterinv += ds.waveforms[k].StationXML
                            except Exception as e:
                                print(e)
                            # end try
                        # end if
                    # end for

                    for k in list(coords_dict.keys()):
                        lon = coords_dict[k]['longitude']
                        lat = coords_dict[k]['latitude']
                        elev_m = coords_dict[k]['elevation_in_m']
                        nc, sc = k.split('.')
                        metadatalist.append([ids, nc, sc, lon, lat, elev_m])
                    # end for
                # end for
                self.conn.executemany('insert into netsta(ds_id, net, sta, lon, lat, elev_m) values '
                                      '(?, ?, ?, ?, ?, ?)', metadatalist)
                self.conn.execute('insert into masterinv(inv) values(?)',
                                  [cPickle.dumps(masterinv, cPickle.HIGHEST_PROTOCOL)])
                self.conn.commit()
                self.conn.close()
            # end if

            tagsCount = 0
            for ids, ds in enumerate(self.asdf_datasets):
                if(self.rank==0): print(('Indexing %s..' % (os.path.basename(self.asdf_file_names[ids]))))

                keys = list(ds.get_all_coordinates().keys())
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

                for irank in np.arange(self.nproc):
                    if(irank == self.rank):
                        if(len(data)):
                            self.conn = sqlite3.connect(self.db_fn)
                            self.conn.executemany('insert into wdb(ds_id, net, sta, loc, cha, st, et, tag) values '
                                                  '(?, ?, ?, ?, ?, ?, ?, ?)', data)
                            print(('\tInserted %d entries on rank %d'%(len(data),
                                                                       self.rank)))
                            tagsCount += len(data)
                            self.conn.commit()
                            self.conn.close()
                        # end if
                    # end if

                    self.comm.Barrier()
                # end for
            # end for

            if(self.rank==0):
                print('Creating table indices..')
                self.conn = sqlite3.connect(self.db_fn)
                self.conn.execute('create index allindex on wdb(net, sta, loc, cha, st, et)')
                self.conn.execute('create index netstaindex on netsta(ds_id, net, sta)')
                self.conn.commit()
                self.conn.close()
                print('Done..')
            # end if
            self.comm.Barrier()
            self.conn = sqlite3.connect(self.db_fn)
        # end if

        # Load metadata
        rows = self.conn.execute('select * from netsta').fetchall()
        for row in rows:
            ds_id, net, sta, lon, lat, elev_m = row
            self.asdf_station_coordinates[ds_id]['%s.%s' % (net.strip(), sta.strip())] = [lon, lat, elev_m]
        # end for

        # Load master inventory
        row = self.conn.execute('select * from masterinv').fetchall()
        self.masterinv = cPickle.loads(row[0][0])
    # end func

    def get_global_time_range(self, network, station, location=None, channel=None):
        query = "select min(st), max(et) from wdb where net='%s' and sta='%s' "%(network, station)

        if (location):
            query += "and loc='%s' "%(location)
        if (channel):
            query += "and cha='%s' "%(channel)

        row = self.conn.execute(query).fetchall()[0]

        min = MAX_DATE
        max = MIN_DATE

        if(len(row)):
            if(row[0] is not None): min = UTCDateTime(row[0])
            if(row[1] is not None): max = UTCDateTime(row[1])
        # end if

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
        results = set()
        for row in rows:
            ds_id, net, sta, loc, cha, st, et, tag = row

            # [net, sta, loc, cha, lon, lat, elev_m]
            rv = (net, sta, loc, cha, *self.asdf_station_coordinates[ds_id]['%s.%s' % (net, sta)])
            results.add(rv)
        # end for

        return list(results)
    # end func

    def get_waveform_count(self, network, station, location, channel, starttime, endtime):

        starttime = UTCDateTime(starttime).timestamp
        endtime = UTCDateTime(endtime).timestamp

        query = "select count(*) from wdb where net='%s' and sta='%s' and loc='%s' and cha='%s' " \
                %(network, station, location, channel) + \
                "and et>=%f and st<=%f" \
                 % (starttime, endtime)

        num_traces = self.conn.execute(query).fetchall()[0][0]

        return num_traces
    # end func

    def get_waveforms(self, network, station, location, channel, starttime,
                      endtime, trace_count_threshold=200):

        starttime = UTCDateTime(starttime)
        endtime = UTCDateTime(endtime)

        query = "select * from wdb where net='%s' and sta='%s' and loc='%s' and cha='%s' " \
                %(network, station, location, channel) + \
                "and et>=%f and st<=%f" \
                 % (starttime.timestamp, endtime.timestamp)

        #print(query)

        rows = self.conn.execute(query).fetchall()
        s = Stream()

        if(len(rows) > trace_count_threshold):
            print("Trace Count exceeds threshold", len(rows))
            return s

        for row in rows:
            ds_id, net, sta, loc, cha, st, et, tag = row
            station_data = self.asdf_datasets[ds_id].waveforms['%s.%s'%(net, sta)]

            '''
            Obspy currently reads all data for a given 'tag' and then trims them as needed. However,
            the read operation fails if data for a given 'tag' exceeds 'single_item_read_limit_in_mb',
            regardless of the timespan indicated by starttime and endtime, if provided. In such 
            instances, we retry reading the data by expanding the read-buffer in each attempt. The 
            current max retry attempts is set to 2 and 'single_item_read_limit_in_mb' is finally reset
            to its original value.
            '''
            numAttempts = 0
            while(1):
                try:
                    data_segment = station_data.get_item(tag, starttime, endtime)
                    s += data_segment
                    break
                except Exception as e:
                    if(isinstance(e, ASDFValueError)): # read failed due to the data buffer being too small
                        self.asdf_datasets[ds_id].single_item_read_limit_in_mb *= 2
                        numAttempts += 1
                        if self.logger:
                            self.logger.warning("Failed to get data between {} -- {} for {}.{} due to:\n{}. "
                                                "Retrying with expanded data buffer."
                                                .format(str(starttime), str(endtime), net, sta, str(e)))
                        # end if
                        if(numAttempts > 2):
                            self.logger.error("Failed to get data between {} -- {} for {}.{} with error:\n{}"
                                              .format(str(starttime), str(endtime), net, sta, str(e)))
                            break
                        # end if
                    else:
                        if self.logger:
                            self.logger.error("Failed to get data between {} -- {} for {}.{} with error:\n{}"
                                              .format(str(starttime), str(endtime), net, sta, str(e)))
                        # end if
                        break
                    # end if
                # end try
            # end while
            if (numAttempts>0):
                self.asdf_datasets[ds_id].single_item_read_limit_in_mb /= 2**numAttempts
            # end if
        # end for

        # Trim traces
        for t in s:
            t.trim(starttime=starttime,
                   endtime=endtime)
        # end for

        # apply corrections if available
        if(self.corrections_enabled):
            s = self._apply_correction(s)
        # end if

        return s
    # end func

    def stations_iterator(self, network_list=[], station_list=[]):
        workload = None
        if(self.rank==0):
            workload = []
            for i in np.arange(self.nproc):
                workload.append(defaultdict(partial(defaultdict, list)))
            # end for

            nets = self.conn.execute('select distinct net from wdb').fetchall()
            if(len(network_list)): # filter networks
                nets = [net for net in nets if net[0] in network_list]
            # end if

            for net in nets:
                net = net[0]
                stas = self.conn.execute("select distinct sta from wdb where net='%s'"%(net)).fetchall()

                if (len(station_list)):  # filter stations
                    stas = [sta for sta in stas if sta[0] in station_list]
                # end if

                for sta in stas:
                    sta = sta[0]

                    # trace-count, min(st), max(et)
                    attribs = self.conn.execute("select count(st), min(st), max(et) from wdb where net='%s' and sta='%s'"
                                                %(net, sta)).fetchall()

                    if(len(attribs)==0): continue
                    tcount, min_st, max_et = np.array(attribs).flatten()

                    # create start and end times for each rank
                    r = np.linspace(min_st, max_et, self.nproc + 1)
                    rank_spans = np.vstack([r[0:-1], r[1:]]).T

                    # reproducibly shuffle rank-spans to balance load across ranks
                    np.random.seed(int(tcount))
                    rank_spans = split_list(rank_spans, self.nproc)
                    np.random.shuffle(rank_spans)

                    for iproc in np.arange(self.nproc):
                        if (len(rank_spans[iproc])):
                            workload[iproc][net][sta] = rank_spans[iproc].flatten()
                        # end for
                    # end for
                # end for
            # end for
        # end if

        workload = self.comm.scatter(workload, root=0)
        for (nk, nv) in iteritems(workload):
            for (sk, sv) in iteritems(nv):
                start_time = None
                end_time = None
                try:
                    start_time = UTCDateTime(workload[nk][sk][0])
                    end_time = UTCDateTime(workload[nk][sk][1])
                except Exception:
                    if self.logger:
                        self.logger.warning("Failed to convert start and end times for keys {}, {}".format(nk, sk))
                    continue
                # end try

                yield nk, sk, start_time, end_time
            # end for
        # end for
    # end func

    def get_inventory(self, network=None, station=None):
        inv = self.masterinv.select(network=network, station=station)

        return inv
    # end func

    def find_gaps(self, network=None, station=None, location=None,
                  channel=None, start_date_ts=None, end_date_ts=None,
                  min_gap_length=86400):

        clause_added = 0
        query = 'select net, sta, loc, cha, st, et from wdb '
        if (network or station or location or channel or (start_date_ts and end_date_ts)): query += " where "

        if (network):
            query += ' net="{}" '.format(network)
            clause_added += 1
        # end if

        if (station):
            if (clause_added):
                query += ' and sta="{}" '.format(station)
            else:
                query += ' sta="{}" '.format(station)
            clause_added += 1
        # end if

        if (location):
            if (clause_added):
                query += ' and loc="{}" '.format(location)
            else:
                query += ' loc="{}" '.format(location)
            clause_added += 1
        # end if

        if (channel):
            if (clause_added):
                query += ' and cha="{}" '.format(channel)
            else:
                query += ' cha="{}" '.format(channel)
            clause_added += 1
        # end if

        if (start_date_ts and end_date_ts):
            if (clause_added):
                query += ' and st>={} and et<={}'.format(start_date_ts, end_date_ts)
            else:
                query += ' st>={} and et<={}'.format(start_date_ts, end_date_ts)
        # end if
        query += ' order by st, et'

        rows = self.conn.execute(query).fetchall()

        array_dtype = [('net', 'U10'), ('sta', 'U10'),
                       ('loc', 'U10'), ('cha', 'U10'),
                       ('st', 'float'), ('et', 'float')]
        rows = np.array(rows, dtype=array_dtype)

        # Process rows
        tree = lambda: defaultdict(tree)
        nested_dict = tree()
        for i in np.arange(rows.shape[0]):
            net = rows['net'][i]
            sta = rows['sta'][i]
            loc = rows['loc'][i]
            cha = rows['cha'][i]
            st = rows['st'][i]
            et = rows['et'][i]

            if (type(nested_dict[net][sta][loc][cha]) == defaultdict):
                nested_dict[net][sta][loc][cha] = []
            # end if

            nested_dict[net][sta][loc][cha].append([st, et])
        # end for

        result = []
        for net in nested_dict.keys():
            for sta in nested_dict[net].keys():
                for loc in nested_dict[net][sta].keys():
                    for cha in nested_dict[net][sta][loc].keys():
                        arr = nested_dict[net][sta][loc][cha]
                        if (len(arr)):
                            arr = np.array(arr)
                            st = arr[:, 0]
                            et = arr[:, 1]
                            assert np.allclose(np.array(sorted(st)), st), 'Start-times array not sorted!'
                            gaps = np.argwhere((st[1:] - et[:-1]) >= min_gap_length)

                            if (len(gaps)):
                                for i, idx in enumerate(gaps):
                                    idx = idx[0]

                                    result.append((net, sta, loc, cha, et[idx], st[idx + 1]))
                                # end for
                            # end if
                        # end if
                    # end for
                # end for
            # end for
        # end for
        result = np.array(result, dtype=array_dtype)

        return result
    # end func

    def cleanup(self):
        for i, ds in enumerate(self.asdf_datasets):
            # if self.logger:
            #     self.logger.info('Closing ASDF file %s..'%(self.asdf_file_names[i]))
            del ds
        # end for

        self.conn.close()
    # end func
# end class

if __name__=="__main__":
    fn = os.path.join('/tmp', 'test.log')
    logger = setup_logger('main', fn)

    # fds = _FederatedASDFDataSetImpl('/g/data/ha3/rakib/tmp/a.txt', logger)
    fds = _FederatedASDFDataSetImpl("/g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt", logger)
    s = fds.get_waveforms('AU', 'QIS', '', 'BHZ', '2015-06-01T00:00:00', '2015-06-02T00:06:00')  # cannot use wildcard *
    print(s)
    # select * from wdb where net='AU' and sta='QIS' and loc='' and cha='BHZ' and et>=1433116800.000000 and st<=1433203560.000000
    # 2 Trace(s) in Stream:
    # AU.QIS..BHZ | 2015-05-31T23:59:59.994500Z - 2015-06-01T23:59:59.994500Z | 40.0 Hz, 3456001 samples
    # AU.QIS..BHZ | 2015-06-01T23:59:55.769500Z - 2015-06-02T00:05:59.994500Z | 40.0 Hz, 14570 samples

    print(s)
