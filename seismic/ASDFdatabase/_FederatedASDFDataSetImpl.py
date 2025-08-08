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
from ordered_set import OrderedSet as set
import numpy as np

from obspy.core import Stream, UTCDateTime
import pyasdf
from pyasdf.exceptions import ASDFValueError
from collections import defaultdict
import sqlite3
import hashlib
from functools import partial
from seismic.ASDFdatabase.utils import MIN_DATE, MAX_DATE, cleanse_inventory, \
    InventoryAggregator, get_file_signature
from seismic.misc import split_list, setup_logger
import pickle as cPickle
import pandas as pd
from rtree import index
import traceback
import gc

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
    def __init__(self, asdf_source, force_reindex=False, logger=None,
                 single_item_read_limit_in_mb=1024,
                 single_threaded_access=True):
        """
        :param asdf_source: path to a text file containing a list of ASDF files:
               Entries can be commented out with '#'
        :param force_reindex: Force reindex even if a preexisting db file is found
        :param logger: logger instance
        :param single_item_read_limit_in_mb: buffer size for Obspy reads
        :param single_threaded_access: By default, data are read via unthreaded MPI-processes.
               This can be relaxed for threaded GUI applications, though data access will still
               remain single-threaded.
        """

        self.comm = MPI.COMM_WORLD
        self.nproc = self.comm.Get_size()
        self.rank = self.comm.Get_rank()

        self.logger = logger
        self.single_threaded_access = single_threaded_access
        self.asdf_source = None
        self.asdf_file_names = []
        self.history_fn = None
        self.previous_db_fn = None
        self.asdf_station_coordinates = []
        self._unique_coordinates = defaultdict(list)

        if isinstance(asdf_source, str):
            self.asdf_source = asdf_source
            self.source_sha1 = hashlib.sha1(open(self.asdf_source).read().encode('utf-8')).hexdigest()
            self.db_fn = os.path.join(os.path.dirname(self.asdf_source), self.source_sha1 + '.db')
            self.history_fn = os.path.join(os.path.dirname(self.asdf_source), '.fasdf_history')

            fileContents = list(filter(len, open(self.asdf_source).read().splitlines()))

            # collate file names
            for i in range(len(fileContents)):
                if(fileContents[i][0]=='#'): continue # filter commented lines

                fn = fileContents[i].strip(' \t\n\r\n')
                if(os.path.exists(fn)):
                    self.asdf_file_names.append(os.path.abspath(fn))
                else:
                    print("Warning: file {} not found. Moving along..".format(fn))
                # end if
            # end for
            self.asdf_file_names = list(set(self.asdf_file_names)) # drop duplicates if present
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

        # Remove preexisting db if force_reindex is True
        if (force_reindex):
            if(self.rank == 0):
                if(os.path.exists(self.db_fn)):
                    os.remove(self.db_fn)
                # end if
            # end if
        # end if

        if(self.rank == 0):
            # retrieve an earlier version of the database, if available
            self.previous_db_fn = self._get_previous_db()
        # end if

        self.comm.Barrier()

        # Create database
        self.conn = None
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

    def _update_history(self):
        fh = open(self.history_fn, 'a+')
        fh.write('{}\n'.format(os.path.abspath(self.db_fn)))
        fh.close()
    # end func

    def _get_previous_db(self):
        result = None
        if(os.path.exists(self.history_fn)):
            fh = open(self.history_fn, 'r')
            lines = fh.readlines()
            if(len(lines) > 0):
                fn = lines[-1].strip() if len(lines[-1]) > 0 else None

                if(fn is not None and os.path.exists(fn)):
                    result = fn
                # end if
            # end if
        # end if
        return result
    # end func

    def _copy_from_previous_db(self, table_name: str, new_ds_id: int):
        """
        Copy entries from an earlier database if the associated asdf file has not changed
        @param table_name:
        @param new_ds_id:
        @return: boolean success/failure
        """

        rval = False
        cur = self.conn.cursor()

        # Attach source
        cur.execute("attach database ? as src", (self.previous_db_fn,))

        # Get ds_id from previous database where file-signature matches with that of given ds_id
        cur.execute(f"select old_ds.ds_id from src.ds as old_ds, ds as new_ds where \
                    old_ds.abs_path='{self.asdf_file_names[new_ds_id]}' and \
                    old_ds.abs_path=new_ds.abs_path and \
                    old_ds.st_size=new_ds.st_size and old_ds.st_mtime=new_ds.st_mtime \
                    and old_ds.st_ctime=new_ds.st_ctime and old_ds.st_ino=new_ds.st_ino \
                    and old_ds.st_dev=new_ds.st_dev")
        r = cur.fetchall()

        if(len(r) == 1):
            # found a matching entry in the previous database where the file signature matches
            old_ds_id = r[0][0]
            # Get column names from target table
            cur.execute(f"pragma table_info({table_name})")
            columns = [row[1] for row in cur.fetchall()]

            # all columns remain the same except ds_id, which is replaced by new_ds_id
            column_mods = {'ds_id': str(new_ds_id)}
            select_exprs = [
                column_mods[col] if col in column_mods else col
                for col in columns
            ]

            print(f"Copying entries in table '{table_name}' for '{self.asdf_file_names[new_ds_id]}'"
                  f" from earlier database ({self.previous_db_fn})")
            query = f"""
            insert into {table_name} ({', '.join(columns)})
            select {', '.join(select_exprs)}
            from src.{table_name} as src_table where src_table.ds_id={old_ds_id}
            """
            cur.execute(query)

            rval = True
        # end if

        self.conn.commit()
        cur.execute(f"detach database src")

        return rval
    # end func

    def create_database(self):
        def decode_tag(tag, type='raw_recording'):
            """
            Tags are expected in the form: {NET}.{STA}.{LOC}.{CHA}__{ST}__{ET}__{TAG}, where
            ST and ET are expected as YYYY-MM-DDThh:mm:ss.s.

            @param tag: tag
            @param type: str
            @return: network, station, location, channel, starttime, endtime
            """
            if (type not in tag): return None
            try:
                nslc, st, et, _ = tag.split('__')
                nc, sc, lc, cc = nslc.split('.')

                starttime = UTCDateTime(st).timestamp
                endtime  = UTCDateTime(et).timestamp

                if((endtime - starttime) == 0):
                    # Filtering out zero-length traces
                    return None
                else:
                    return nc, sc, lc, cc, starttime, endtime
                # end if
            except Exception:
                if self.logger:
                    self.logger.warning("Failed to decode tag {}".format(tag))
                return None
            # end try
        # end func

        dbFound = os.path.exists(self.db_fn)
        self.comm.Barrier()

        if(dbFound):
            print('Found database: %s'%(self.db_fn))
            self.conn = sqlite3.connect(self.db_fn,
                                        check_same_thread=self.single_threaded_access)
        else:
            if(self.rank==0):
                ia = InventoryAggregator()

                self.conn = sqlite3.connect(self.db_fn,
                                            check_same_thread=self.single_threaded_access)
                self.conn.execute('create table ds(ds_id smallint, abs_path text, '
                                  'st_size UNSIGNED BIG INT, st_mtime double, st_ctime double, '
                                  'st_ino UNSIGNED BIG INT, st_dev UNSIGNED BIG INT)')
                self.conn.execute('create table wtag(ds_id smallint, net varchar(6), sta varchar(6), loc varchar(6), '
                                  'cha varchar(6), st double, et double, tag text)')
                self.conn.execute('create table meta(ds_id smallint, net varchar(6), sta varchar(6), lon double, '
                                  'lat double, elev_m double)')
                self.conn.execute('create table masterinv(inv blob)')

                metadatalist = []
                for ids, ds in enumerate(self.asdf_datasets):
                    sig = get_file_signature(self.asdf_file_names[ids])
                    self.conn.execute('insert into ds(ds_id, abs_path, st_size, st_mtime, st_ctime, st_ino, st_dev) '
                                      'values(?, ?, ?, ?, ?, ?, ?)',
                                      [ids, sig['abs_path'], sig['st_size'], sig['st_mtime'],
                                       sig['st_ctime'], sig['st_ino'], sig['st_dev']])

                    coords_dict = ds.get_all_coordinates()

                    # report any missing metadata
                    wsta = set(list(ds.waveforms.list()))
                    msta = set(list(coords_dict.keys()))
                    if (len(wsta) != len(msta)):
                        missing = wsta - msta
                        print('WARNING: {} stations with missing metadata found in {}..'.\
                              format(len(missing), self.asdf_file_names[ids]))
                    # end if

                    # aggregate inventories
                    for k in coords_dict.keys():
                        inv = cleanse_inventory(ds.waveforms[k].StationXML)
                        ia.append(inv)
                    # end for

                    if(self.previous_db_fn is not None and \
                       self._copy_from_previous_db('meta', ids)):
                        # copied entries from previous database for ASDF files that
                        # have not changed
                        pass
                    else:
                        # failed to copy required entries from an earlier database
                        for k in coords_dict.keys():
                            # we keep coordinates from all ASDF files to be able to track
                            # potential discrepancies
                            lon = coords_dict[k]['longitude']
                            lat = coords_dict[k]['latitude']
                            elev_m = coords_dict[k]['elevation_in_m']
                            nc, sc = k.split('.')
                            metadatalist.append([ids, nc, sc, lon, lat, elev_m])
                        # end for
                    # end if
                # end for

                masterinv = ia.summarize()
                if(len(metadatalist) > 0):
                    self.conn.executemany('insert into meta(ds_id, net, sta, lon, lat, elev_m) values '
                                          '(?, ?, ?, ?, ?, ?)', metadatalist)
                # end if
                self.conn.execute('insert into masterinv(inv) values(?)',
                                  [cPickle.dumps(masterinv, cPickle.HIGHEST_PROTOCOL)])

                # clean up memory bloat caused by aggregated inventory
                del masterinv
                gc.collect()
                self.conn.commit()
                self.conn.close()
            # end if
            self.comm.Barrier()

            tagsCount = 0
            for ids, ds in enumerate(self.asdf_datasets):
                has_copied_entries = False
                if(self.rank==0):
                    if(self.previous_db_fn is not None):
                        self.conn = sqlite3.connect(self.db_fn,
                                                    check_same_thread=self.single_threaded_access)
                        has_copied_entries = self._copy_from_previous_db('wtag', ids)
                        self.conn.close()
                    # end if

                    if(not has_copied_entries):
                        print('Indexing %s..' % (os.path.basename(self.asdf_file_names[ids])))
                    # end if
                # end if
                has_copied_entries = self.comm.bcast(has_copied_entries, root=0)
                self.comm.Barrier()

                if(has_copied_entries): continue

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
                            self.conn = sqlite3.connect(self.db_fn,
                                                        check_same_thread=self.single_threaded_access)
                            self.conn.executemany('insert into wtag(ds_id, net, sta, loc, cha, st, et, tag) values '
                                                  '(?, ?, ?, ?, ?, ?, ?, ?)', data)
                            print('\tInserted %d entries on rank %d'%(len(data),
                                                                      self.rank))
                            tagsCount += len(data)
                            self.conn.commit()
                            self.conn.close()
                        # end if
                    # end if

                    self.comm.Barrier()
                # end for
            # end for

            if(self.rank==0):
                self.conn = sqlite3.connect(self.db_fn,
                                            check_same_thread=self.single_threaded_access)
                print('Creating table indices..')
                self.conn.execute('create index all_wtag_index on wtag(ds_id, net, sta, loc, cha, st, et)')
                self.conn.execute('create index all_meta_index on meta(ds_id, net, sta)')
                self.conn.execute('create index fast_wtag_index on wtag(net, sta, loc, cha, st, et)')
                self.conn.execute('create index fast_meta_index on meta(net, sta)')
                self.conn.commit()

                print('Creating convenience table with start-/end-times..')
                self.conn.execute('create table nslc as select net, sta, loc, cha, min(st) as st, max(et) as et '
                                  'from wtag group by net, sta, loc, cha order by net, sta, loc, cha')
                self.conn.execute('create index all_nslc_index on nslc(net, sta, loc, cha, st, et)')
                self.conn.commit()

                print('Creating convenience table containing total recording durations in seconds..')
                self.conn.execute('create table recording_time as select net, sta, loc, cha, sum(et-st) '
                                  'as duration_seconds from wtag group by net, sta, loc, cha '
                                  'order by net, sta, loc, cha;')
                self.conn.execute('create index all_recording_time_index on '
                                  'recording_time(net, sta, loc, cha, duration_seconds)')

                print('Creating convenience table containing timespans of continuous recordings')
                # use sqlite windowing to generate records of contiguous blocks of recordings where gaps
                # less than a day are ignored to keep the final row-count reasonable
                self.conn.execute("""create table coverage as WITH ordered AS (
                            SELECT 
                                net, sta, loc, cha, st, et,
                                LAG(et) OVER (PARTITION BY net, sta, loc, cha ORDER BY st, et) AS prev_et
                            FROM wtag
                        ),
                        segment_marks AS (
                            SELECT 
                                net, sta, loc, cha, st, et,
                                CASE 
                                    WHEN prev_et IS NULL OR st - prev_et >= 86400 THEN 1 
                                    ELSE 0 
                                END AS new_segment
                            FROM ordered
                        ),
                        segmented AS (
                            SELECT 
                                net, sta, loc, cha, st, et,
                                SUM(new_segment) OVER (PARTITION BY net, sta, loc, cha ORDER BY st, et) AS segment_id
                            FROM segment_marks
                        )
                        SELECT 
                            s.net, s.sta, s.loc, s.cha,  
                            MIN(st) AS block_st,
                            MAX(et) AS block_et
                        FROM segmented as s
                        GROUP BY s.net, s.sta, s.loc, s.cha, segment_id 
                        ORDER BY s.net, s.sta, s.loc, s.cha, block_st;""")
                self.conn.execute('create index all_coverage_index on coverage '
                                  '(net, sta, loc, cha, block_st, block_et)')

                self.conn.close()
                self._update_history() # update the history file with latest db_fn
                print('Done..')
            # end if
            self.comm.Barrier()
            self.conn = sqlite3.connect(self.db_fn,
                                        check_same_thread=self.single_threaded_access)
        # end if

        # Load metadata
        rows = self.conn.execute('select * from meta').fetchall()
        for row in rows:
            ds_id, net, sta, lon, lat, elev_m = row
            self.asdf_station_coordinates[ds_id]['%s.%s' % (net.strip(), sta.strip())] = [lon, lat, elev_m]
        # end for

        # Populate unique coordinates dict
        for ds_dict in self.asdf_station_coordinates:
            for key in list(ds_dict.keys()):
                lon, lat, _ = ds_dict[key]
                self._unique_coordinates[key] = [lon, lat]
            # end for
        # end for

        # Load master inventory
        row = self.conn.execute('select * from masterinv').fetchall()
        self.masterinv = cPickle.loads(row[0][0])
    # end func

    def get_recording_timespan(self, network, station=None, location=None, channel=None):
        query = "select min(st), max(et) from nslc where net='%s' " % (network)

        if (station is not None):
            query += "and sta='%s' " % (station)
        if (location is not None):
            query += "and loc='%s' " % (location)
        if (channel is not None):
            query += "and cha='%s' " % (channel)

        row = self.conn.execute(query).fetchall()[0]

        min = MAX_DATE
        max = MIN_DATE

        if (len(row)):
            if (row[0] is not None): min = UTCDateTime(row[0])
            if (row[1] is not None): max = UTCDateTime(row[1])
        # end if

        return min, max
    # end func

    def get_all_recording_timespans(self):
        query = "select net, sta, loc, cha, st, et from nslc"
        rows = self.conn.execute(query).fetchall()

        fields = {'names': ['net', 'sta', 'loc', 'cha', 'min_st', 'max_et'],
                  'formats': ['U10', 'U10', 'U10', 'U10', 'f8', 'f8']}
        result = np.zeros(len(rows), dtype=fields)

        for i, row in enumerate(rows): result[i] = row

        return result
    # end if

    def get_stations(self, starttime, endtime, network=None, station=None, location=None, channel=None):
        starttime = UTCDateTime(starttime).timestamp
        endtime = UTCDateTime(endtime).timestamp

        query = 'select ds_id, net, sta, loc, cha from wtag where '
        if (network is not None): query += " net='%s' "%(network)
        if (station is not None):
            if(network is not None): query += "and sta='%s' "%(station)
            else: query += "sta='%s' "%(station)
        if (location is not None):
            if((network is not None) or
               (station is not None)): query += "and loc='%s' "%(location)
            else: query += "loc='%s' "%(location)
        if (channel is not None):
            if((network is not None) or
               (station is not None) or
               (location is not None)): query += "and cha='%s' "%(channel)
            else: query += "cha='%s' "%(channel)
        if ((network is not None) or
            (station is not None) or
            (location is not None) or
            (channel is not None)): query += ' and '
        query += ' et>=%f and st<=%f' \
                 % (starttime, endtime)
        query += ' group by net, sta, loc, cha order by net, sta, loc, cha'

        rows = self.conn.execute(query).fetchall()
        results = set()
        for row in rows:
            ds_id, net, sta, loc, cha = row

            # [net, sta, loc, cha, lon, lat, elev_m]
            rv = (net, sta, loc, cha, *self.asdf_station_coordinates[ds_id]['%s.%s' % (net, sta)])
            results.add(rv)
        # end for

        return list(results)
    # end func

    def get_waveform_count(self, network, station, location, channel, starttime, endtime):

        starttime = UTCDateTime(starttime).timestamp
        endtime = UTCDateTime(endtime).timestamp

        query = "select count(*) from wtag where net='%s' and sta='%s' and loc='%s' and cha='%s' " \
                %(network, station, location, channel) + \
                "and et>=%f and st<=%f" \
                 % (starttime, endtime)

        num_traces = self.conn.execute(query).fetchall()[0][0]

        return num_traces
    # end func

    def get_waveforms(self, network, station, location, channel, starttime,
                      endtime, trace_count_threshold=200, nearest_sample=True):

        starttime = UTCDateTime(starttime)
        endtime = UTCDateTime(endtime)

        query = "select * from wtag where net='%s' and sta='%s' and loc='%s' and cha='%s' " \
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

        # Trim stream
        s.trim(starttime=starttime, endtime=endtime, nearest_sample=nearest_sample)

        # apply corrections if available
        if(self.corrections_enabled):
            s = self._apply_correction(s)
        # end if

        return s
    # end func

    def get_location_codes(self, network, station, starttime=None, endtime=None):
        st, et = self.get_recording_timespan(network, station)

        if(starttime):
            starttime = UTCDateTime(starttime)
            if(starttime > st): st = starttime
        # end if

        if(endtime):
            endtime = UTCDateTime(endtime)
            if(endtime < et): et = endtime
        # end if

        rows = self.get_stations(st, et, network=network, station=station)
        uniqueLocCodes = set()
        for row in rows:
            uniqueLocCodes.add(row[2])
        # end for

        return sorted(list(uniqueLocCodes))
    # end func

    def stations_iterator(self, network_list=[], station_list=[]):
        workload = None
        if(self.rank==0):
            workload = []
            for i in np.arange(self.nproc):
                workload.append(defaultdict(partial(defaultdict, list)))
            # end for

            nets = self.conn.execute('select distinct net from wtag').fetchall()
            if(len(network_list)): # filter networks
                nets = [net for net in nets if net[0] in network_list]
            # end if

            for net in nets:
                net = net[0]
                stas = self.conn.execute("select distinct sta from wtag where net='%s'"%(net)).fetchall()

                if (len(station_list)):  # filter stations
                    stas = [sta for sta in stas if sta[0] in station_list]
                # end if

                for sta in stas:
                    sta = sta[0]

                    # trace-count, min(st), max(et)
                    attribs = self.conn.execute("select count(st), min(st), max(et) from wtag where net='%s' and sta='%s'"
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
                  channel=None, starttime=None, endtime=None,
                  min_gap_length=86400):

        if(starttime is not None): starttime = UTCDateTime(starttime).timestamp
        if(endtime is not None): endtime= UTCDateTime(endtime).timestamp

        clause_added = 0
        query = 'select net, sta, loc, cha, st, et from wtag '
        if (network or station or location or channel or (starttime and endtime)): query += " where "

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

        if (starttime):
            if (clause_added):
                query += ' and st>={} '.format(starttime)
            else:
                query += ' st>={} '.format(starttime)
            clause_added += 1
        # end if

        if (endtime):
            if (clause_added):
                query += ' and et<={}'.format(endtime)
            else:
                query += ' et<={} '.format(endtime)
            clause_added += 1
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

    def get_recording_duration(self, network=None, station=None, location=None, channel=None,
                               starttime=None, endtime=None, cumulative=False):

        if(starttime is not None): starttime = UTCDateTime(starttime).timestamp
        if(endtime is not None): endtime= UTCDateTime(endtime).timestamp

        clause_added = 0
        query = """
            select net, sta, loc, cha, """

        if(starttime is not None and endtime is not None):
            if(cumulative):
                query += """
                    max(block_st, {}), min(block_et, {}), 
                    sum(
                        max(0, 
                            min(block_et, {}) - max(block_st, {})
                        )
                    ) as duration """.format(starttime, endtime, endtime, starttime)
            else:
                query += """ max(block_st, {}), min(block_et, {})
                        """.format(starttime, endtime)
            # end if
        else:
            if(cumulative):
                query += " min(block_st), max(block_et), sum (block_et - block_st) as duration "
            else:
                query += " block_st, block_et "
            # end if
        # end if

        query += " from coverage "

        if (network or station or location or channel or (starttime and endtime)): query += " where "

        if(starttime is not None and endtime is not None):
            query += """
              block_et > {}
              and block_st < {}
            """.format(starttime, endtime)
            clause_added += 1
        # end if

        if(network is not None):
            if(clause_added): query += "and net='{}' ".format(network)
            else: query += "net='{}' ".format(network)
            clause_added += 1
        # end if
        if(station is not None):
            if (clause_added): query += " and sta='{}' ".format(station)
            else: query += " sta='{}' ".format(station)
            clause_added += 1
        # end if
        if(location is not None):
            if(clause_added): query += " and loc='{}' ".format(location)
            else: query += " loc='{}' ".format(location)
            clause_added += 1
        # end if
        if(channel is not None):
            if(clause_added): query += " and cha='{}' ".format(channel)
            else: query += " cha='{}' ".format(channel)
            clause_added += 1
        # end if

        if(cumulative):
            query += """ 
                group by net, sta, loc, cha
                order by net, sta, loc, cha;
                """
        else:
            query += """ 
                group by net, sta, loc, cha, block_st, block_et
                order by net, sta, loc, cha, block_st, block_et;
                """
        # end if

        print('\n{}\n'.format(query))
        
        rows = self.conn.execute(query).fetchall()

        array_dtype = None
        if(cumulative):
            array_dtype = [('net', 'U10'), ('sta', 'U10'),
                           ('loc', 'U10'), ('cha', 'U10'),
                           ('min_st', 'float'), ('max_et', 'float'),
                           ('duration_seconds', 'float')]
        else:
            array_dtype = [('net', 'U10'), ('sta', 'U10'),
                           ('loc', 'U10'), ('cha', 'U10'),
                           ('block_st', 'float'), ('block_et', 'float')]
        # end if

        result = np.array(rows, dtype=array_dtype)

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
