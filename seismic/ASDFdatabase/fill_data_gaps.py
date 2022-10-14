"""
Description:
    Script for filling gaps in a FederatedASDFDataSet archive

References:

CreationDate:   05/10/22
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     05/10/22   RH
"""

import numpy as np
from obspy import UTCDateTime
import click
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
from seismic.ASDFdatabase.utils import MseedIndex
from seismic.ASDFdatabase.utils import MIN_DATE, MAX_DATE
import os
from mpi4py import MPI
from obspy.core import Stream
import pyasdf
import shutil
import tempfile
from ordered_set import OrderedSet as set

class DataPool:
    def __init__(self, sources_txt_file:str, min_gap_length:float, scratch_path:str, mseed_pattern='*.mseed'):
        self.sources_txt_file = sources_txt_file
        self.min_gap_length = min_gap_length
        self.scratch_path = scratch_path
        self.mseed_pattern = mseed_pattern
        self.sources = []
        self.comm = MPI.COMM_WORLD
        self.nproc = self.comm.Get_size()
        self.rank = self.comm.Get_rank()
        self.asdf_source_fn = None
        self.tempdir = None

        lines = set(open(self.sources_txt_file).readlines())
        h5_lines = set([line for line in lines if ('h5' or 'H5') in line])
        mseed_lines = lines.difference(h5_lines)

        h5_lines = list(h5_lines)
        mseed_lines = list(mseed_lines)

        # Create temporary folder for storing index
        if (self.rank == 0):
            self.tempdir = tempfile.mkdtemp(dir=self.scratch_path)
            self.asdf_source_fn = os.path.join(self.tempdir, 'asdf_files.txt')

            open(self.asdf_source_fn, 'w').writelines(h5_lines)
        # end if
        self.asdf_source_fn = self.comm.bcast(self.asdf_source_fn, root=0)
        self.comm.barrier()

        self.sources.append(FederatedASDFDataSet(self.asdf_source_fn))
        for mseed_folder in mseed_lines:
            self.sources.append(MseedIndex(mseed_folder.strip(), self.mseed_pattern))
        # end for
    # end func

    def get_time_range(self, net, sta, loc, cha):
        result_min = MAX_DATE
        result_max = MIN_DATE

        for src in self.sources:
            min = max = None
            if(type(src)==MseedIndex): min, max = src.get_time_range(net, sta, loc, cha)
            else: min, max = src.get_global_time_range(net, sta, loc, cha)

            if(min < result_min): result_min = min
            if(max > result_max): result_max = max
        # end for

        return result_min, result_max
    # end fun

    def get_stations(self, st:UTCDateTime, et:UTCDateTime, net=None, sta=None, loc=None, cha=None):
        result = []
        for src in self.sources:
            rows = src.get_stations(st, et, net, sta, loc, cha)
            if(len(rows)):
                for row in rows: result.append(tuple(row[0:4]))
            # end if
        # end for
        return result
    # end func

    def get_waveforms(self, net, sta, loc, cha, st:UTCDateTime, et:UTCDateTime):
        assert self.rank == 0, 'This function is only accessible from Rank 0. Aborting..'

        def get_gaps(test_stream):
            if(len(test_stream) == 0): return [[net, sta, loc, cha, st, et]]
            _gaps = test_stream.get_gaps(min_gap=self.min_gap_length)

            min_st = np.min(np.array([tr.stats.starttime for tr in test_stream]))
            max_et = np.max(np.array([tr.stats.endtime for tr in test_stream]))

            if (st < min_st): _gaps.append([net, sta, loc, cha, st, min_st])
            if (et > max_et): _gaps.append([net, sta, loc, cha, max_et, et])

            return _gaps
        # end func

        wd = Stream([])
        gaps = []
        for i, src in enumerate(self.sources):
            if(i == 0):
                before = len(wd)
                wd += src.get_waveforms(net, sta, loc, cha, st, et)
                after = len(wd)

                #if(after>before):
                #    print('Found traces in {}'.format(i))
                #    print(wd)

                gaps = get_gaps(wd)
                #print("----1", gaps)
            elif(len(gaps)):
                before = len(wd)
                for gap in gaps:
                    wd += src.get_waveforms(net, sta, loc, cha, gap[4], gap[5])
                # end for
                after = len(wd)

                #if(after>before):
                #    print('Found traces in {}'.format(i))
                #    print(wd)

                gaps = get_gaps(wd)
                #print("----2", gaps)
            # end if
        # end for

        return wd
    # end func

    def waveform_iterator(self, net, sta, loc, cha, st:UTCDateTime, et:UTCDateTime, step=86400):
        wst = st
        wet = st
        while (wet < et):
            if (wst + step > et): wet = et
            else: wet = wst + step

            wd = self.get_waveforms(net, sta, loc, cha, wst, wet)
            yield wd

            wst += step
        # wend
    # end func

    def __del__(self):
        if(self.rank == 0): shutil.rmtree(self.tempdir)
    # end func
# end class

def augment_and_fill(rds:FederatedASDFDataSet, dp:DataPool, output_filename:str, dry_run=False):
    # get list of stations from DataPool to augment data from
    dp_stations = set(dp.get_stations(MIN_DATE, MAX_DATE))

    output_ds = None
    if (not dry_run):
        if (os.path.splitext(output_filename)[1].lower() != '.h5'): output_filename += '.h5'
        if (os.path.exists(output_filename)): os.remove(output_filename)
        output_ds = pyasdf.ASDFDataSet(output_filename, compression='gzip-3', mpi=False)
    # end if

    updated_netsta = set()
    for net, sta, loc, cha in dp_stations:
        # querying the reference dataset with net and sta only allows for
        # augmenting data from missing locations and channels
        ref_st, ref_et = rds.get_global_time_range(net, sta)
        dp_st, dp_et = dp.get_time_range(net, sta, loc, cha)

        if(ref_st == MAX_DATE and ref_et == MIN_DATE):
            # skip stations not in the reference dataset
            continue
        # end if

        updated_netsta.add('{}.{}'.format(net, sta))

        if(not dry_run):
            wdata_added = False

            # ========================================================
            # augment before
            # ========================================================
            if(dp_st < ref_st):
                print('Augmenting before {}.{}.{}.{}:[{} - {}]'.format(net, sta, loc, cha, dp_st, ref_st),
                      end='', flush=True)

                count = 0
                for wd in dp.waveform_iterator(net, sta, loc, cha, dp_st, ref_st):
                    output_ds.add_waveforms(wd, tag='raw_recording')
                    wdata_added = True
                    count += 1
                # end for
                print('. Added {} traces'.format(count))
            # end if

            # ========================================================
            # fill gaps
            # ========================================================
            print('Filling gaps {}.{}.{}.{}:[{} - {}]'.format(net, sta, loc, cha, ref_st, ref_et),
                  end='', flush=True)
            gaps = rds.find_gaps(network=net, station=sta,
                                 location=loc, channel=cha,
                                 start_date_ts=ref_st.timestamp,
                                 end_date_ts=ref_et.timestamp,
                                 min_gap_length=dp.min_gap_length)
            count = 0
            for gap in gaps:
                for wd in dp.waveform_iterator(gap[0], gap[1], gap[2], gap[3],
                                               UTCDateTime(gap[4]),
                                               UTCDateTime(gap[5])):
                    output_ds.add_waveforms(wd, tag='raw_recording')
                    wdata_added = True
                    count += 1
                # end for
            # end for
            print('. Added {} traces'.format(count))

            # ========================================================
            # augment after
            # ========================================================
            if(dp_et > ref_et):
                print('Augmenting after {}.{}.{}.{}:[{} - {}]'.format(net, sta, loc, cha, ref_et, dp_et),
                      end='', flush=True)

                count = 0
                for wd in dp.waveform_iterator(net, sta, loc, cha, ref_et, dp_et):
                    output_ds.add_waveforms(wd, tag='raw_recording')
                    wdata_added = True
                    count += 1
                # end for
                print('. Added {} traces'.format(count))
            # end if

            if(wdata_added):
                inventory = rds.get_inventory(network=net, station=sta)
                output_ds.add_stationxml(inventory)
            # end if

        # end if
    # end for

    if(not dry_run):
        del output_ds
    # end if

    return updated_netsta
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('asdf-ref', required=True,
                type=click.Path(exists=True))
@click.argument('data-src', required=True,
                type=click.Path(exists=True))
@click.argument('output-filename', required=True,
                type=click.Path(dir_okay=False))
@click.option('--min-gap-length', default=120, type=float, show_default=True,
              help="Minimum length of gaps in seconds to fill")
@click.option('--dry-run', default=False, is_flag=True, show_default=True,
              help="Dry run only reports stations whose data can be augmented/filled")
def process(asdf_ref, data_src, output_filename, min_gap_length=600, dry_run=False):
    """
    ASDF_REF: Text file containing paths to ASDF files, the data from which are to be used as the reference
              for finding and filling gaps\n

    DATA_SRC: Text file containing paths to ASDF files and mseed folders the data from
              which are to be used to fill gaps in the data in ASDF_REF\n

    OUTPUT_FILENAME: Output H5 file-name\n

    Example usage:
        mpirun -np 48 python fill_data_gaps.py /g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt \
        sources.txt filled.h5
    """
    assert min_gap_length > 0, '--min-gap-length must be > 0'

    rds = FederatedASDFDataSet(asdf_ref)
    dp = DataPool(data_src, min_gap_length, os.path.dirname(data_src))

    if(rds.fds.rank == 0):
        u = augment_and_fill(rds, dp, output_filename, dry_run=dry_run)
        if(dry_run):
            print('\n\nData can be augmented/filled for the following stations: {}\n\n'.format(u))
        else:
            print('\n\nData augmented/filled for the following stations: {}\n\n'.format(u))
        # end if
    # end if
# end func

if (__name__ == '__main__'):
    process()
# end if
