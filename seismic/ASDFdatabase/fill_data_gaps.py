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
import os
from mpi4py import MPI
from obspy.core import Stream

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

        # instantiate sources
        self.asdf_source_fn = os.path.join(self.scratch_path, 'asdf_files.txt')

        lines = set(open(self.sources_txt_file).readlines())
        h5_lines = set([line for line in lines if ('h5' or 'H5') in line])
        mseed_lines = lines.difference(h5_lines)

        h5_lines = list(h5_lines)
        mseed_lines = list(mseed_lines)

        if (self.rank == 0):
            open(self.asdf_source_fn, 'w').writelines(h5_lines)
        # end if
        self.comm.barrier()

        self.sources.append(FederatedASDFDataSet(self.asdf_source_fn))
        for mseed_folder in mseed_lines:
            self.sources.append(MseedIndex(mseed_folder.strip(), self.mseed_pattern))
        # end for
    # end func

    def get_waveforms(self, net, sta, loc, cha, st:UTCDateTime, et:UTCDateTime):
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
                if(after>before):
                    print('Found traces in 1')
                    print(wd)

                gaps = get_gaps(wd)
                print("----1", gaps)
            elif(len(gaps)):
                before = len(wd)
                for gap in gaps:
                    wd += src.get_waveforms(net, sta, loc, cha, gap[4], gap[5])
                # end for
                after = len(wd)
                if(after>before):
                    print('Found traces in 2')
                    print(wd)

                gaps = get_gaps(wd)
                print("----2", gaps)
            # end if
        # end for
        return wd
    # end func
# end class

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('asdf-ref', required=True,
                type=click.Path(exists=True))
@click.argument('data-src', required=True,
                type=click.Path(exists=True))
@click.argument('output-filename', required=True,
                type=click.Path(dir_okay=False))
@click.option('--min-gap-length', default=600, type=float, show_default=True,
              help="Minimum length of gaps in seconds to fill")
def process(asdf_ref, data_src, output_filename, min_gap_length):
    """
    ASDF_REF: Text file containing paths to ASDF files, the data from which are to be used as the reference
              for finding and filling gaps\n

    DATA_SRC: Text file containing paths to ASDF files and mseed folders the data from
              which are to be used to fill gaps in the data in ASDF_REF\n

    OUTPUT_FILENAME: Output H5 file-name\n

    Example usage:

    """
    assert min_gap_length > 0, '--min-gap-length must be > 0'

    rds = FederatedASDFDataSet(asdf_ref)
    dp = DataPool(data_src, min_gap_length, os.path.dirname(data_src))

    if(rds.fds.rank == 0):
        gaps = rds.find_gaps(network='AU', station='AXCOZ', min_gap_length=min_gap_length)
        for gap in gaps:
            print("******", gap[0], gap[1], gap[2], gap[3], UTCDateTime(gap[4]), UTCDateTime(gap[5]))
            wd = dp.get_waveforms(gap[0], gap[1], gap[2], gap[3],
                                  UTCDateTime(gap[4]),
                                  UTCDateTime(gap[5]))
# end func

if (__name__ == '__main__'):
    process()
# end if
