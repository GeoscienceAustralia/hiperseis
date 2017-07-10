import pyasdf
import time
import os
from xcorqc.xcorqc import IntervalStackXCorr, saveXCorrPlot
from obspy import Stream, Trace
import matplotlib.pyplot as plt
from obspy import UTCDateTime


code_start_time = time.time()


# =========================== User Input Required =========================== #

#Path to the ASDF file for xcor
data_path = '/g/data/ha3/g.h5'
out = '/g/data/ha3/out.h5'

temp_net = "XX"

# =========================================================================== #


if os.path.exists(out):
    os.remove(out)

ds = pyasdf.ASDFDataSet(data_path)
print(ds.waveforms.list())

# dictionary with unique ids as keys and stream objects (for permanent station data) as values
id_st_dict = {}

# first extract obspy stream data from the ASDF file that is from a permanent station (i.e. not the temporary network)
for net_sta in ds.waveforms.list():
    iter_net = net_sta.split('.')[0]

    if not temp_net == iter_net:
        # get the station accessor
        iter_sta_accessor = ds.waveforms[net_sta]

        # get a list of the waveforms
        iter_waveforms = iter_sta_accessor.list()

        for wave in iter_waveforms:
            if not wave == "StationXML":
                iter_st = iter_sta_accessor[wave]
                print(iter_st)

                # add it to the dict
                tag = wave.split('__')[3]

                id_st_dict[tag] = iter_st


print(id_st_dict)

print('')

def test_process(st, inv):

    xcor_st = Stream()

    tr1 = Stream([st[0]])
    tr2 = Stream([st[1]])

    y, x, comp = IntervalStackXCorr(tr1, tr2)

    for day_stack_xcor in y:
        print(type(day_stack_xcor))
        print(day_stack_xcor[0])

        # fill in headers
        stats = {'network': tr1[0].stats.network, 'station': tr1[0].stats.station, 'location': '',
                 'channel': tr1[0].stats.channel, 'npts': len(day_stack_xcor[0]), 'sampling_rate': 100,
                 'mseed': {'dataquality': 'D'}}

        stats["starttime"] = tr1[0].stats.starttime

        xcor_st += Trace(data=day_stack_xcor[0], header=stats)


    #
    # for tr in st:
    #     temp_st = Stream(traces=[tr])
    #     print('')
    #
    #     # get the uid label
    #     uid_label = tr.stats.asdf.labels[1]
    #     perm_sta = tr.stats.asdf.labels[2]
    #
    #     print(uid_label, perm_sta)
    #
    #
    #     perm_st = id_st_dict[uid_label]
    #
    #     print(temp_st)
    #     print(perm_st)
    #     print("DO XCOR......")
    #
    #
    #     y, x, comp = IntervalStackXCorr(perm_st, temp_st)
    #
    #
    #
    #     for day_stack_xcor in y:
    #         print(type(day_stack_xcor))
    #         print(day_stack_xcor[0])
    #
    #         # fill in headers
    #         stats = {'network': tr.stats.network, 'station': tr.stats.station, 'location': '',
    #                  'channel': tr.stats.channel, 'npts': len(day_stack_xcor[0]), 'sampling_rate': 100,
    #                  'mseed': {'dataquality': 'D'}}
    #
    #         stats["starttime"] = tr.stats.starttime
    #
    #         xcor_st += Trace(data=day_stack_xcor[0], header=stats)


    return xcor_st



ds.process(process_function=test_process,
           output_filename=out,
           tag_map={"raw_recording": "test_write"})


del ds