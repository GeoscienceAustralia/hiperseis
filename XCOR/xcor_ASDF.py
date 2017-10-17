import pyasdf
import time
import os
from xcorqc.xcorqc import xcorr2
from obspy import Stream, Trace
import matplotlib.pyplot as plt
from obspy import UTCDateTime


code_start_time = time.time()


# =========================== User Input Required =========================== #

#Path to the ASDF file for xcor
data_path = '/g/data1/ha3/Passive/_AusArray/OA/processing/test_xcor_ashby/two_stns_xcor_test.h5'
out = '/g/data1/ha3/Passive/_AusArray/OA/processing/test_xcor_ashby/two_stns_xcor_out.h5'

temp_net = "OA"

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
                # print(iter_st)

                # add it to the dict
                tag = wave.split('__')[3]

                id_st_dict[tag] = iter_st


# print(id_st_dict)

print('')

def xcor_process(st, inv):

    xcor_st = Stream()

    # tr1 = Stream([st[0]])
    # tr2 = Stream([st[1]])
    #
    # y, x, comp = IntervalStackXCorr(tr1, tr2)
    #
    # for day_stack_xcor in y:
    #     print(type(day_stack_xcor))
    #     print(day_stack_xcor[0])
    #
    #     # fill in headers
    #     stats = {'network': tr1[0].stats.network, 'station': tr1[0].stats.station, 'location': '',
    #              'channel': tr1[0].stats.channel, 'npts': len(day_stack_xcor[0]), 'sampling_rate': 100,
    #              'mseed': {'dataquality': 'D'}}
    #
    #     stats["starttime"] = tr1[0].stats.starttime
    #
    #     xcor_st += Trace(data=day_stack_xcor[0], header=stats)



    for tr in st:
        temp_st = Stream(traces=[tr])
        print('')
        print(tr.stats.asdf.labels)

        # get the uid label
        uid_label = tr.stats.asdf.labels[1]

        print(uid_label)


        perm_st = id_st_dict[uid_label]

        temp_tr = temp_st[0]
        ref_tr = perm_st[0]

        stationPair = ref_tr.stats.station + '.' + temp_tr.stats.station

        print(temp_st)
        print(perm_st)
        print("DO XCOR......")

        xcl, winsPerInterval = xcorr2(ref_tr, temp_tr)

        if (xcl is None):
            print("\t\tWarning: no cross-correlation results returned for station-pair %s, " %
                  (stationPair) + " due to gaps in data.")
            continue


        print(xcl)
        print(winsPerInterval)

        # saveXCorrPlot(y, x, '/g/data/ha3/', 'test_plot', comp)


        # fill in headers for new xcor trace
        stats = {'network': tr.stats.network, 'station': tr.stats.station, 'location': "",
                 'channel': uid_label[2:], 'npts': len(xcl[0]), 'sampling_rate': tr.stats.sampling_rate,
                 'mseed': {'dataquality': 'D'}, "asdf": {}}


        stats["starttime"] = tr.stats.starttime

        temp_tr = Trace(data=xcl[0], header=stats)

        temp_tr.stats.asdf.labels = tr.stats.asdf.labels

        xcor_st += temp_tr


        # break


    return xcor_st



ds.process(process_function=xcor_process,
           output_filename=out,
           tag_map={"raw_recording": "xcor"})


# load in new ASDF that has been written
new_ds = pyasdf.ASDFDataSet(out)

# now just copy over waveforms from input ASDF file to new ASDF
for sta_id in ds.waveforms.list():
    sta_acc = ds.waveforms[sta_id]
    for asdf_st_id in sta_acc.list():
        if not asdf_st_id == "StationXML":
            new_ds.add_waveforms(sta_acc[asdf_st_id],tag=asdf_st_id.split('__')[3])
        else:
            new_ds.add_stationxml(sta_acc[asdf_st_id])




del new_ds

del ds