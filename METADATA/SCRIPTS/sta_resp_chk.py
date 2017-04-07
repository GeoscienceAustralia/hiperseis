#!/usr/bin/env python
import os
import argparse
from obspy import UTCDateTime, read_inventory, read
import numpy as np
from scipy import signal as sg
import matplotlib.pyplot as plt
from obspy.xseed import Parser
from scipy import signal
from obspy.neic.client import Client
from obspy.xseed.utils import SEEDParserException
from PIL import Image, ImageDraw, ImageMath, ImageOps
from matplotlib.backends.backend_agg import FigureCanvasAgg
import sys
from matplotlib import gridspec

# This is a script to analyse seismic noise in order to check correctness of the seismic station response file aka "metadata"
# it reads metadata XML file and seismic waveforms as a file or from NEIC server (CWB)
# file location is taken from ROOT_DIR or taken from the command line via args
# There are two hard-wired files HNM.model and LNM.model used by this code in order to compare results agains McNamara standard noise models.

# Standard High and Low noise models from McNamara
std_noise = os.path.join(os.path.dirname(__file__), "./DATA/std_global_noise")


parser = argparse.ArgumentParser(
    description="Calculation of spectral density plots for one station. It outputs by default a jpg image with station name and its state as a RED or GREEN prefix. ",
    epilog='All parameters are defined at the beggining of the code. Please look inside to be sure that they satisfy your requirements.'
)

parser.add_argument(
    '--s',
    dest='socket_list',
    default='http://rhe-eqm-event-dev.dev.lan:8080',
    type=str,
    help='HTTP address of CWB.')
parser.add_argument(
    '--r',
    dest='resp_file',
    default='',
    type=str,
    help='File name of response file in XML format.')
parser.add_argument(
    '--f',
    dest='file_name',
    default='\n',
    type=str,
    help='File name of waveform data.')
parser.add_argument(
    '--start',
    dest='t_start',
    default='2016-5-14T0:0:0',
    help='Start time as 2016-5-14T0:0:0')
parser.add_argument(
    '--end',
    dest='t_end',
    default='3000-12-12T0:0:0',
    help='End time as 2016-5-14T0:0:0')
parser.add_argument(
    '--tw', dest='t_w', default='0.25', help='Running window in hours.')
parser.add_argument('--net', dest='net', default='AU', help='Network name.')
parser.add_argument('--sta', dest='sta', default='CNB..', help='Station name.')
parser.add_argument('--loc', dest='loc', default='..', help='Location code.')
parser.add_argument(
    '--ch',
    dest='chan',
    default='BHZ',
    help='Location code.\n Remember that network,station,loc, and chan must compose all together 12 character string in REGEXP format.\n There is no .* allowed.'
)
parser.add_argument(
    '--plot', dest='plot', default=False, help='Plot window, default is False')

args = parser.parse_args()

file_name = args.file_name
net = args.net
sta = args.sta
loc = args.loc
chan = args.chan
plot = args.plot
t_start = UTCDateTime(args.t_start)
t_end = UTCDateTime(args.t_end)
resp_file = args.resp_file
if not resp_file:
    dbmaster="./DATA/"+net+".xml"
# XML response file
#print ROOT_DIR
    resp_file = os.path.join(os.path.dirname(__file__), dbmaster)

if len(resp_file) > 3:
    if os.path.isfile(resp_file):
        inv = read_inventory(resp_file)
    else:
        print "File ", resp_file, " does not exist."
        exit(-1)
else:

    exit(-1)

if len(file_name) > 2:
    station = read(file_name)
    t_start = station[0].stats.starttime
else:

    # checking times and number of time windows to process
    if UTCDateTime(t_end) > UTCDateTime.utcnow():
        t_end = UTCDateTime.utcnow() - 60.
        t_start = t_end - 1. * 60. * 60.

    parser.print_help()

    # initialize the cwb port
    #client=Client(host='10.7.161.60',port=2061,debug=False,nonice=True)
    client = Client(host='10.7.161.60', port=2061, debug=False)
    #client=Client(host='136.177.24.70',port=2061,debug=False)

    # read waveform
    #
    print net, sta, loc, chan, t_start, '-d ', int(t_end - t_start)
    station = client.getWaveform(net, sta, loc, chan, t_start, t_end)

print station
if len(station) == 0:
    exit(-1)
# length of running window in seconds
# McNamara uses 15 minutes windows
time_win = np.float(args.t_w) * 60. * 60.

# length of running window in samples
t_w = np.ceil(time_win * station[0].stats.sampling_rate)
sampling = station[0].stats.sampling_rate
# Total length of the seismic record
time_diff = station[0].stats.endtime - station[0].stats.starttime

# number of running windows to process
# McNamara uses 75% overlap
time_shift=time_win*0.25
total_win = np.int(np.ceil(time_diff / (time_shift)))
#print '---', total_win,time_diff,t_w
station.merge(method=0,fill_value='interpolate',interpolation_samples=0)
if isinstance(station[0].data,np.ma.masked_array):
    station[0].data=station[0].data.filled()

station.detrend()
seedid = station[0].getId()

# remove instrument
err = station[0].attach_response(
    inv.select(
        network=station[0].stats.network,
        station=station[0].stats.station,
        channel=station[0].stats.channel))
if err:
    print "Can not find responce for ", seedid
    exit - 1
err = station[0].remove_response(output='ACC')



# read ambient noise models
lnm = []
x_lnm = []
hnm = []
x_hnm = []

with open(std_noise+"/LNM.model") as f:
    for line in f:
        x, y = line.split()
        x_lnm.append(float(x))
        lnm.append(float(y))

with open(std_noise+"/HNM.model") as f:
    for line in f:
        x, y = line.split()
        x_hnm.append(float(x))
        hnm.append(float(y))

frq_max = 0
amp_max = 0
fig1 = plt.figure(figsize=(10, 10), dpi=200, facecolor='white')
i_fig = 0
stack = []
hnm_cntrl = np.array(hnm)[::-1]
hnm_per = np.array(x_hnm)[::-1]
lnm_cntrl = np.array(lnm)[::-1]
lnm_per = np.array(x_lnm)[::-1]
hr = np.int(time_diff / 3600.)
mn = np.int((time_diff - hr * 3600) / 60.)
sec = time_diff - hr * 3600 - mn * 60
print "Loaded ", hr, " hr ", mn, " min ", sec, " sec of waveform data."
print "It will analyze ", total_win - 1, " windows of ", time_win, " sec length"
ax = plt.gca()
for i in xrange(total_win - 1):
    sys.stdout.write("Processing %d window\r" % (i))
    sys.stdout.flush()
    amp = []
    data_slice = 0
    data_slice = station[0].slice(t_start + i * time_shift, t_start + (i*time_shift + time_win)).copy()
    data_slice=sg.detrend(data_slice.data)
    num = np.size(data_slice)
    n_total=np.power(2,np.int(np.log2(num))+1)
    if num > 0:
        taper = np.hanning(np.size(data_slice))
        data_slice=(data_slice-np.mean(data_slice))*taper
        amp = np.fft.rfft(data_slice - np.mean(data_slice) * taper,n_total)
        num=n_total
        amp = amp[0:int(num / 2) + 1] * np.conjugate(amp[0:int(num / 2) + 1])
        norm = 2.0 * station[0].stats.delta / float(n_total)
        amp = np.abs(amp[0:int(num / 2) + 1]) * norm

        frq=np.fft.rfftfreq(n_total,station[0].stats.delta)
        frq=frq[0:int(num/2)+1]
        frq[1:] = 1. / frq[1:]
        frq[0]=time_win
        p_ma = num / sampling
        p_mi = 1 / (sampling / 2.)

        amp = signal.medfilt(amp, kernel_size=99)
        ''' convert to dB '''
        amp = 10. * np.log10(amp)
        plt.clf()

        frq = frq[::-1]
        amp = amp[::-1]
        amp_int = np.interp(hnm_per, frq, amp.flatten())
        stack.append(amp_int)

        bandwidth = (frq >= p_mi) & (frq <= p_ma)
        frq = frq[bandwidth]
        amp = amp[bandwidth]

        spectr = plt.plot(frq, amp, 'bo')
        l_mod = plt.plot(x_lnm, lnm, 'r')
        h_mod = plt.plot(x_hnm, hnm, 'r')
        plt.setp(l_mod, color='r', linewidth='4.')
        plt.setp(h_mod, color='r', linewidth='4.')
        plt.xlim(0.1, 100)
        plt.ylim(-200, -20)
        xl = plt.xlabel("Period (sec)", fontsize=30)
        yl = plt.ylabel("Power Db", fontsize=30)
        plt.xscale('log')
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(40)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(40)
        plt.tight_layout()
        canvas = plt.get_current_fig_manager().canvas

        agg = canvas.switch_backends(FigureCanvasAgg)
        agg.draw()
        s = agg.tostring_rgb()
        # get the width and the height to resize the matrix if required
        l, b, w, h = agg.figure.bbox.bounds
        w, h = int(w), int(h)

        X = np.fromstring(s, np.uint8)
        X.shape = h, w, 3

        try:
            img = Image.fromstring("RGB", (w, h), s)
        except Exception:
            img = Image.frombytes("RGB", (w, h), s)
        img = ImageOps.grayscale(img)
        if i > 0:
            out = ImageMath.eval("(a+b)", a=out, b=img)

        else:
            out = img

        i_fig = i_fig + 1

#median_amp = np.median(np.array(stack), axis=0)
median_amp = np.average(np.array(stack), axis=0)
# lets use period band from 1sec to 10 sec for comparison coz it contains higher energy and consequently more stable results
band = (hnm_per > 1) & (hnm_per < 10)
res1 = np.sum(hnm_cntrl[band] - median_amp[band])
res2 = np.sum(lnm_cntrl[band] - median_amp[band])
flag = ''
if res1 < 0. or res2 > 0.:
    flag = 'RED-'
else:
    flag = 'GREEN-'

out = ImageMath.eval("a/i_fig", a=out, i_fig=i_fig)
plt.plot(hnm_per, median_amp, color='g', linewidth='4.',)
plt.setp(spectr, color='w', marker='o', mec='w')
plt.setp(l_mod, color='r', linewidth='4.')
plt.setp(h_mod, color='r', linewidth='4.')
plt.xlim(0.1, 100)
plt.ylim(-200, -20)
plt.xlabel("Period (sec)", fontsize=30)
plt.ylabel("Power Db", fontsize=30)
plt.xscale('log')
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(40)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(40)
plt.tight_layout()
canvas = plt.get_current_fig_manager().canvas
agg = canvas.switch_backends(FigureCanvasAgg)
agg.draw()
s = agg.tostring_rgb()
l, b, w, h = agg.figure.bbox.bounds
w, h = int(w), int(h)

X = np.fromstring(s, np.uint8)
X.shape = h, w, 3

try:
    img = Image.fromstring("RGB", (w, h), s)
except Exception:
    img = Image.frombytes("RGB", (w, h), s)

img = img.convert('RGBA')
out_w, out_h = out.size
img_w, img_h = img.size

offset = ((out_w - img_w) / 2, (out_h - img_h) / 2)

unpack_out = np.array(out)
unpack_img = np.array(img)

red, green, blue, alpha = unpack_img.T
white_areas = (red == 255) & (blue == 255) & (green == 255)
red_lines = (red == 255) & (blue == 0) & (green == 0)
red = np.zeros_like(unpack_out[white_areas.T])
green = np.zeros_like(unpack_out[white_areas.T])

dark_blue = np.min(unpack_out[white_areas.T])
saturated_blue = (unpack_img[..., 2] > dark_blue) & (
    unpack_img[..., 2] <= dark_blue * 1.1)

unpack_img[..., 0][white_areas.T] = unpack_out[white_areas.T] - dark_blue
unpack_img[..., 1][white_areas.T] = unpack_out[white_areas.T] - dark_blue
unpack_img[..., 2][white_areas.T] = unpack_out[white_areas.T]

gray_mask = (unpack_out.T == 255) & white_areas
unpack_img[..., 0][gray_mask.T] = 255
unpack_img[..., 1][gray_mask.T] = 255
unpack_img[..., 2][gray_mask.T] = 255
unpack_img[..., 0][red_lines.T] = 255
unpack_img[..., 1][red_lines.T] = 0
unpack_img[..., 2][red_lines.T] = 0

# this part supposed to find highest saturation and paint it yellow
unpack_img[..., 0][saturated_blue] = 255
unpack_img[..., 1][saturated_blue] = 240
unpack_img[..., 2][saturated_blue] = 40

out = Image.fromarray(unpack_img, mode='RGBA')

plt.close(fig1)
fig2 = plt.figure(facecolor='white')
inv_stat = inv.select(
    network=station[0].stats.network,
    station=station[0].stats.station,
    channel=station[0].stats.channel)
resp = inv_stat[0][0][0].response
''' although we can use subplot axes in response file plotting it seems to be hard to scale both images properly, therefore we create just resp image '''

resp.plot(0.001, output="VEL", show=False)

canvas = plt.get_current_fig_manager().canvas
agg = canvas.switch_backends(FigureCanvasAgg)
agg.draw()
s = agg.tostring_rgb()
# get the width and the height to resize the matrix
l, b, w, h = agg.figure.bbox.bounds
w, h = int(w), int(h)
X = np.fromstring(s, np.uint8)
X.shape = h, w, 3
try:
    img = Image.fromstring("RGB", (w, h), s)
except Exception:
    img = Image.frombytes("RGB", (w, h), s)

plt.close(fig2)
plt.close('all')
''' plot two images together '''
fig_l = plt.figure(dpi=100, facecolor='white')
gs = gridspec.GridSpec(2, 1)
ax1 = fig_l.add_subplot(gs[0])
plt.imshow(out)
ax1.axis('off')
ax2 = fig_l.add_subplot(gs[1])
plt.imshow(img)
ax2.axis('off')
plt.tight_layout()
plt.savefig(flag + seedid + ".jpg", facecolor="white")
if plot:
    plt.show(fig_l)
