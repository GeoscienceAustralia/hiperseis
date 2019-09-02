#!/bin/bash

#################################################################################################################################
#performs H-K stacking ##
#a=1.0###RRF############
# ### -P, Argument = Average crustal velocity#
###### -Z, Argument = Moho Depth search. start/end/increment#
########-K, Argument = Vp/Vs ratio search. start/end/increment#
#########-W, Argument = Weighing factor for 0p1s/2p1s/1p2s. Default = 0.7/0.2/0.1#
##########-T, Argument = Smoothing time window in sec. Default = 0.1#
############-N, Argument = n-th stacking. Default = 1#
#############-X, Argument = XC weight. Default = 0#
################# 0: No XC weight; #
##################1: All phase; 2: 0p1s+2p1s; 3: 0p1s+1p2s. As a function of Vp/Vs.#
##################4: All phase; 5: 0p1s+2p1s; 6: 0p1s+1p2s. As a function of Moho depth.#
##################7: All phase; 8: 0p1s+2p1s; 9: 0p1s+1p2s. As a function of Vp/Vs and Moho.#
##################10: Do x = 0, 1, 2, 3, 7, 8, 9 one by one.#
#############-R, Argument = BAZ range. Default = 0/360.#
############-C, Argument = Dpeak/mh/rv. The values control peak search.#
################Search for peaks p[i] whose (p_max-p[i])/p_max < Dpeak, #
################but the location of p[i] must be outside of an area defined by mh and rv. #
################Usually use default = 0.2/3.0/0.15.#
#############-I, Argument = dir of RF data.#
#############-L, Argument = RF data file list: The file name of list of time traces.#
#############-D, Argument = dirout: The dir for output file.#
#############-O, Argument = output files: sfunc and mh_rv2d.#
#############-o, No argument. To use ASCII output for sfunctions.#
#############-s, No argument. To unset Time shift.#
#############-v, No argument. To output some sfunctions to show cross correlation weight.#
##################################################################################################################################

#######X=0#########################################################


./hk2d -P 6.459 -Z 15/60/0.05 -K 1.55/1.90/0.001 -W 0.60/0.25/0.15 -T 0.1 -N 4 -X 0 \
     -I /Users/nita/GEOPHYSICS/H_K/RF_HK \
     -L /Users/nita/GEOPHYSICS/H_K/RF_HK/details_event.table \
     -D /Users/nita/GEOPHYSICS/H_K/RF_HK/Vp_tomo/BAZ_w0.360 -O sfr2d.AU.KMBL/hkr2d.AU.KMBL \

#####################X=10##########################################

./hk2d -P 6.459 -Z 15/60/0.05 -K 1.55/1.90/0.001 -W 0.60/0.25/0.15 -T 0.1 -N 4 -X
     -I /Users/nita/GEOPHYSICS/H_K/RF_HK \
     -L /Users/nita/GEOPHYSICS/H_K/RF_HK/details_event.table \
     -D /Users/nita/GEOPHYSICS/H_K/RF_HK/Vp_tomo/BAZ_w0.360 -O sfr2d.AU.KMBL/hkr2d.AU.KMBL \

./hk2d -P 6.459 -Z 15/60/0.05 -K 1.55/1.90/0.001 -W 0.60/0.25/0.15 -T 0.1 -N 4 -X 2 \
     -I /Users/nita/GEOPHYSICS/H_K/RF_HK \
     -L /Users/nita//GEOPHYSICS/H_K/RF_HK/details_event.table \
     -D /Users/nita//GEOPHYSICS/H_K/RF_HK/Vp_tomo/BAZ_w0.360 -O sfr2d.AU.KMBL/hkr2d.AU.KMBL \

./hk2d -P 6.459 -Z 15/60/0.05 -K 1.55/1.90/0.001 -W 0.60/0.25/0.15 -T 0.1 -N 4 -X 3 \
     -I /Users/nita/GEOPHYSICS/H_K/RF_HK \
     -L /Users/nita//GEOPHYSICS/H_K/RF_HK/details_event.table \
     -D /Users/nita//GEOPHYSICS/H_K/RF_HK/Vp_tomo/BAZ_w0.360 -O sfr2d.AU.KMBL/hkr2d.AU.KMBL \

./hk2d -P 6.459 -Z 15/60/0.05 -K 1.55/1.90/0.001 -W 0.60/0.25/0.15 -T 0.1 -N 4 -X 7 \
     -I /Users/nita/GEOPHYSICS/H_K/RF_HK \
     -L /Users/nita//GEOPHYSICS/H_K/RF_HK/details_event.table \
     -D /Users/nita//GEOPHYSICS/H_K/RF_HK/Vp_tomo/BAZ_w0.360 -O sfr2d.AU.KMBL/hkr2d.AU.KMBL \

./hk2d -P 6.459 -Z 15/60/0.05 -K 1.55/1.90/0.001 -W 0.60/0.25/0.15 -T 0.1 -N 4 -X 8 \
     -I /Users/nita/GEOPHYSICS/H_K/RF_HK \
     -L /Users/nita//GEOPHYSICS/H_K/RF_HK/details_event.table \
     -D /Users/nita//GEOPHYSICS/H_K/RF_HK/Vp_tomo/BAZ_w0.360 -O sfr2d.AU.KMBL/hkr2d.AU.KMBL \

./hk2d -P 6.459 -Z 15/60/0.05 -K 1.55/1.90/0.001 -W 0.60/0.25/0.15 -T 0.1 -N 4 -X 9 \
     -I /Users/nita/GEOPHYSICS/H_K/RF_HK \
     -L /Users/nita//GEOPHYSICS/H_K/RF_HK/details_event.table \
     -D /Users/nita//GEOPHYSICS/H_K/RF_HK/Vp_tomo/BAZ_w0.360 -O sfr2d.AU.KMBL/hkr2d.AU.KMBL \

./hk2d -P 6.459 -Z 15/60/0.05 -K 1.55/1.90/0.001 -W 0.60/0.25/0.15 -T 0.1 -N 4 -X 10 \
    -I /Users/nita/GEOPHYSICS/H_K/RF_HK \
    -L /Users/nita/GEOPHYSICS/H_K/RF_HK/details_event.table \
    -D /Users/nita/GEOPHYSICS/H_K/RF_HK/Vp_tomo/BAZ_w0.360 -O sfr2d.AU.KMBL/hkr2d.AU.KMBL \

##############According to RF under analysis & gaussian width can change the parameters or add similar command line#################




