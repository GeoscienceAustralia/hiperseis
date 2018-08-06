#!/bin/env bash

# This is gmt-6 modern syntax.
# NOT for older gmt version 4 or 5
# can run in VDI if gmt6 installed

# Example Usage: ./make_aus_ttt.bash 02

# User, please set P or S wave below
P_S='S'

ttt_figure="${P_S}_ttt_horizontal_slice_$1"

# gmt begin $ttt_figure png,jpeg
gmt begin $ttt_figure png

# User set variables
gmtdir="/g/data/ha3/fxz547/travel_time_tomography/inversion_${P_S}1x1"
ddir=$gmtdir'/DATA'

cptdir=$gmtdir'/cpt'
bounf=$gmtdir'/boundary'

fhead='Plocsol.1x1.'

# spacing will affect images look
#spacing=1/1
spacing=2/2  # original 2/2 degree

amin=-50.000
amax=10.00
omin=95.00
omax=190.0
XCO=m.20
annotation=a20
annotationy=a20
amid=22.5
omid=130.000

#*********************************************
# gmt coast -R100/180/-50/0 -JM6i -B -W0.5p -Gchocolate 

ifile=$fhead$1

# Model Layers 01--23: 

KM_DEPTH=('5.0' '22.5'  '52.5'  '90.0' '135.0' '185.0' '235.0'  '285.0' '335.0' '385.0' '435.0' '485.0' '535.0'  '585.0' '635.0'  '685.0'  '760.0' 
 '860.0'  '960.0'  '1060.0'  '1180.0'  '1325.0'  '1500.0');

index=$((10#$1 - 1))  #base 10 the number like 02, offset 1
indepth=${KM_DEPTH[$index]}

echo The index is $index at $indepth

# user change the title lable S or P wvae?
depth="${P_S} - wave velocity profile  depth=$indepth KM"

echo plotting the $ifile

#............................................... make grid file.......
gmt xyz2grd $ddir/$ifile -G$ddir/$ifile.bin -I$spacing -R$omin/$omax/$amin/$amax -r -V -Ddegree/degree/%/1/1
gmt grdsample $ddir/$ifile.bin -G$ddir/$ifile.grd -V -r -I0.1/0.1 -R$omin/$omax/$amin/$amax
gmt grdimage $ddir/$ifile.grd -R$omin/$omax/$amin/$amax -J$XCO -C$cptdir/palT2.cpt -X1.5  -Y5 -V  
gmt coast  -J$XCO  -Dl -W0.5p #-Gchocolate

gmt plot $bounf -J$XCO  -B${annotation}f2nseW -W0.5/80/255/0

gmt text -R0/21/0/27 -Jx1   << EOF

0.1 0.2  12 0. 4 5  $depth

EOF


#awk '{print $8 " " $9}' < $gmtdir/sorted_region_P.csv | sort -n | uniq > stations.txt
#awk '{print $8 " " $9}' < $gmtdir/region_P.csv | sort -n | uniq > stations.txt
gmt psxy stations.txt -R$omin/$omax/$amin/$amax -J$XCO -B${annotation}f2nsew -St.1  

#**************************************************  Left Bottom (next figure)    **

# gmt colorbar -Ct.cpt -Dx8c/1c+w12c/0.5c+jTC+h -Bxaf+l"topography" -By+lkm

#gmt colorbar -C$cptdir/palT2.cpt -D10.0/-1.2/7/0.25h -Ba  # original
gmt colorbar -C$cptdir/pal3-T2.cpt -D10.0/-1.2/7/0.25h -Ba  # [-3 +3]


gmt end
