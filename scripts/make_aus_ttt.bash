#!/bin/env bash

# This is gmt-version6 modern syntax. Not for older gmt version 4 or 5

# Example Usage: ./make_aus_ttt.bash 02

ttt_figure="australia_ttt_horizontal_slice_$1"

gmt begin $ttt_figure png,jpeg

# gmt coast -R100/180/-50/0 -JM6i -B -W0.5p -Gchocolate 

gmtdir='/g/data/ha3/fxz547/travel_time_tomography/inversion_S1x1'
cptdir=$gmtdir'/cpt'
# set psdir  = $gmtdir'/ps2'
ddir=$gmtdir'/DATA1x1'
bounf=$gmtdir'/boundary'

fhead='Plocsol.1x1.'
spacing=2/2
# set psfile = "TTT_Australia_P_1x1_$1.ps"
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

ifile=$fhead$1

#Layers 01--23: 
#Depth KM: 5.0   22.5   52.5   90.0  135.0  185.0  235.0  285.0  335.0  385.0  435.0  485.0  535.0  585.0  635.0  685.0  760.0  860.0  960.0 1060.0 1180.0 1325.0 1500.0
depth="Horizontal Slice at Depth 22KM"
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
gmt psxy stations4202.txt -R$omin/$omax/$amin/$amax -J$XCO -B${annotation}f2nsew -St.1  

#**************************************************  Left Bottom (next figure)    **

# gmt colorbar -Ct.cpt -Dx8c/1c+w12c/0.5c+jTC+h -Bxaf+l"topography" -By+lkm


gmt colorbar -C$cptdir/palT2.cpt -D10.0/-1.2/7/0.25h -Ba 


gmt end
