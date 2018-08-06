#!/usr/bin/env bash
# Fei Zhang
# to convert the csh gmt script into bash and using newer version of GMT (V5 and V6)

GMT=gmt
gmtdir='/g/data/ha3/fxz547/travel_time_tomography/inversion_P_wave'
cptdir=$gmtdir'/cpt'
psdir=$gmtdir'/ps2'
ddir=$gmtdir'/DATA'
bounf=$gmtdir'/boundary'
fhead='Plocsol.2x2.'
spacing=2/2
psfile='TTT_Australia_P_2x2.ps'
amin=-50.000
amax=10.00
omin=95.00
oma=190.0
XCO=m.06
annotation=a20
annotationy=a20
amid=22.5
omid=130.000

rm $psdir/$psfile

#**************************************************     Left Top     **
ifile=$fhead'02'
depth='35-70'
echo $ifile
#............................................... make grid file.......
xyz2grd $ddir/$ifile -G$ddir/$ifile.bin -I$spacing -R$omin/$omax/$amin/$amax -F -V -Ddegree/degree/%/1/1 -V

grdsample $ddir/$ifile.bin -G$ddir/$ifile.grd -V -F -I0.1/0.1 -R$omin/$omax/$amin/$amax

grdimage $ddir/$ifile.grd -R$omin/$omax/$amin/$amax -J$XCO -C$cptdir/palT2.cpt -X1.5  -Y22 -K -V -P >$psdir/$psfile

pscoast -R -J$XCO -O -K  -Dl -W >>$psdir/$psfile

psxy $bounf -R -J$XCO -O -K -B${annotation}f2nseW -W0.5/80/255/0 -M >>$psdir/$psfile

pstext -R0/21/0/27 -Jx1 -G0 -O -K  <<END >>$psdir/$psfile

0.1 0.2  12 0. 4 5  $depth

END
#awk '{print $8 " " $9}' < ../sorted_region_P.csv | sort -n | uniq > stations.txt
awk '{print $8 " " $9}' < $gmtdir/region_P.csv | sort -n | uniq > stations.txt

psxy stations.txt -R$omin/$omax/$amin/$amax -J$XCO  -O -K -B${annotation}f2nsew -St.1  >>$psdir/$psfile

#**************************************************     Left Bottom     **
#ifile  = $fhead'03'
#depth  = '70-110'
#echo $ifile