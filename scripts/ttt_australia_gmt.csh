#!/bin/csh

# this script works for GMT version 4.5.11
# module load gmt/4.5.11
# ./australia2.gmt

set gmtdir = '/g/data/ha3/fxz547/travel_time_tomography/inversion_P_wave'

set cptdir = $gmtdir'/cpt'
set psdir  = $gmtdir'/ps2'
set ddir   = $gmtdir'/DATA'
set bounf  = $gmtdir'/boundary'
set fhead  = 'Plocsol.2x2.'
set spacing = 2/2
set psfile = 'TTT_Australia_P_2x2.ps'
set amin   = -50.000
set amax   =  10.00
set omin   =  95.00
set omax   = 190.0
set XCO    = m.20
set annotation = a20
set annotationy = a20
set amid   = 22.5
set omid   = 130.000
\rm $psdir/$psfile
#**************************************************     Left Top     **
set ifile  = $fhead'02'
set depth  = '35-70'
echo $ifile
#............................................... make grid file.......
xyz2grd $ddir/$ifile -G$ddir/$ifile.bin -I$spacing -R$omin/$omax/$amin/$amax -F -V -Ddegree/degree/%/1/1 -V
grdsample $ddir/$ifile.bin -G$ddir/$ifile.grd -V -F -I0.1/0.1 -R$omin/$omax/$amin/$amax
grdimage $ddir/$ifile.grd -R$omin/$omax/$amin/$amax -J$XCO -C$cptdir/palT2.cpt -X1.5  -Y5 -K -V -P >$psdir/$psfile
pscoast -R -J$XCO -O -K  -Dl -W >>$psdir/$psfile
psxy $bounf -R -J$XCO -O -K -B${annotation}f2nseW -W0.5/80/255/0 -M >>$psdir/$psfile
pstext -R0/21/0/27 -Jx1 -G0 -O -K  <<END >>$psdir/$psfile
0.1 0.2  12 0. 4 5  $depth
END

#get correct stations awk '{print $8 " " $9}' < $gmtdir/sorted_region_P.csv | sort -n | uniq > stations.txt

awk '{print $8 " " $9}' < $gmtdir/region_P.csv | sort -n | uniq > stations.txt

psxy stations.txt -R$omin/$omax/$amin/$amax -J$XCO  -O -K -B${annotation}f2nsew -St.1  >>$psdir/$psfile
