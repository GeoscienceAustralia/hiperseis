#!/usr/bin/perl
#
#
# Plot GT residuals with ak135 and RSTT predictions
# Istvan Bondar, 2014
#
#
$|=1;
use strict;
use PDL;
use PDL::NiceSlice;
#use Copula;
use GMTfuncs;
#use MyDistributions;

my $pi = 4 * atan2(1,1);
my $DEG2RAD = $pi / 180;
my $RAD2DEG = 180 / $pi;
my $RADIUS_EARTH = 6371.0;
my $DEG2KM = $DEG2RAD * $RADIUS_EARTH;

if (@ARGV < 1) {
   print "  Usage: AlexeiMockupFigs.pl queryresult\n";
   print "Example: ./AlexeiMockupFigs.pl mockup.txt\n";
   exit -1;
}

my ($infil) = @ARGV;
my $n =0;
my @w = ();
#
# read infil
#
open(IN, $infil) || die "cannot open $infil!\n";
while (<IN>) {
    push @w,  [ split ];  # splits on spaces and stores result in a 2D array
}
close IN;
$n = @w;

gmtdefaults();
`gmtset PAGE_ORIENTATION landscape X_ORIGIN 1.5 Y_ORIGIN 1.5`;

#
# plot locations
#
PlotMap($n, \@w);



# eventid,agency,OT,microsec,OT_uncertainty,lat,lon,sminax,smajax,strike,depth,depth_uncertainty,depthType,nass,ndef,nsta,rms,gap,sgap

sub PlotMap {
    my ($n, $w) = @_;
    my ($i, $evid, $galat, $galon, $aklat, $aklon, $rslat, $rslon);
    my ($r, $kmbar, $xtic, $ytic, $psfil, $bx, $smin, $smaj, $strk, $x);
#print "$n\n$w[1][6],$w[1][7]\n";
    for ($i = 0; $i < $n; $i += 3) {
        $evid = $w[$i][0];
        ($galat, $galon) = ($w[$i][6], $w[$i][7]);
        ($aklat, $aklon) = ($w[$i+1][6], $w[$i+1][7]);
        ($rslat, $rslon) = ($w[$i+2][6], $w[$i+2][7]);
#
#       projection box
#
       (undef, undef, $psfil) = split '\/', $evid;
       $evid =~ s/\//\-/g;
       $evid =~ s/\:/\-/g;
       $psfil .= ".ps";
print "$i $evid\n";
       ($r, $kmbar, $xtic, $ytic, undef) = projbox(1, 4.5, $galat, $galon);
       $bx = "-Ba$xtic"."g$xtic/a$ytic"."g$ytic:.\"$evid\":WeSn";
       `pscoast $r $bx $kmbar -Di -N1 -Wthin -K > $psfil`;
       `psxy $r -Sa0.1 -Ggreen -Wthinnest -O -K << END >> $psfil\n$galon $galat\nEND\n`;
       `psxy $r -Sa0.1 -Gblue -Wthinnest -O -K << END >> $psfil\n$aklon $aklat\nEND\n`;
       `psxy $r -Sa0.1 -Gred -Wthinnest -O -K << END >> $psfil\n$rslon $rslat\nEND\n`;
       ($smin, $smaj, $strk) = ($w[$i][8], $w[$i][9], $w[$i][10]);
       $x = sprintf "%f %f %f", $strk, 2 * $smin, 2 * $smaj;
       `psxy $r -SE -Wthin,green -O -K << END >> $psfil\n$galon $galat $x\nEND\n`;
       ($smin, $smaj, $strk) = ($w[$i+1][8], $w[$i+1][9], $w[$i+1][10]);
       $x = sprintf "%f %f %f", $strk, 2 * $smin, 2 * $smaj;
       `psxy $r -SE -Wthin,blue -O -K << END >> $psfil\n$aklon $aklat $x\nEND\n`;
       ($smin, $smaj, $strk) = ($w[$i+2][8], $w[$i+2][9], $w[$i+2][10]);
       $x = sprintf "%f %f %f", $strk, 2 * $smin, 2 * $smaj;
       `psxy $r -SE -Wthin,red -O -K << END >> $psfil\n$rslon $rslat $x\nEND\n`;
       `psxy /dev/null $r -O >> $psfil`;
    }
    $psfil = "all.ps";
    $r = "-R113/172/-47/-4 -JM7";
    $bx = "-Ba10g10WeSn";
    `pscoast $r $bx -Dc -Wthin -K > $psfil`;
    for ($i = 0; $i < $n; $i+=3) {
        ($galat, $galon) = ($w[$i][6], $w[$i][7]);
        ($aklat, $aklon) = ($w[$i+1][6], $w[$i+1][7]);
        ($rslat, $rslon) = ($w[$i+2][6], $w[$i+2][7]);
       `psxy $r -Sc0.07 -Ggreen -Wthinnest -O -K << END >> $psfil\n$galon $galat\nEND\n`;
       `psxy $r -Sc0.07 -Gblue -Wthinnest -O -K << END >> $psfil\n$aklon $aklat\nEND\n`;
       `psxy $r -Sc0.07 -Gred -Wthinnest -O -K << END >> $psfil\n$rslon $rslat\nEND\n`;
    }
    `psxy /dev/null $r -O >> $psfil`;
}



#
# Pg/Pb/Pn, Sg/Sb/Sn GT residuals mapview
#     picks: sta phase delta esaz timeres timedef hypid
#
sub map_gtres {
    my ($name, $clat, $clon, $n, $m, $inpicks, $outpicks, $net, $gthypo, $cpt) = @_;
    my ($lat, $lon, $sta, $i, $k, $tres, $hypid, $phase);
    my ($psfil, $pipe, $r, $bx, $kmbar, $xtic, $ytic);
    my (%pg, %pb, %pn, $x);
    my $mpar = "-I1/thinner,LightSkyBlue -CLightSkyBlue -W1/thin -Givory -Sazure";
#
#   projection box
#
    ($r, $kmbar, $xtic, $ytic, undef) = projbox(15, 4.5, $clat, $clon);
#
#   Pb/Pg/Pn GT residuals, ak135
#
    $psfil = "$name.PgPnResMap.ps";
    %pb = %pg = %pn = ();
    $bx = "-Ba$xtic"."g$xtic/a$ytic"."g$ytic:.\"$name, ak135\":WeSn";
    `pscoast $r $bx $kmbar $mpar -Di -N1 -K > $psfil`;
    $pipe = "psxy $r -Sa0.1 -Ggreen -Wthinnest -O -K >> $psfil";
    open (PIPE,"| $pipe") || die "ERROR : Cannot open pipe!\n";
    for $hypid (keys %$gthypo) {
        next if ($gthypo->{$hypid}[5] ne $name);
        printf PIPE "%f %f\n", $gthypo->{$hypid}[3], $gthypo->{$hypid}[2];
        for ($i = 0; $i < $n; $i++) {
            next unless ($inpicks->[$i][6] == $hypid);
            ($sta, $tres) = ($inpicks->[$i][0], $inpicks->[$i][4]);
            $phase = $inpicks->[$i][1];
            push @{$pb{$sta}}, $tres if ($phase eq "Pb");
            push @{$pg{$sta}}, $tres if ($phase eq "Pg");
            push @{$pn{$sta}}, $tres if ($phase eq "Pn");
        }
    }
    close PIPE;
#
#   Plot Pb GT residuals as circles
#
    $pipe = "psxy $r -Sc0.12 -C$cpt -Wthinnest -L -O -K >> $psfil";
    open (PIPE,"| $pipe") || die "ERROR : Cannot open pipe!\n";
    for $sta (keys %pb) {
        ($lon, $lat, undef) = split ' ', $net->{$sta};
        $k = @{$pb{$sta}};
        $x = pdl(@{$pb{$sta}});
        printf PIPE "%f %f %f\n", $lon, $lat, median($x);
    }
    close PIPE;
#
#   Plot Pg GT residuals as inverted triangles
#
    $pipe = "psxy $r -Si0.13 -C$cpt -Wthinnest -L -O -K >> $psfil";
    open (PIPE,"| $pipe") || die "ERROR : Cannot open pipe!\n";
    for $sta (keys %pg) {
        ($lon, $lat, undef) = split ' ', $net->{$sta};
        $k = @{$pg{$sta}};
        $x = pdl(@{$pg{$sta}});
        printf PIPE "%f %f %f\n", $lon, $lat, median($x);
    }
    close PIPE;
#
#   Plot Pn GT residuals as triangles
#
    $pipe = "psxy $r -St0.13 -C$cpt -Wthinnest -L -O -K >> $psfil";
    open (PIPE,"| $pipe") || die "ERROR : Cannot open pipe!\n";
    for $sta (keys %pn) {
        ($lon, $lat, undef) = split ' ', $net->{$sta};
        $k = @{$pn{$sta}};
        $x = pdl(@{$pn{$sta}});
        printf PIPE "%f %f %f\n", $lon, $lat, median($x);
    }
    close PIPE;
#
#   Pb/Pg/Pn GT residuals RSTT
#
    %pb = %pg = %pn = ();
    $bx = "-Ba$xtic"."g$xtic/a$ytic"."g$ytic:.\"$name, RSTT\":wESn";
    `pscoast $r $bx $kmbar $mpar -Di -N1 -X5 -O -K >> $psfil`;
    $pipe = "psxy $r -Sa0.1 -Ggreen -Wthinnest -O -K >> $psfil";
    open (PIPE,"| $pipe") || die "ERROR : Cannot open pipe!\n";
    for $hypid (keys %$gthypo) {
        next if ($gthypo->{$hypid}[5] ne $name);
        printf PIPE "%f %f\n", $gthypo->{$hypid}[3], $gthypo->{$hypid}[2];
        for ($i = 0; $i < $m; $i++) {
            next unless ($outpicks->[$i][6] == $hypid);
            ($sta, $tres) = ($outpicks->[$i][0], $outpicks->[$i][4]);
            $phase = $outpicks->[$i][1];
            push @{$pb{$sta}}, $tres if ($phase eq "Pb");
            push @{$pg{$sta}}, $tres if ($phase eq "Pg");
            push @{$pn{$sta}}, $tres if ($phase eq "Pn");
        }
    }
    close PIPE;
#
#   Plot Pb GT residuals as circles
#
    $pipe = "psxy $r -Sc0.12 -C$cpt -Wthinnest -L -O -K >> $psfil";
    open (PIPE,"| $pipe") || die "ERROR : Cannot open pipe!\n";
    for $sta (keys %pb) {
        ($lon, $lat, undef) = split ' ', $net->{$sta};
        $k = @{$pb{$sta}};
        $x = pdl(@{$pb{$sta}});
        printf PIPE "%f %f %f\n", $lon, $lat, median($x);
    }
    close PIPE;
#
#   Plot Pg GT residuals as inverted triangles
#
    $pipe = "psxy $r -Si0.13 -C$cpt -Wthinnest -L -O -K >> $psfil";
    open (PIPE,"| $pipe") || die "ERROR : Cannot open pipe!\n";
    for $sta (keys %pg) {
        ($lon, $lat, undef) = split ' ', $net->{$sta};
        $k = @{$pg{$sta}};
        $x = pdl(@{$pg{$sta}});
        printf PIPE "%f %f %f\n", $lon, $lat, median($x);
    }
    close PIPE;
#
#   Plot Pn GT residuals as triangles
#
    $pipe = "psxy $r -St0.13 -C$cpt -Wthinnest -L -O -K >> $psfil";
    open (PIPE,"| $pipe") || die "ERROR : Cannot open pipe!\n";
    for $sta (keys %pn) {
        ($lon, $lat, undef) = split ' ', $net->{$sta};
        $k = @{$pn{$sta}};
        $x = pdl(@{$pn{$sta}});
        printf PIPE "%f %f %f\n", $lon, $lat, median($x);
    }
    close PIPE;
    `psscale -C$cpt -E -D-0.25/-0.6/9/0.12h -Ba1:\"Pg (\@~\321\@~), Pb (\@~o\@~), Pn (\@~D\@~) GT median time residual [s]\": -O >> $psfil`;

#
#   Sb/Sg/Sn GT residuals, ak135
#
    $psfil = "$name.SgSnResMap.ps";
    %pb = %pg = %pn = ();
    $bx = "-Ba$xtic"."g$xtic/a$ytic"."g$ytic:.\"$name, ak135\":WeSn";
    `pscoast $r $bx $kmbar $mpar -Di -N1 -K > $psfil`;
    $pipe = "psxy $r -Sa0.1 -Ggreen -Wthinnest -O -K >> $psfil";
    open (PIPE,"| $pipe") || die "ERROR : Cannot open pipe!\n";
    for $hypid (keys %$gthypo) {
        next if ($gthypo->{$hypid}[5] ne $name);
        printf PIPE "%f %f\n", $gthypo->{$hypid}[3], $gthypo->{$hypid}[2];
        for ($i = 0; $i < $n; $i++) {
            next unless ($inpicks->[$i][6] == $hypid);
            ($sta, $tres) = ($inpicks->[$i][0], $inpicks->[$i][4]);
            $phase = $inpicks->[$i][1];
            push @{$pb{$sta}}, $tres if ($phase eq "Sb");
            push @{$pg{$sta}}, $tres if ($phase eq "Sg" or $phase eq "Lg");
            push @{$pn{$sta}}, $tres if ($phase eq "Sn");
        }
    }
    close PIPE;
#
#   Plot Sb GT residuals as circles
#
    $pipe = "psxy $r -Sc0.12 -C$cpt -Wthinnest -L -O -K >> $psfil";
    open (PIPE,"| $pipe") || die "ERROR : Cannot open pipe!\n";
    for $sta (keys %pb) {
        ($lon, $lat, undef) = split ' ', $net->{$sta};
        $k = @{$pb{$sta}};
        $x = pdl(@{$pb{$sta}});
        printf PIPE "%f %f %f\n", $lon, $lat, median($x);
    }
    close PIPE;
#
#   Plot Sg GT residuals as inverted triangles
#
    $pipe = "psxy $r -Si0.13 -C$cpt -Wthinnest -L -O -K >> $psfil";
    open (PIPE,"| $pipe") || die "ERROR : Cannot open pipe!\n";
    for $sta (keys %pg) {
        ($lon, $lat, undef) = split ' ', $net->{$sta};
        $k = @{$pg{$sta}};
        $x = pdl(@{$pg{$sta}});
        printf PIPE "%f %f %f\n", $lon, $lat, median($x);
    }
    close PIPE;
#
#   Plot Sn GT residuals as triangles
#
    $pipe = "psxy $r -St0.13 -C$cpt -Wthinnest -L -O -K >> $psfil";
    open (PIPE,"| $pipe") || die "ERROR : Cannot open pipe!\n";
    for $sta (keys %pn) {
        ($lon, $lat, undef) = split ' ', $net->{$sta};
        $k = @{$pn{$sta}};
        $x = pdl(@{$pn{$sta}});
        printf PIPE "%f %f %f\n", $lon, $lat, median($x);
    }
    close PIPE;
#
#   Sb/Sg/Sn GT residuals RSTT
#
    %pb = %pg = %pn = ();
    $bx = "-Ba$xtic"."g$xtic/a$ytic"."g$ytic:.\"$name, RSTT\":wESn";
    `pscoast $r $bx $kmbar $mpar -Di -N1 -X5 -O -K >> $psfil`;
    $pipe = "psxy $r -Sa0.1 -Ggreen -Wthinnest -O -K >> $psfil";
    open (PIPE,"| $pipe") || die "ERROR : Cannot open pipe!\n";
    for $hypid (keys %$gthypo) {
        next if ($gthypo->{$hypid}[5] ne $name);
        printf PIPE "%f %f\n", $gthypo->{$hypid}[3], $gthypo->{$hypid}[2];
        for ($i = 0; $i < $m; $i++) {
            next unless ($outpicks->[$i][6] == $hypid);
            ($sta, $tres) = ($outpicks->[$i][0], $outpicks->[$i][4]);
            $phase = $outpicks->[$i][1];
            push @{$pb{$sta}}, $tres if ($phase eq "Sb");
            push @{$pg{$sta}}, $tres if ($phase eq "Sg" or $phase eq "Lg");
            push @{$pn{$sta}}, $tres if ($phase eq "Sn");
        }
    }
    close PIPE;
#
#   Plot Sb GT residuals as circles
#
    $pipe = "psxy $r -Sc0.12 -C$cpt -Wthinnest -L -O -K >> $psfil";
    open (PIPE,"| $pipe") || die "ERROR : Cannot open pipe!\n";
    for $sta (keys %pb) {
        ($lon, $lat, undef) = split ' ', $net->{$sta};
        $k = @{$pb{$sta}};
        $x = pdl(@{$pb{$sta}});
        printf PIPE "%f %f %f\n", $lon, $lat, median($x);
    }
    close PIPE;
#
#   Plot Sg GT residuals as inverted triangles
#
    $pipe = "psxy $r -Si0.13 -C$cpt -Wthinnest -L -O -K >> $psfil";
    open (PIPE,"| $pipe") || die "ERROR : Cannot open pipe!\n";
    for $sta (keys %pg) {
        ($lon, $lat, undef) = split ' ', $net->{$sta};
        $k = @{$pg{$sta}};
        $x = pdl(@{$pg{$sta}});
        printf PIPE "%f %f %f\n", $lon, $lat, median($x);
    }
    close PIPE;
#
#   Plot Sn GT residuals as triangles
#
    $pipe = "psxy $r -St0.13 -C$cpt -Wthinnest -L -O -K >> $psfil";
    open (PIPE,"| $pipe") || die "ERROR : Cannot open pipe!\n";
    for $sta (keys %pn) {
        ($lon, $lat, undef) = split ' ', $net->{$sta};
        $k = @{$pn{$sta}};
        $x = pdl(@{$pn{$sta}});
        printf PIPE "%f %f %f\n", $lon, $lat, median($x);
    }
    close PIPE;
    `psscale -C$cpt -E -D-0.25/-0.6/9/0.12h -Ba1:\"Sg (\@~\321\@~), Sb (\@~o\@~), Sn (\@~D\@~) GT median time residual [s]\": -O >> $psfil`;
}

#
# GT residuals vs delta
#
sub gtresiduals_delta {
    my ($name,  $inpicks, $outpicks, $gthypo, $pha) = @_;
    my ($tlo, $thi, $tic, $tic2, $dhi, $dtic, $dtic2, $tstep, $hi);
    my ($psfil, $r, $bx, $pipe, $ytit, $xtit, $i, $k, $n, $m, $hypid);
    my ($delta, $tres, $phase, @pg, @pb, @pn);
    my ($d, $t, $x);
    my ($xsize, $yshift) = ("5.5/3.2", 3.4);
    $hi = 15;
#
#   ak135 residuals vs distance
#
    @pg = @pb = @pn = ();
    $n = @$inpicks;
    ($tlo, $thi, $tic, $tic2, $tstep) = (-5, 5, 1, 0.5, 0.1);
    ($tlo, $thi, $tic, $tic2, $tstep) = (-10, 10, 2, 1, 0.1) if ($pha eq "S");
    ($dhi, $dtic, $dtic2) = (15, 1, 0.5);
    $ytit = "GT time residual [s]";
    $xtit = "Delta [degree]";
    $r = "-R0/$dhi/$tlo/$thi -JX$xsize";
    $psfil = "$name.delta.tres.$pha.ps";
    $bx = "-Ba$dtic"."g$dtic2:\"$xtit\":/a$tic"."g$tic:\"$ytit\":WeSn";
    `psbasemap $r $bx -K > $psfil`;
    for $hypid (keys %$gthypo) {
        if ($name ne "All") {
            next if ($gthypo->{$hypid}[5] ne $name);
        }
        for ($i = 0; $i < $n; $i++) {
            next unless ($inpicks->[$i][6] == $hypid);
            ($delta, $tres) = ($inpicks->[$i][2], $inpicks->[$i][4]);
            $phase = $inpicks->[$i][1];
            if ($pha eq "P") {
                push @pb, [ $delta, $tres ] if ($phase eq "Pb");
                push @pg, [ $delta, $tres ] if ($phase eq "Pg");
                push @pn, [ $delta, $tres ] if ($phase eq "Pn");
            }
            else {
                push @pb, [ $delta, $tres ] if ($phase eq "Sb");
                push @pg, [ $delta, $tres ] if ($phase eq "Sg");
                push @pn, [ $delta, $tres ] if ($phase eq "Sn");
            }
        }
    }
#
#   residuals vs distance plots
#
    $n = @pg;
    if ($n > 10) {
        $x = pdl(@pg);
        plot_dist($psfil, $r, $x, "salmon", "red");
    }
    $n = @pb;
    if ($n > 10) {
        $x = pdl(@pb);
        plot_dist($psfil, $r, $x, "green", "DarkSeaGreen");
    }
    $n = @pn;
    if ($n > 10) {
        $x = pdl(@pn);
        plot_dist($psfil, $r, $x, "DodgerBlue", "blue");
    }
    $i = "0.3 -4.8 14 0 1 BL ak135\: \@;red;Pg \@;;\@;DarkSeaGreen;Pb \@;;\@;blue;Pn \@;;";
    $i = "0.3 -9.8 14 0 1 BL ak135\: \@;red;Lg \@;;\@;DarkSeaGreen;Sb \@;;\@;blue;Sn \@;;";
    `pstext $r -O -K << END >> $psfil\n$i\nEND\n`;
#
#   histogram of residuals
#
    $bx = "-Ba$tic"."g$tic2:\"GT time residual [s]\":/a2g1:Percent:wESn";
    $r = "-R$tlo/$thi/0/$hi -JX3.2";
    `psbasemap $r $bx -X5.8 -O -K >> $psfil`;
    $n = @pg;
    if ($n > 10) {
        $x = pdl(@pg);
        $t = flat($x(1, ));
        $m = $t->dim(0);
        plot_histo($psfil, $r, "red", $m, $t, $tlo, $thi, $tstep);
    }
    $n = @pb;
    if ($n > 10) {
        $x = pdl(@pb);
        $t = flat($x(1, ));
        $m = $t->dim(0);
        plot_histo($psfil, $r, "DarkSeaGreen", $m, $t, $tlo, $thi, $tstep);
    }
    $n = @pn;
    if ($n > 10) {
        $x = pdl(@pn);
        $t = flat($x(1, ));
        $m = $t->dim(0);
        plot_histo($psfil, $r, "blue", $m, $t, $tlo, $thi, $tstep);
    }
#
#   RSTT residuals vs distance
#
    @pg = @pb = @pn = ();
    $n = @$outpicks;
    $ytit = "GT time residual [s]";
    $xtit = "Delta [degree]";
    $r = "-R0/$dhi/$tlo/$thi -JX$xsize";
    $bx = "-Ba$dtic"."g$dtic2/a$tic"."g$tic:\"$ytit\":Wesn";
    `psbasemap $r $bx -Y$yshift -X-5.8 -O -K >> $psfil`;
    for $hypid (keys %$gthypo) {
        if ($name ne "All") {
            next if ($gthypo->{$hypid}[5] ne $name);
        }
        for ($i = 0; $i < $n; $i++) {
            next unless ($outpicks->[$i][6] == $hypid);
            ($delta, $tres) = ($outpicks->[$i][2], $outpicks->[$i][4]);
            $phase = $outpicks->[$i][1];
            if ($pha eq "P") {
                push @pg, [ $delta, $tres ] if ($phase eq "Pg");
                push @pn, [ $delta, $tres ] if ($phase eq "Pn");
            }
            else {
                push @pg, [ $delta, $tres ] if ($phase eq "Sg" or $phase eq "Lg");
                push @pn, [ $delta, $tres ] if ($phase eq "Sn");
            }
        }
    }
#
#   residuals vs distance plots
#
    $n = @pg;
    if ($n > 10) {
        $x = pdl(@pg);
        plot_dist($psfil, $r, $x, "salmon", "red");
    }
    $n = @pn;
    if ($n > 10) {
        $x = pdl(@pn);
        plot_dist($psfil, $r, $x, "DodgerBlue", "blue");
    }
    $i = "RSTT\: \@;red;Pg \@;;\@;blue;Pn \@;;";
    $i = "0.3 -9.8 14 0 1 BL ak135\: \@;red;Lg \@;;\@;blue;Sn \@;;";
    `pstext $r -O -K << END >> $psfil\n$i\nEND\n`;
#
#   histogram of residuals
#
    $bx = "-Ba$tic"."g$tic2/a2g1:Percent:wEsn";
    $r = "-R$tlo/$thi/0/$hi -JX3.2";
    `psbasemap $r $bx -X5.8 -O -K >> $psfil`;
    $n = @pg;
    if ($n > 10) {
        $x = pdl(@pg);
        $t = flat($x(1, ));
        $m = $t->dim(0);
        plot_histo($psfil, $r, "red", $m, $t, $tlo, $thi, $tstep);
    }
    $n = @pn;
    if ($n > 10) {
        $x = pdl(@pn);
        $t = flat($x(1, ));
        $m = $t->dim(0);
        plot_histo($psfil, $r, "blue", $m, $t, $tlo, $thi, $tstep);
    }
    `psxy /dev/null $r -O >> $psfil`;
}


sub plot_dist {
    my ($psfil, $r, $x, $cp, $cl) = @_;
    my ($pipe, $i, $n, $m, $pct, $tlo, $thi);
    my ($d, $t, $u, $xc, $med);
    $pct = 5;
#
#   outliers: cut out upper and lower 1 percentile
#
    $tlo = min(pdl(pct($x(1,), 0.008)));
    $thi = max(pdl(pct($x(1,), 0.992)));
    $u = $x(, which($x(1,) >= $tlo & $x(1,) <= $thi))->copy;
    $d = flat($u(0,));
    $t = flat($u(1,));
    $n = $d->dim(0);
    if ($n) {
        $pipe = "psxy $r -Sx0.025 -Wthin,$cp -O -K >> $psfil";
        open (PIPE, "| $pipe") || die "ERROR : Cannot open pipe!\n";
        for ($i = 0; $i < $n; $i++) {
            printf PIPE "%f %f\n", $d->at($i), $t->at($i);
        }
        close PIPE;
        ($xc, $med, undef) = pdl_percentile_stats($d, $t, $pct);
        $m = $xc->dim(0);
        $pipe = "psxy $r -Wthick,$cl -O -K >> $psfil";
        open (PIPE, "| $pipe") || die "Cannot open pipe: $!\n";
        for ($i = 0; $i < $m - 1; $i++) {
            next if ($med->at($i) == -999);
            printf PIPE "%f %f\n", $xc->at($i), $med->at($i);
        }
        close PIPE;
    }
}

sub plot_histo {
    my ($psfil, $r, $col, $n, $t, $tlo, $thi, $tstep) = @_;
    my ($pipe, $i, $m, $med, $smad, $step, $u);
    $pipe = "pshistogram $r -W$tstep -Z1 -F -L$col -O -K >> $psfil";
    open (PIPE,"| $pipe") || die "Cannot open pipe: $!\n";
    for ($i = 0; $i < $n; $i++) {
        printf PIPE "%f\n", $t->at($i);
    }
    close PIPE;
    $med = median($t);
    $smad = 1.4826 * median(abs($t - $med));
    $step = $tstep / 10;
    ($u, undef) = hist($t, $tlo, $thi, $step);
    $m = $u->dim(0);
    $pipe = "psxy $r -Wthicker,$col -O -K >> $psfil";
    open (PIPE, "| $pipe") || die "Cannot open pipe: $!\n";
    for ($i = 0; $i < $m; $i++) {
        printf PIPE "%f %f\n", $u->at($i), 10 * gauss_pdf($u->at($i), $med, $smad);
    }
    close PIPE;
}

