#
#
# Routines for GMT4.
#
# Istvan Bondar, 2005/2/16 - 2005/2/22
#
# Science Application International Corporation
# M/S 2-1 1953 Gallows Road, Suite 260
# Vienna, VA 22182, USA
# tel: (+1) 703 610-4170,  fax: (+1) 703 610-4163
# e-mail: Istvan.K.Bondar@saic.com
#
# Ben Kohl 2005/3/9 - added routines to plot an arrival on the waveform
#       see:  optimal_pick_line(), plot_pick()
# Istvan Bondar 2005/5/18 - added optimal_pick_line_pdl().
#
#    gmtdefaults()
#        Sets default GMT settings.
#
#    ($lo, $hi, $tic) = gettic($xmin, $xmax)
#        Returns nice plot limits lo,hi and tic interval.
#
#    ($r, $kmbar, $xtic, $ytic, $west, $south, $east, $north) =
#    projbox($radius, $isize, $midlat, $midlon)
#        Calculates equidistant projection box encompassing a a circle with
#        a given radius around the center point. Returns projection box
#        (-R, -JE) and scale bar (-L) flags for GMT as well as wesn coordinates
#        and annotation/gridline intervals.
#
#    $rj = wfmbasemap($psfil, $stime, $etime, $ymin, $ymax, $xsize, $ysize,
#                     $ytitle, $title, $isnew, $istannot, $xshift, $yshift,
#                     $ylimit_set_by_user)
#        Creates basemap for a time series (waveform) plot. Returns GMT -R -J
#        parameters.
#
#    plot_waveform($psfil, $r, $tstart, $sps, $color, \@wfm)
#        Plots time series using $rj returned by wfmbasemap.
#
#    plot_waveform_pdl($psfil, $r, $tstart, $sps, $color, $wfm_pdl)
#        Plots time series piddle using $rj returned by wfmbasemap.
#
#    ($p_ymin, $p_ymax) = optimal_pick_line($ptime, $tstart, $sps, \@wfm)
#        Calculates ymin, ymax values for plotting a phase pick.
#
#    ($p_ymin, $p_ymax) = optimal_pick_line_pdl($ptime, $tstart, $sps, $wfm_pdl)
#        Calculates ymin, ymax values for plotting a phase pick.
#
#    plot_pick($psfil, $rj, $pick_time, $pick_min, $pick_max, $color,
#              $fontsize, $font, $phase)
#        Plots an arrival (flag with phase label).
#

package GMTfuncs;

require Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw(gmtdefaults gettic projbox wfmbasemap
                 plot_waveform plot_waveform_pdl
                 optimal_pick_line optimal_pick_line_pdl plot_pick);

use PDL;
use Math::Trig qw(deg2rad);
use strict;

#
#
# GMT defaults (GMT4.+)
#
#
sub gmtdefaults {
#
#   paper media
#
    `gmtset PAGE_COLOR 255/255/255 PAGE_ORIENTATION landscape`;
    `gmtset PAPER_MEDIA A4`;
#
#   basemap annotations
#
    `gmtset ANNOT_MIN_ANGLE 45 ANNOT_MIN_SPACING 0.25`;
    `gmtset ANNOT_FONT_PRIMARY Helvetica-Bold ANNOT_FONT_SIZE_PRIMARY 14p`;
    `gmtset ANNOT_FONT_SECONDARY Helvetica-Bold ANNOT_FONT_SIZE_SECONDARY 16p`;
    `gmtset ANNOT_OFFSET_PRIMARY 0.075i ANNOT_OFFSET_SECONDARY 0.075i`;
    `gmtset DEGREE_SYMBOL ring OBLIQUE_ANNOTATION 1`;
    `gmtset HEADER_FONT Helvetica-Bold HEADER_FONT_SIZE 18p`;
    `gmtset LABEL_FONT Helvetica-Bold LABEL_FONT_SIZE 16p`;
    `gmtset HEADER_OFFSET 0.12i LABEL_OFFSET 0.1125i`;
    `gmtset PLOT_CLOCK_FORMAT hh:mm:ss PLOT_DATE_FORMAT yyyy-mm-dd`;
    `gmtset PLOT_DEGREE_FORMAT D  Y_AXIS_TYPE hor_text`;
#
#   basemap layout
#
    `gmtset BASEMAP_AXES WeSn BASEMAP_FRAME_RGB 0/0/0 BASEMAP_TYPE fancy`;
    `gmtset FRAME_PEN 1.25p FRAME_WIDTH 0.05i`;
    `gmtset GRID_CROSS_SIZE_PRIMARY 0 GRID_CROSS_SIZE_SECONDARY 0`;
    `gmtset GRID_PEN_PRIMARY 0.25p GRID_PEN_SECONDARY 0.5p`;
    `gmtset MAP_SCALE_HEIGHT 0.075i POLAR_CAP 85/90`;
    `gmtset TICK_LENGTH 0.075i TICK_PEN 0.5p`;
    `gmtset X_AXIS_LENGTH 9i Y_AXIS_LENGTH 6i`;
    `gmtset X_ORIGIN 1i Y_ORIGIN 1i`;
    `gmtset UNIX_TIME FALSE UNIX_TIME_POS -075i/-0.75i`;
#
#   color system
#
    `gmtset COLOR_BACKGROUND 0/0/0 COLOR_FOREGROUND 255/255/255`;
    `gmtset COLOR_NAN 128/128/128 COLOR_IMAGE adobe COLOR_MODEL rgb`;
    `gmtset HSV_MIN_SATURATION 1 HSV_MAX_SATURATION 0.1`;
    `gmtset HSV_MIN_VALUE 0.3 HSV_MAX_VALUE 1`;
#
#   postscript
#
    `gmtset CHAR_ENCODING Standard+ DOTS_PR_INCH 300 N_COPIES 1`;
    `gmtset PS_COLOR rgb PS_IMAGE_COMPRESS none PS_IMAGE_FORMAT bin`;
    `gmtset PS_LINE_CAP butt PS_LINE_JOIN miter PS_MITER_LIMIT 0`;
    `gmtset PS_VERBOSE FALSE GLOBAL_X_SCALE 1 GLOBAL_Y_SCALE 1`;
#
#   I/O formats
#
    `gmtset D_FORMAT %lg FIELD_DELIMITER tab`;
    `gmtset GRIDFILE_SHORTHAND FALSE GRID_FORMAT nf`;
    `gmtset INPUT_CLOCK_FORMAT hh:mm:ss INPUT_DATE_FORMAT yyyy-mm-dd`;
    `gmtset IO_HEADER FALSE N_HEADER_RECS 1`;
    `gmtset OUTPUT_CLOCK_FORMAT hh:mm:ss OUTPUT_DATE_FORMAT yyyy-mm-dd`;
    `gmtset OUTPUT_DEGREE_FORMAT D XY_TOGGLE FALSE`;
#
#   projection
#
    `gmtset ELLIPSOID WGS-84 MAP_SCALE_FACTOR default MEASURE_UNIT inch`;
#
#   calendar/time
#
    `gmtset TIME_FORMAT_PRIMARY full TIME_FORMAT_SECONDARY full`;
    `gmtset TIME_IS_INTERVAL OFF TIME_INTERVAL_FRACTION 0.5`;
    `gmtset TIME_EPOCH 1970-01-01T00:00:00 TIME_SYSTEM UNIX`;
    `gmtset TIME_LANGUAGE us Y2K_OFFSET_YEAR 1950`;
    `gmtset TIME_UNIT s TIME_WEEK_START Sunday`;
#
#   miscellaneous
#
    `gmtset HISTORY FALSE INTERPOLANT akima`;
    `gmtset LINE_STEP 0.01i VECTOR_SHAPE 1 VERBOSE FALSE`;
}


#
#
# gettic
#   Returns tic interval and adjusted limits for an input range. The input
#   limits are adjusted so that the output limits are nice round numbers
#   and they can be divided into three or five round intervals.
#   Input:  input limits (min, max)
#   Output: output limits (min, max) and tic interval for annotation
#   Usage:  ($xlo, $xhi, $xtic) = gettic($xmin, $xmax);
#
#
sub gettic {
    my ($xmin, $xmax) = @_;
    my ($x, $y, $z, $u, $tic, $lo, $hi);
    $x = ($xmax - $xmin) / 3;
    $y = sclr(pow(10.0, floor(log10($x))));
    $z = $x / $y;
    $u = 2;
    $u = 1 if ($z < 1.5);
    $u = 5 if ($z > 3.5);
    $tic = $u * $y;
    $lo = sclr(floor($xmin / $tic) * $tic);
    $hi = sclr(ceil ($xmax / $tic) * $tic);
    return ($lo, $hi, $tic);
}


#
#
# projbox
#   Calculates equidistant projection box encompassing a a circle with
#   a given radius around the center point.
#   Input:  radius around center point in degrees
#           size of plot in inches
#           center's latitude
#           center's longitude
#   Output: projection parameters in GMT "-R -J" format
#           km scale bar in GMT "-L" format
#           longitude, latitude annotation intervals
#           lower left and upper right coordinates of bounding box
#   Usage:
#           ($r, $kmbar, $xtic, $ytic, $west, $south, $east, $north) =
#                   projbox($radius, $isize, $clat, $clon);
#           `pscoast $r -Ba$xtic/a$ytic $kmbar -Givory -Wthin -Dc > $psfil`;
#
#
sub projbox {
    my ($radius, $isize, $midlat, $midlon) = @_;
    my ($scale, $kmbar, $barlength, $r, $xtic, $ytic);
    my ($top, $bot, $left, $right, $ds, $x, $y, $z, $ispole);
    my ($ll_x, $ll_y, $ur_x, $ur_y, $west, $south, $east, $north);
#
#   scale
#
    $scale = $radius * 11112000. / ($isize * 2.54);
#
#   compute rough number of degrees north/south, east/west at the current scale
#
    $ds = 2.54 * $scale / 11112000.;
    $top = $midlat + $isize * $ds;
    $bot = $midlat - $isize * $ds;
    if (abs($midlat) < 90.) {
        $left = $midlon - $isize * $ds / cos(deg2rad($midlat));
        $right = $midlon + $isize * $ds / cos(deg2rad($midlat));
        $ispole = 0;
    }
    else {
        $left = -180.;
        $right = 180.;
        $ispole = 1;
    }
#
#   sanity checks
#
    if ($right - $left > 360) {
        $left = $midlon - 180;
        $right = $midlon + 180;
    }
    $top = 90  if ($top > 90);
    $bot = -90 if ($bot < -90);
#
#   rough plot range  based on scale and plot size
#
    $r = sprintf "-R%.3f/%.3f/%.3f/%.3fr -Je%.3f/%.3f/1:%.0f",
                 $left, $bot, $right, $top, $midlon, $midlat, $scale;
#
#   compute the inches of the midlon, midlat point
#
    chomp($z = `echo $midlon $midlat | mapproject $r -Di`);
    ($x, $y) = split ' ',$z;
#
#   use the xy-inches of the midpoint to inverse project the inches to
#   get a rectangular plot of a particular size
#
    $ll_x = $x - $isize;
    $ll_y = $y - $isize;
    $ur_x = $x + $isize;
    $ur_y = $y + $isize;
#
#   after normalizing to get the center point, compute the lower-left and
#   upper-right points, guarentees a nice square plot.
#
    chomp($z = `echo $ll_x $ll_y | mapproject $r -Di -I`);
    ($west, $south) = split ' ',$z;
    chomp($z = `echo $ur_x $ur_y | mapproject $r -Di -I`);
    ($east, $north) = split ' ',$z;

    $east += 360 if ($east < $west);
#
#   projection
#
    $r = sprintf "-R%.3f/%.3f/%.3f/%.3fr -JE%.3f/%.3f/%.2f",
                  $west, $south, $east, $north, $midlon, $midlat, $isize;
#
#   gridline intervals
#
    $xtic = gettic($west, $east);
    if ($ispole) {
        $ytic = 5.;
    }
    else {
        $ytic = gettic($south, $north);
    }
#
#   scalebar
#
    $barlength = $isize * 0.5;
    $barlength *= 2.54 * $scale / 100000;
    $barlength = gettic(0., 3. * $barlength);
    $barlength += 1 if ($barlength < 1);
    $kmbar = sprintf "-Lfx%.3f/%.3f/%.3f/%d", $isize / 4, $isize / 8,
                                              $midlat, int($barlength);
#
#   return projection, kmbar flags, gridline intervals and coordinates
#   of lower left and upper right corners
#
    return ($r, $kmbar, $xtic, $ytic, $west, $south, $east, $north);
}

#
#
# basemap for waveform plot
#   Creates basemap for waveform plots
#   Notes:  time axis is annotated in accordance with the PLOT_DATE_FORMAT and
#           PLOT_CLOCK_FORMAT settings. Use gmtset to change the defaults.
#   Input:  postscript file name
#           time and amplitude limits
#           X, Y size (inches)
#           y title
#           plot title
#           new plot flag (if true, uses -K >  otherwise -O -K >> )
#           time annotation flag (if true time axis will be annotated)
#           X, Y shifts
#           amplitude limits are set by user (0/1)
#   Output: GMT4 -R, -J parameters
#   Usage:
#           $rj = wfmbasemap($psfil, $stime, $etime, $ylo, $yhi, $xsize, $ysize,
#                            "PSZ bz [nm/s]", $title, 0, 1, 0, $yshift, $yuser);
#
sub wfmbasemap {
    my ($psfil, $stime, $etime, $ymin, $ymax, $xsize, $ysize,
        $ytitle, $title, $isnew, $istannot, $xshift, $yshift,
        $ylimit_set_by_user) = @_;
    my ($ylo, $yhi, $ytic, $tlo, $thi, $t, $r, $s, $bx);
    ($ylo, $yhi, $ytic) = gettic($ymin, $ymax);
    ($ylo, $yhi) = ($ymin, $ymax) if ($ylimit_set_by_user);
#
#   -R -J for epoch times
#
    $tlo = sprintf "%.0f", $stime;
    $thi = sprintf "%.0f", $etime;
    $t = $thi - $tlo;
    $r = "-R$tlo/$thi/$ylo/$yhi -JX$xsize"."t/$ysize";
#
#   time axis annotation
#
    if ($istannot) {
        if   ($t > 180 * 86400) { $bx = "-Bsa1Y/0WeSn -Bpa3of1o";   }
        elsif ($t > 30 * 86400) { $bx = "-Bsa2O/0WeSn -Bpa30df10d"; }
        elsif ($t >  7 * 86400) { $bx = "-Bsa7D/0WeSn -Bpa2df1d";   }
        elsif ($t >      86400) { $bx = "-Bsa1D/0WeSn -Bpa12hf6h";  }
        elsif ($t >   6 * 3600) { $bx = "-Bsa6H/0WeSn -Bpa3hf1h";   }
        elsif ($t >       3600) { $bx = "-Bsa1H/0WeSn -Bpa30mf10m"; }
        elsif ($t >       1800) { $bx = "-Bsa30M/0WeSn -Bpa10mf5m"; }
        elsif ($t >        600) { $bx = "-Bsa10M/0WeSn -Bpa5mf1m";  }
        elsif ($t >        300) { $bx = "-Bsa2M/0wSen -Bpa1mf30s";  }
        elsif ($t >         60) { $bx = "-Bsa2M/0WeSn -Bpa1mf30s";  }
        elsif ($t >         30) { $bx = "-Bsa1M/0WeSn -Bpa30sf10s"; }
        elsif ($t >         10) { $bx = "-Bsa1M/0WeSn -Bpa5sf1s";   }
        else                    { $bx = "-Bsa1M/0WeSn -Bpa2sf1s";   }
    }
    else {
        if   ($t > 180 * 86400) { $bx = "-Bsf1Y/0WeSn -Bpf1o";  }
        elsif ($t > 30 * 86400) { $bx = "-Bsf2O/0WeSn -Bpf10d"; }
        elsif ($t >  7 * 86400) { $bx = "-Bsf7D/0WeSn -Bpf1d";  }
        elsif ($t >      86400) { $bx = "-Bsf1D/0WeSn -Bpf6h";  }
        elsif ($t >   6 * 3600) { $bx = "-Bsf6H/0WeSn -Bpf1h";  }
        elsif ($t >       3600) { $bx = "-Bsf1H/0WeSn -Bpf10m"; }
        elsif ($t >       1800) { $bx = "-Bsf30M/0WeSn -Bpf5m"; }
        elsif ($t >        600) { $bx = "-Bsf10M/0WeSn -Bpf1m"; }
        elsif ($t >        300) { $bx = "-Bsf2M/0wSen -Bpf30s"; }
        elsif ($t >         60) { $bx = "-Bsf2M/0WeSn -Bpf30s"; }
        elsif ($t >         30) { $bx = "-Bsf1M/0WeSn -Bpf10s"; }
        elsif ($t >         10) { $bx = "-Bsf1M/0WeSn -Bpf1s";  }
        else                    { $bx = "-Bsf1M/0WeSn -Bpf1s";  }
    }
    $bx .= "/a$ytic:\"$ytitle\"::.\"$title\":WeSn";
#
#   shift of origin
#
    $s = "";
    $s .= "-X$xshift " if ($xshift);
    $s .= "-Y$yshift " if ($yshift);
#
#   basemap
#
    if ($isnew) {
       `psbasemap $r $bx $s -K > $psfil`;
    }
    else {
       `psbasemap $r $bx $s -O -K >> $psfil`;
    }
    return $r;
}


#
#
# plot_waveform
#    Plots a waveform
#    Notes: waveform is written in append mode (GMT's -O -K).
#           No gaps are assumed.
#    Input: postscript file name
#           GMT range and projection flags
#           start time of waveform
#           samples per second
#           color
#           waveform
#   Usage:
#           $rj = wfmbasemap($psfil, $stime, $etime, $ylo, $yhi, $xsize, $ysize,
#                            $xtit, "", 1, 1, 0, 0, $isyuser);
#           plot_waveform($psfil, $rj, $tstart, $sps, $color{$chan}, \@wfm);
#           ... (set wfm, chan, ylo, yhi, xtit) ...
#           $rj = wfmbasemap($psfil, $stime, $etime, $ylo, $yhi, $xsize, $ysize,
#                            $xtit, $title, 0, 0, 0, $ysize + 0.2, $isyuser);
#           plot_waveform($psfil, $rj, $tstart, $sps, $color{$chan}, \@wfm);
#           `psxy /dev/null $rj -O >> $psfil`;
#
sub plot_waveform {
    my ($psfil, $r, $tstart, $sps, $color, $wfm) = @_;
    my ($t, $n, $i, $pipe);
    $n = @$wfm;
    $sps = 1 / $sps;
    $pipe = "psxy $r -Wthin/$color -O -K >> $psfil";
    open (PIPE,"| $pipe") || die "ERROR : Cannot open pipe!\n";
    for ($i = 0, $t = $tstart; $i < $n; $i++, $t += $sps) {
        print PIPE "$t $wfm->[$i]\n";
    }
    close PIPE;
}


#
#
# plot_waveform_pdl
#    Plots a waveform piddle
#    Notes: waveform is written in append mode (GMT's -O -K).
#           No gaps are assumed.
#    Input: postscript file name
#           GMT range and projection flags
#           start time of waveform
#           samples per second
#           color
#           waveform piddle
#   Usage:
#           $rj = wfmbasemap($psfil, $stime, $etime, $ylo, $yhi, $xsize, $ysize,
#                            $xtit, "", 1, 1, 0, 0, $isyuser);
#           plot_waveform_pdl($psfil, $rj, $tstart, $sps, $color{$chan}, $wfm);
#           ... (set wfm, chan, ylo, yhi, xtit) ...
#           $rj = wfmbasemap($psfil, $stime, $etime, $ylo, $yhi, $xsize, $ysize,
#                            $xtit, $title, 0, 0, 0, $ysize + 0.2, $isyuser);
#           plot_waveform_pdl($psfil, $rj, $tstart, $sps, $color{$chan}, $wfm);
#           `psxy /dev/null $rj -O >> $psfil`;
#
sub plot_waveform_pdl {
    my ($psfil, $r, $tstart, $sps, $color, $wfm) = @_;
    my ($t, $n, $i, $pipe);
    $n = $wfm->dim(0);
    $sps = 1 / $sps;
    $pipe = "psxy $r -Wthin/$color -O -K >> $psfil";
    open (PIPE,"| $pipe") || die "ERROR : Cannot open pipe!\n";
    for ($i = 0, $t = $tstart; $i < $n; $i++, $t += $sps) {
        print PIPE "$t ", $wfm->at($i), "\n";
    }
    close PIPE;
}


#
#
# Gets lower and upper y coordinates of a phase pick flag (PDL version)
#     Given a waveform piddle and a pick this routine computes the ymin, ymax
#     values for plotting the pick
#
#     Basic logic is, starts at the waveform, and draws upward an amount of
#            0.25 * ( $wfmax - $wfmin );
#
#    if it doesn't overlap it draws from
#           ( $wfmax + $wfmin ) / 2.;
#
sub optimal_pick_line_pdl
{
    my ($ptime, $wftime, $sps, $wfm) = @_;
    my ($pdl, $i, $n, $wfmin, $wfmax, $pick_min, $pick_max);
    $n = $wfm->dim(0);
    ($wfmin, $wfmax) = minmax($wfm);

    $i = int(($ptime - $wftime) * $sps);
    if ($i < 0 || $i > $n - 2) {
        $pick_min = ($wfmax + $wfmin) / 2.;
        $pick_max = $pick_min + 0.25 * ($wfmax - $wfmin);
        return($pick_min, $pick_max);
    }

    $pick_min = $wfm->at($i + 1);
    $pick_min = $wfm->at($i) if ($wfm->at($i) > $wfm->at($i+1));
    $pick_max = $pick_min + 0.25 * ($wfmax - $wfmin);

    return($pick_min, $pick_max);
}


#
#
# Gets lower and upper y coordinates of a phase pick flag
#     Given a waveform and a pick (which presumable overlays the waveform)
#     this routine computes the ymin, ymax values for plotting the pick
#
#     Basic logic is, starts at the waveform, and draws upward an amount of
#            0.25 * ( $wfmax - $wfmin );
#
#    if it doesn't overlap it draws from
#           ( $wfmax + $wfmin ) / 2.;
#
sub optimal_pick_line
{
    my ($ptime, $wftime, $sps, $wf_ref) = @_;
    my ($pdl, $i, $wfmin, $wfmax, $pick_min, $pick_max);

    $pdl = pdl(@{$wf_ref});
    ($wfmin, $wfmax) = minmax($pdl);

    $i = int(($ptime - $wftime) * $sps);
    if ($i < 0 || $i > $#$wf_ref - 1) {
        $pick_min = ($wfmax + $wfmin) / 2.;
        $pick_max = $pick_min + 0.25 * ($wfmax - $wfmin);
        return($pick_min, $pick_max);
    }

    if ($wf_ref->[$i] > $wf_ref->[$i+1]) {
        $pick_min = $wf_ref->[$i];
    } else {
        $pick_min = $wf_ref->[$i+1];
    }
    $pick_max = $pick_min + 0.25 * ($wfmax - $wfmin);

    return($pick_min, $pick_max);
}


#
#
# plots an arrival
#    Notes: arrival is drawn in append mode (GMT's -O -K).
#           No gaps are assumed.
#    Input: postscript file name
#           GMT range and projection flags
#           time of pick
#           min and max values in amplitude units for drawing pick line
#           color
#           font size
#           font id
#           phase name
#   Usage:
#           $rj = wfmbasemap($psfil, $stime, $etime, $ylo, $yhi, $xsize, $ysize,
#                            $xtit, "", 1, 1, 0, 0, $isyuser);
#           plot_waveform_pdl($psfil, $rj, $tstart, $sps, $color{$chan}, $wfm);
#              ($pick_min, $pick_max) = optimal_pick_line_pdl($pick_time, $tstart,
#                                                         $sps, $wfm );
#              plot_pick($psfil, $rj, $pick_time, $pick_min, $pick_max, $color,
#                    $fontsize, $font, $phase);
#           `psxy /dev/null $rj -O >> $psfil`;
#
sub plot_pick {
    my ($psfil, $r, $tstart, $ymin, $ymax, $color, $fontsize, $font, $text) = @_;
    my ($tmin, $tmax, $s);
    foreach my $element_of_r (split(/\s+/,$r)) {
        next unless ( $element_of_r =~ /^-R/ );
        $element_of_r =~ s/^-R//;
        if ( $element_of_r =~ /r$/ ) {
            $element_of_r =~ s/r$//;
            ($tmin,$tmax) = (split(/\//,$element_of_r))[0,2];
        } else {
            ($tmin,$tmax) = (split(/\//,$element_of_r))[0,1];
        }
    }
    unless ( defined($tmin) ) {
        print STDERR "Error: can not determine tmin,tmax from: $r\n";
        return;
    }

    if ( $tstart < $tmin || $tstart > $tmax ) { return; }

    $s = "$tstart $ymin\n$tstart $ymax";
    `psxy $r -Wthin/$color -O -K << END >> $psfil\n$s\nEND\n`;
    $s = "$tstart $ymax $fontsize 0 $font BL $text";
    `pstext $r -W$color -Gwhite -O -K << END >> $psfil\n$s\nEND\n`;
}


1;

