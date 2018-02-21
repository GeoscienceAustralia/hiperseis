#!/usr/bin/perl
#
#
$|=1;
use strict;
use DBI;
use Connect;
use MYSQLselect;
use Gregion;

my ($dbh, $eventids, $eventid, $evid, $orid, $mid, $fid, $fname, $fe);
my ($etype, $grname, $ecom, $block, $mblock, $key, $n, $m, $ht, $ilocroot);
my ($startday, $endday, $db, $outfile);

if (@ARGV < 3) {
   print "  Usage: GetISF2Bulletin.pl startday endday outfile\n";
   print "Example: ./GetISF2Bulletin.pl 2015-11-01 2015-12-01 out.isf\n";
   exit -1;
}
($startday, $endday, $outfile) = @ARGV;

#
# establish DB connection
#
$dbh = Connect::mysql_logon("seiscomp3", "sysop", "sysop", "localhost");

#
# get FlinnEngdahl regionalization file
#
$ilocroot = $ENV{ILOCROOT};
$fname = "$ilocroot/etc/FlinnEngdahl/FE.dat";
if ($ilocroot) {
    $fname = "$ilocroot/etc/FlinnEngdahl/FE.dat";
}
else {
    print "Warning: ILOCROOT env variable not set! ";
    print "Looking for FE.dat in the current directory instead.\n";
}
$fe = read_FlinnEngdahl($fname);

#
# Get event list
#
($n, $eventids) = get_event_list($startday, $endday, $db, $dbh);
print "$n events between $startday and $endday\n";
$n = 0;
open(OU, "> $outfile") || die "cannot open $outfile!\n";
for $eventid (@$eventids) {
    print "    $eventid...\n";
#
#   Get event from DB
#
    ($m, $evid, $orid, $mid, $fid, $ht, $etype,
     $grname, $ecom) = get_event($eventid, $fe, $db, $dbh);
    next if ($m);
#
#   open output ISF2 file
#
    print OU "BEGIN IMS2.0\n";
    print OU "DATA_TYPE BULLETIN IMS1.0:short with ISF2.0 extensions\n";
    printf OU "\nEvent %s %-s\n\n", $eventid, $grname;
#
#   Get origins and magnitudes from DB
#
    ($block, $mblock) = get_origins($evid, $orid, $mid, $fid,
                                    $etype, $db, $dbh);
    print OU "   Date       Time        Err   RMS Latitude Longitude  ";
    print OU "Smaj  Smin  Az Depth   Err Ndef Nsta Gap  mdist  Mdist ";
    print OU "Qual   Author      OrigID    Rep   DPdep   Err Sgp\n";
    print OU "$block";
    for $key (keys %$ecom) {
        printf OU " ($key: %s)\n", $ecom->{$key};
    }
    if ($mid ne "X") {
        print OU "\nMagnitude  Err Nsta Author      OrigID\n";
        print OU "$mblock";
    }
#
#   Get phases and amplitudes from DB
#
    $block = get_phases($orid, $ht, $db, $dbh);
    print OU "\nSta     Dist  EvAz Phase        Time      TRes  ";
    print OU "Azim AzRes   Slow   SRes Def   SNR       Amp   Per ";
    print OU "Qual Magnitude    ArrID    Agy   Deploy   Ln Auth  ";
    print OU "Rep   PCh ACh L\n";
    print OU "$block\n";

    print OU "\nSTOP\n";
    $n++;
}
close OU;
$dbh->Connect::mysql_logoff;
print "$n events between $startday and $endday\n";
print "done\n";


#
# Get event list
#
sub get_event_list {
    my ($startday, $endday, $db, $dbh) = @_;
    my ($query, @eventids);
    $query = "select ep.publicID
                from Event e, Origin o, PublicObject ep, PublicObject op
               where ep._oid = e._oid
                 and o.time_value between '$startday' and '$endday'
                 and op.publicID = e.preferredOriginID
                 and op._oid = o._oid
               order by o.time_value";
    chomp(@eventids = mysql_select($dbh, $query));
    $n = @eventids;
    return $n, \@eventids;
}

#
# Get evid, orid, etype, lat, lon, depth of the preferred origin
#
sub get_event {
    my ($eventid, $fe, $db, $dbh) = @_;
    my ($query, @row, $etype, $gregname, $grn);
    my ($evid, $orid, $tc, $et, $lat, $lon, $depth, $date, $time, $ms);
    my ($dtp, $dtx, $mid, $fid, $ht, $tim);
    my %comments = ();
#
#   Get evid, orid, etype, lat, lon, depth of the preferred origin
#
    $query = "select e._oid, o._oid,
                     coalesce(e.typeCertainty, 'k'),
                     coalesce(replace(e.type, ' ', '_'), 'Q'),
                     o.latitude_value, o.longitude_value, o.depth_value,
                     o.time_value, coalesce(time_value_ms, 0),
                     coalesce(e.preferredMagnitudeID, 'X'),
                     coalesce(e.preferredFocalMechanismID, 'X')
                from Event e, Origin o, PublicObject ep, PublicObject op
               where ep.publicID = '$eventid'
                 and ep._oid = e._oid
                 and op.publicID = e.preferredOriginID
                 and op._oid = o._oid";
    chomp(@row = mysql_select($dbh, $query));
    for (@row) {
        ($evid, $orid, $tc, $et, $lat, $lon, $depth, $date, $tim, $ms,
         $mid, $fid) = split;
    }
#
#   event type certainty
#
    $etype = set_etype($tc, $et);
#    if ($etype eq "!!") {
#        print "$eventid: event type is not existing or not locatable ";
#        print "Possibly bogus event, exiting.\n";
#        return 1, $evid, $orid, $mid, $fid, $ht, $etype, $grname, \%comments;
#    }
#
#   Check for valid depth
#
    if ($depth > 700.) {
        printf "$eventid: depth (%.2f) is larger than 700 km! ", $depth;
        print "Possibly bogus event, skipping event.\n";
        return 1, $evid, $orid, $mid, $fid, $ht, $etype, $grname, \%comments;
    }
#
#   event description
#
    $gregname = "";
    $query = "select type, ';', text
                from EventDescription
               where _parent_oid = $evid";
    chomp(@row = mysql_select($dbh, $query));
    for $grn (@row) {
        ($dtp, $dtx) = split ';', $grn;
        if ($dtp =~ /region name/) {
            $gregname = $dtx;
        }
        else {
            $comments{$dtp} = $dtx;
        }
    }
    if ($gregname eq "") {
#
#       get FlinnEngdahl region name
#
        if (not exists $fe->{error}) {
            $grn = gregnum($lat, $lon, $fe);
            $gregname = gregion($grn);
        }
    }
    $time = sprintf("%s.%03d", $tim, $ms / 1000);
    $ht = "$date $time";
    $ht =~ s/\-/\//g;
    return 0, $evid, $orid, $mid, $fid, $ht, $etype, $gregname, \%comments;
}

#
# Get origins, magnitudes and focal mechanisms from DB
#
sub get_origins {
    my ($evid, $orid, $mid, $fid, $etype, $db, $dbh) = @_;
    my ($query, $p, $s, $id, $auth, $prevauth);
    my (@row, @cow, $line, $mine, $mblock, $block);
    $block = $mblock = "";
#
#   get preferred origin
#
    $query = "select time_value,
                     coalesce(time_value_ms, 0),
                     coalesce(timeFixed, -1),
                     coalesce(time_uncertainty, -1),
                     coalesce(quality_standardError, -1),
                     latitude_value,
                     longitude_value,
                     coalesce(epicenterFixed, -1),
                     coalesce(uncertainty_maxHorizontalUncertainty, -1),
                     coalesce(uncertainty_minHorizontalUncertainty, -1),
                     coalesce(uncertainty_azimuthMaxHorizontalUncertainty, -1),
                     coalesce(depth_value, -1),
                     coalesce(replace(depthType, ' ','_'),'Q'),
                     coalesce(depth_uncertainty, -1),
                     coalesce(quality_usedPhaseCount, -1),
                     coalesce(quality_usedStationCount, -1),
                     coalesce(quality_azimuthalGap, -1),
                     coalesce(quality_minimumDistance, -1),
                     coalesce(quality_maximumDistance, -1),
                     coalesce(evaluationMode, 'Q'),
                     coalesce(type, '$etype'),
                     creationInfo_agencyID,
                     _oid,
                     creationInfo_author,
                     coalesce(quality_secondaryAzimuthalGap, 360)
                from Origin
               where _oid = $orid";
    chomp(@row = mysql_select($dbh, $query));
    for $line (@row) {
        ($p, $id, $auth) = originline($line);
    }
#
#   get the rest of origins
#
    $query = "select o.time_value,
                     coalesce(o.time_value_ms, 0),
                     coalesce(o.timeFixed, -1),
                     coalesce(o.time_uncertainty, -1),
                     coalesce(o.quality_standardError, -1),
                     o.latitude_value, o.longitude_value,
                     coalesce(o.epicenterFixed, -1),
                     coalesce(o.uncertainty_maxHorizontalUncertainty, -1),
                     coalesce(o.uncertainty_minHorizontalUncertainty, -1),
                     coalesce(o.uncertainty_azimuthMaxHorizontalUncertainty, -1),
                     coalesce(o.depth_value, -1),
                     coalesce(replace(o.depthType, ' ','_'),'Q'),
                     coalesce(o.depth_uncertainty, -1),
                     coalesce(o.quality_usedPhaseCount, -1),
                     coalesce(o.quality_usedStationCount, -1),
                     coalesce(o.quality_azimuthalGap, -1),
                     coalesce(o.quality_minimumDistance, -1),
                     coalesce(o.quality_maximumDistance, -1),
                     coalesce(o.evaluationMode, 'Q'),
                     coalesce(o.type, '$etype'),
                     o.creationInfo_agencyID,
                     o._oid,
                     o.creationInfo_author,
                     coalesce(o.quality_secondaryAzimuthalGap, 360)
                from Origin o, OriginReference r, PublicObject rp
               where r._parent_oid = $evid
                 and rp.publicID = r.OriginID
                 and rp._oid = o._oid
                 and o._oid != $orid
                 and o.creationInfo_agencyID != '$auth'
            order by o.creationInfo_agencyID,
                     o.creationInfo_creationTime desc";
    chomp(@row = mysql_select($dbh, $query));
    $prevauth = "";
    for $line (@row) {
        ($s, $id, $auth) = originline($line);
        next if ($auth eq $prevauth);
        $block .= $s;
#
#       get focal mechanisms, moment tensors, comments if any
#
#       magnitude block
#
        $query = "select coalesce(nullif(type, ''), 'Q'),
                         coalesce(magnitude_value, -9),
                         coalesce(magnitude_uncertainty, -1),
                         coalesce(stationCount, -1),
                         creationInfo_agencyID
                    from Magnitude
                   where _parent_oid = $id
                   order by type";
        chomp(@cow = mysql_select($dbh, $query));
        for $mine (@cow) {
            $mblock .= magnitudeline($mine, $id);
        }
        $prevauth = $auth;
    }
#
#   preferred origin: get focal mechanisms, moment tensors, comments if any
#
    $block .= $p;
    $block .= " (#PRIME)\n";
#
#   magnitude
#
    $query = "select coalesce(nullif(type, ''), 'Q'),
                     coalesce(magnitude_value, -9),
                     coalesce(magnitude_uncertainty, -1),
                     coalesce(stationCount, -1),
                     creationInfo_agencyID
                from Magnitude m, PublicObject p
               where m._parent_oid = $orid
                 and m._oid = p._oid
                 and p.publicID != '$mid'
               order by type";
    chomp(@cow = mysql_select($dbh, $query));
    for $mine (@cow) {
        $mblock .= magnitudeline($mine, $orid);
    }
#
#   preferred magnitude
#
    $query = "select coalesce(nullif(type, ''), 'Q'),
                     coalesce(magnitude_value, -9),
                     coalesce(magnitude_uncertainty, -1),
                     coalesce(stationCount, -1),
                     creationInfo_agencyID
                from Magnitude m, PublicObject p
               where m._oid = p._oid
                 and p.publicID = '$mid'
               order by type";
    chomp(@cow = mysql_select($dbh, $query));
    for $mine (@cow) {
        $mblock .= magnitudeline($mine, $orid);
    }
    return $block, $mblock;
}

#
# get phase data from DB
#
sub get_phases {
    my ($orid, $ot, $db, $dbh) = @_;
    my ($query, $pickid, @row, @cow, @sow, $line, $block, $s, $stamag, $mdef);
    my ($psta, $i, $k, $n, $m, $sta, $snr, $type, $amp, $per, $ach, $ampmag);
    $block = "";
    $psta = "";
    $i = 0;
#
#   get phases associated to the preferred origin
#
    $query = "select r.pickID,
                     p.waveformID_stationCode,
                     coalesce(r.distance, -1),
                     coalesce(r.azimuth, -1),
                     coalesce(r.phase_code, 'Q'),
                     coalesce(p.time_value, '1970-01-01 00:00:00'),
                     coalesce(p.time_value_ms, 0),
                     coalesce(r.timeResidual, -999),
                     coalesce(p.backazimuth_value, -1),
                     coalesce(r.backazimuthResidual, -999),
                     coalesce(p.horizontalSlowness_value, -1),
                     coalesce(r.horizontalSlownessResidual, -999),
                     coalesce(r.timeUsed, r.weight > 0.5, 0),
                     coalesce(r.backazimuthUsed, 0),
                     coalesce(r.horizontalSlownessUsed, 0),
                     coalesce(nullif(p.polarity, ''), '_'),
                     coalesce(nullif(p.onset, ''), '_'),
                     coalesce(p.waveformID_networkCode, 'IR'),
                     coalesce(nullif(p.waveformID_locationCode, ''), '--'),
                     coalesce(p.waveformID_channelCode, '\?\?\?'),
                     coalesce(r.creationInfo_author, 'BUD'),
                     coalesce(p.creationInfo_agencyID, 'BUD'),
                     r._oid
                from Arrival r, Pick p, PublicObject rp
               where r._parent_oid = $orid
                 and rp.publicID = r.pickID
                 and rp._oid = p._oid
               order by r.distance, p.waveformID_stationCode, 
                     p.time_value, p.time_value_ms";
    chomp(@row = mysql_select($dbh, $query));
    for $line (@row) {
        ($pickid, $sta, undef) = split ' ', $line;
#
#       comment lines for station magnitudes
#
        if ($sta ne $psta and $i) {
            $block .= $s;
            $s = "";
        }
        $snr = $amp = $per = $ampmag = -1;
        $type = "Q";
        $ach = "???";
#
#       get snr
#
        $query = "select amplitude_value
                    from Amplitude
                   where pickID = '$pickid'
                     and type = 'snr'
                     and amplitude_value is not null";
        chomp(@cow = mysql_select($dbh, $query));
        for (@cow) {
            ($snr, undef) = split;
        }
#
#       get amplitudes
#
        $k = 0;
        $query = "select coalesce(nullif(type, ''), 'Q'),
                         coalesce(amplitude_value, -1),
                         coalesce(period_value, -1),
                         coalesce(waveformID_channelCode, '\?\?\?')
                    from Amplitude
                   where pickID = '$pickid'
                     and type != 'snr'
                     and type != 'Mwp'
                     and amplitude_value is not null";
        chomp(@cow = mysql_select($dbh, $query));
        $n = @cow;
        if ($n) {
            for (@cow) {
                ($type, $amp, $per, $ach) = split;
#
#               get ampmags, if any
#
                $query = "select coalesce(m.magnitude_value, -1),
                                 coalesce(a.snr, -1)
                            from Amplitude a,
                                 StationMagnitude m
                           where m._parent_oid = $orid
                             and a.pickID = '$pickid'
                             and a.type = m.type
                             and a.type = '$type'
                             and m.amplitudeID like a.pickID||'%'";
                chomp(@sow = mysql_select($dbh, $query));
                for (@sow) {
                    ($ampmag, $snr) = split;
                }
                $ampmag = -1 if ($ampmag > 9999998);
                $block .= phaseline($k, $ot, $line, $snr, $type, $amp, $per,
                                    $ach, $ampmag);
                if ($k == 0) {
#
#                   get station mags
#
                    $s = "";
                    $query = "select coalesce(nullif(m.type, ''), 'Q'),
                                     coalesce(m.magnitude_value, -1),
                                     coalesce(c.weight, 0)
                                from Magnitude n,
                                     StationMagnitude m,
                                     StationMagnitudeContribution c,
                                     PublicObject cp
                               where n._parent_oid = $orid
                                 and c._parent_oid = n._oid
                                 and cp.PublicID = c.stationMagnitudeID
                                 and m._oid = cp._oid
                                 and m.type is not null
                                 and c.stationMagnitudeID like '%staMag%$sta%'
                               order by m.type",
                    chomp(@sow = mysql_select($dbh, $query));
                    for (@sow) {
                        ($type, $stamag, $mdef) = split;
                        $type = "" if ($type eq "Q");
                        $type =~ s/\(/\_/;
                        $type =~ s/\)//;
                        $s .= sprintf " (Station %-5s: %4.1f $mdef\n",
                                       $type, $stamag;
                    }
                }
                $k++;
            }
        }
        else {
            $block .= phaseline($k, $ot, $line, $snr, $type, $amp, $per,
                                $ach, $ampmag);
        }
        $psta = $sta;
        $i++;
    }
    return $block;
}

#
# set IASPEI event type
#
sub set_etype {
    my ($tc, $et) = @_;
    my ($etype);
#
#   event type certainty
#
    if ($tc eq "suspected") {
        $etype = "s";
    }
    else {
        $etype = "k";
    }
#
#   event type
#
    if ($et eq "nuclear_explosion") {
        $etype .= "n";
    }
    elsif ($et =~ /explosion/ or
           $et eq "quarry_blast" or
           $et eq "rock_burst" or
           $et eq "kx") {
        $etype .= "x";
    }
    elsif ($et =~ /induced/) {
        $etype .= "i";
    }
    elsif ($et eq "landslide") {
        $etype = "ls";
    }
    elsif ($et eq "not_locatable") {
        $etype = "uk";
    }
    elsif ($et eq "not_existing") {
        $etype = "uk";
    }
    else {
        $etype .= "e";
    }
    return $etype;
}

#
# Origin line in ISF2
#
sub originline {
    my $line = shift;
    my ($date, $time, $ms, $tf, $stim, $sdobs, $lat, $lon, $ef, $s, $q);
    my ($smajax, $sminax, $strk, $depth, $dtype, $sdep, $ndef, $nsta);
    my ($gap, $mind, $maxd, $em, $etp, $author, $id, $reporter, $sgap);
    ($date, $time, $ms, $tf, $stim, $sdobs, $lat, $lon, $ef, $smajax, $sminax,
     $strk, $depth, $dtype, $sdep, $ndef, $nsta, $gap, $mind, $maxd,
     $em, $etp, $author, $id, $reporter, $sgap) = split ' ', $line;
    $q = ($tf < 0) ? " " : "f";
    $date =~ s/\-/\//g; 
    $s = sprintf "%s %s.%02d%s ", $date, $time, $ms / 10000, $q;
    if ($stim < 0) {
        $s .= "      ";
    }
    elsif ($stim > 99.99) {
       $s .= "99.99 "
    }
    else {
       $s .= sprintf "%5.2f ", $stim;
    }
    if ($sdobs < 0) {
        $s .= "      ";
    }
    elsif ($sdobs > 99.99) {
       $s .= "99.99 "
    }
    else {
       $s .= sprintf "%5.2f ", $sdobs;
    }
    $q = ($ef < 0) ? " " : "f";
    $s .= sprintf "%8.4f %9.4f%s", $lat, $lon, $q;
    if ($smajax < 0) {
        $s .= "      ";
    }
    elsif ($smajax > 999.9) {
       $s .= "999.9 "
    }
    else {
       $s .= sprintf "%5.1f ", $smajax;
    }
    if ($sminax < 0) {
        $s .= "      ";
    }
    elsif ($sminax > 999.9) {
       $s .= "999.9 "
    }
    else {
       $s .= sprintf "%5.1f ", $sminax;
    }
    if ($strk < 0) {
        $s .= "    ";
    }
    else {
       $s .= sprintf "%3d ", $strk;
    }
    $q = $dtype;
    $q = " " if ($dtype eq "Q" or $dtype eq "NULL");
    $q = "A" if ($dtype =~ /operator/);
    $q = "R" if ($dtype =~ /location/);
    if ($depth < 0) {
        $s .= sprintf "     %s ", $q;
    }
    else {
        $s .= sprintf "%5.1f%s ", $depth, $q;
    }
    if ($sdep < 0) {
        $s .= "     ";
    }
    else {
        $sdep = 99.9 if ($sdep > 99.9);
        $s .= sprintf "%4.1f ", $sdep;
    }
    if ($ndef < 0) {
        $s .= "     ";
    }
    else {
        $s .= sprintf "%4d ", $ndef;
    }
    if ($nsta < 0) {
        $s .= "     ";
    }
    else {
        $s .= sprintf "%4d ", $nsta;
    }
    if ($gap < 0) {
        $s .= "    ";
    }
    else {
        $s .= sprintf "%3d ", $gap;
    }
    if ($mind < 0) {
        $s .= "       ";
    }
    else {
        $s .= sprintf "%6.2f ", $mind;
    }
    if ($maxd < 0) {
        $s .= "       ";
    }
    else {
        $s .= sprintf "%6.2f ", $maxd;
    }
    $etype = set_etype(" ", $etp);
    $q = ($em =~ /manual/) ? "m" : "a";
    $s .= sprintf "%s   %-2s %-9s %-11s %-5s             %3d\n",
                   $q, $etype, $author, $id, $reporter, $sgap;
    return $s, $id, $author;
}

#
# magnitude line in ISF2
#
sub magnitudeline {
    my ($line, $orid) = @_;
    my ($mtype, $magnitude, $smag, $nsta, $auth, $s, $q);
    ($mtype, $magnitude, $smag, $nsta, $auth) =  split ' ', $line;
    $q = $mtype;
    $q = "" if ($mtype eq "Q");
    $q =~ s/\(/\_/;
    $q =~ s/\)//;
    $s = sprintf "%-5s ", $q;
    if ($magnitude < -8) {
        $s .= "     ";
    }
    else {
       $s .= sprintf "%4.1f ", $magnitude;
    }
    if ($smag < 0) {
        $s .= "    ";
    }
    else {
       $s .= sprintf "%3.1f ", $smag;
    }
    if ($nsta < 0) {
        $s .= "     ";
    }
    else {
       $s .= sprintf "%4d ", $nsta;
    }
    $s .= sprintf "%-9s %d\n", $auth, $orid;
    return $s;
}

#
# phase line in ISF2
#
sub phaseline {
    my ($k, $ht, $line, $snr, $type, $amp, $per, $ach, $ampmag) = @_;
    my ($pickid, $sta, $delta, $esaz, $phase, $date, $time, $ms, $tim);
    my ($tres, $azim, $ares, $slow, $sres, $tdef, $adef, $sdef);
    my ($pol, $onset, $net, $loc, $pch, $auth, $rep, $arid, $at, $tt);
    my ($s, $p, $hh, $mm, $ss);
    $s = "";
#
#   process phase line
#
    ($pickid, $sta, $delta, $esaz, $phase, $date, $tim, $ms, $tres,
     $azim, $ares, $slow, $sres, $tdef, $adef, $sdef, $pol, $onset,
     $net, $loc, $pch, $rep, $auth, $arid) = split ' ', $line;
    $date =~ s/\-/\//g; 
    $s .= sprintf "%-5s ", $sta;
    if ($delta < 0) {
        $s .= "       ";
    }
    else {
        $s .= sprintf "%6.2f ", $delta;
    }
    if ($esaz < 0) {
        $s .= "      ";
    }
    else {
        $s .= sprintf "%5.1f ", $esaz;
    }
    $phase = "" if ($phase eq "NULL" or $phase eq "Q");
    $s .= sprintf "%-8s ", $phase;
    if ($k) {
        $s .= "             ";
        $s .= "      ";
        $s .= "            ";
        $s .= "              ";
        $s .= "___ ";
        $s .= "      ";
    }
    else {
        if ($tim gt "00:00:00") {
            ($hh, $mm, $ss) = split '\:', $tim;
            ($p, undef) = split ' ', $ht;
            $hh += 24 if ($date gt $p);
            $s .= sprintf "%02d:%02d:%02d.%03d ", $hh, $mm, $ss, $ms / 1000;
        }
        else {
            $s .= "             ";
        }
        $tres = 99 if ($tres > 998);
        $tres = -99 if ($tres < -998);
        $ares = 99 if ($ares > 998);
        $ares = -99 if ($ares < -998);

        $s .= sprintf "%5.1f ", $tres;
        if ($azim < 0) {
            $s .= "            ";
        }
        else {
            $s .= sprintf "%5.1f %5.1f ", $azim, $ares;
        }
        if ($slow < 0) {
            $s .= "              ";
        }
        else {
            $s .= sprintf "%6.2f %6.2f ", $slow, $sres;
        }
        if ($tdef) {
            $s .= "T";
        }
        else {
            $s .= "_";
        }
        if ($adef) {
            $s .= "A";
        }
        else {
            $s .= "_";
        }
        if ($sdef) {
            $s .= "S ";
        }
        else {
            $s .= "_ ";
        }
        if ($snr < 0) {
            $s .= "      ";
        }
        else {
            $snr = 999.9 if ($snr > 1000.);
            $s .= sprintf "%5.1f ", $snr;
        }
    }
    if ($amp < 0) {
        $s .= "          ";
    }
    else {
        $s .= sprintf "%9.1f ", $amp;
    }
    if ($per < 0) {
        $s .= "      ";
    }
    else {
        $s .= sprintf "%5.2f ", $per;
    }
    if ($k) {
        $s .= "___ ";
    }
    else {
        $p = "_";
        if ($pol eq "positive" or $pol eq "c" or $pol eq "u" or $pol eq "+") {
            $p .= "c";
        }
        elsif ($pol eq "negative" or $pol eq "d" or $pol eq "-") {
            $p .= "d";
        }
        else {
            $p .= "_";
        }
        if ($onset eq "e" or $pol eq "i" or $pol eq "q") {
            $p .= $onset;
        }
        else {
            $p .= "_";
        }
        $s .= "$p ";
    }
    $type = "" if ($type eq "Q" or $type eq "NULL");
    $type =~ s/\(/\_/;
    $type =~ s/\)//;
    $s .= sprintf "%-5s ", $type;
    if ($ampmag < 0) {
        $s .= "     ";
    }
    else {
        $s .= sprintf "%4.1f ", $ampmag;
    }
    $s .= sprintf "%-11s %-5s %-8s %-2s %-5s %-5s %3s %3s _\n",
           $arid, "FDSN", $net, $loc, $auth, $rep, $pch, $ach;
    return $s;
}
