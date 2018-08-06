#
#
# MYSQLselect.pm: Convenience SQL routines for select, insert/update statements.
#
# Istvan Bondar, 2015/03/02
#
# Usage:
#      @rows = mysql_select($db_handle, $query);
#      mysql_do($db_handle, $query);
#      $nextid = mysql_newid($dbh, $tablename);
#
#
package MYSQLselect;

require Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw(mysql_select mysql_select2 mysql_do mysql_newid PublicObject);

use DBI;
use strict;

#
#
# mysql_select($dbh, $query);
#
#   executes sql select statement
#   returns each row (with newline) as a separate element of @rows
#   (list of strings)
#
#   example:
#      use DBI;
#      use Connect;
#      use MYSQLselect;
#      ...
#      $dbh = Connect::mysql_logon($user, $password);
#      $query = "select * from Event";
#      chomp(@rows = mysql_select($dbh, $query));
#      for (@rows) {
#          ($year, $n) = split;
#          ...
#      }
#      $dbh->Connect::mysql_logoff;
#
sub mysql_select {
    my ($dbh, $query) = @_;
    my ($csr, @row, $line);
    my @rows = ();
    $csr = $dbh->prepare($query) || die $dbh->errstr;
    $csr->execute || die "$query\n".$dbh->errstr;
    while (@row = $csr->fetchrow_array) {
        push @rows, "@row\n";
    }
    $csr->finish;
    return @rows;
}

#
#
# my_select2($dbh, $query);
#
#   executes sql select statement
#   returns each row (with newline) as a separate element of @rows
#   (list of lists)
#
#
sub mysql_select2 {
    my ($dbh, $query) = @_;
    my ($csr, @row, $line);
    my @rows = ();
    $csr = $dbh->prepare($query) || die $dbh->errstr;
    $csr->execute || die "$query\n".$dbh->errstr;
    while (@row = $csr->fetchrow_array) {
        push @rows, [ @row ];
    }
    $csr->finish;
    return @rows;
}

#
#
# mysql_do($dbh, $query);
#
#   executes sql update/insert statement
#   returns 0 on error, 1 otherwise
#
#
sub mysql_do {
    my ($dbh, $query) = @_;
    $dbh->do($query) or die "$query\n", $dbh->errstr;
    $dbh->commit()   or return 0;
    return 1;
}

#
#
# mysql_newid($dbh, $tablename)
#   $nextid = mysql_newid($dbh, "seiscomp3.Object");
#
#
sub mysql_newid {
    my ($dbh, $tablename) = @_;
    my ($query, $nextid, $csr);
    $query = "insert into $tablename (_timestamp) values (now())";    
    $csr = $dbh->prepare($query) || die $dbh->errstr;
    $csr->execute || die "$query\n".$dbh->errstr;
    $nextid = $dbh->{'mysql_insertid'};
    $csr->finish; 
    $dbh->commit();   
    return $nextid;
}
    
#
#
# PublicObject($dbh, $tablename)
#   PublicObject($dbh, "seiscomp3.PublicObject", $oid, $publicid);
#
#
sub PublicObject {
    my ($dbh, $tablename, $oid, $publicid) = @_;
    my ($query);
    $query = "replace into $tablename (_oid, publicID) 
               values ($oid, '$publicid')";
    mysql_do($dbh, $query) or die "$query failed!\n";    
}
    
1;

