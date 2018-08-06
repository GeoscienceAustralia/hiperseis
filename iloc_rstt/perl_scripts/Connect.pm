#
#DOC    File:         Connect_dbi.pm                                                     #
#DOC
#DOC    Package:      Connect_dbi;
#DOC
#DOC    This packages provides some standard routine for using DBI to
#DOC    log onto and from the ISC Oracle database.
#DOC
#DOC    subroutine:  oracle_logon      no arguments, returns the database
#DOC                                   handle or dies.
#DOC
#DOC    subroutine:  oracle_logoff     requires the database handle and
#DOC                                   returns true or false
#DOC
#DOC    subroutine:  mysql_logon       no arguments, returns the database
#DOC                                   handle or dies.
#DOC
#DOC    subroutine:  mysql_logoff      requires the database handle and
#DOC                                   returns true or false
#DOC
#DOC    subroutine:  psql_logon        no arguments, returns the database
#DOC                                   handle or dies.
#DOC
#DOC    subroutine:  psql_logoff       requires the database handle and
#DOC                                   returns true or false
#DOC
#

package Connect;
require 5.000;
require Exporter;
require DBI;
use Carp;

@ISA = qw(Exporter);
@EXPORT = qw(oracle_logon oracle_logoff mysql_logon mysql_logoff psql_logon psql_logoff);

#
#DOC
#DOC    subroutine:  oracle_logon
#DOC
#DOC                 provides a method of logging onto the ISC
#DOC                 database. The database handle is returned
#DOC                 or the program will be excited.
#DOC
#DOC                 The database driver and name are hardcoded
#DOC                 as for the forseeable time these are the
#DOC                 only options.
#DOC
#DOC                 Requires that a user has the ORACLE_USER
#DOC                 and ORACLE_PASS environment variables set.
#DOC
#DOC        Notes:   Possible, though currently not necessary
#DOC                 options would be to pass arguments for
#DOC                 different databases and none standard users.
#

sub oracle_logon {
    my $oracle_user = $ENV{ORACLE_USER};
    my $oracle_pass = $ENV{ORACLE_PASS};
    my $database_driver = "Oracle";
    my $database_name   = "ISCDB";
    my $dbh =  DBI->connect("DBI:$database_driver:$database_name",
                             $oracle_user,
                             $oracle_pass,
                             { AutoCommit=>0 }
                           )
            or croak "unable to connect to the database: $DBI::errstr";
    return $dbh;
}

#
#DOC
#DOC    subroutine:  oracle_logoff
#DOC
#DOC                 provides a method of logging off a session.
#DOC                 A value of 1 (succeeded) or 0 (failed) is
#DOC                 returned.
#DOC
#DOC                 Requires the database handle to be passed
#DOC                 as an argument.
#

sub oracle_logoff {
    my $dbh = $_[0];
    my $rc = $dbh->disconnect;
    if ( $rc == 1 ) {
        return 1;
    } else {
        warn "Unable to disconnect from the database: $DBI::errstr";
        return 0;
    }
}

#
#DOC
#DOC    subroutine:  mysql_logon
#DOC
#DOC                 provides a method of logging onto the ISC
#DOC                 database. The database handle is returned
#DOC                 or the program will be excited.
#DOC
#DOC                 The database driver and name are hardcoded
#DOC                 as for the forseeable time these are the
#DOC                 only options.
#DOC
#DOC                 Requires that a user has the MYSQL_USER
#DOC                 and MYSQL_PASS environment variables set.
#DOC
#DOC        Notes:   Possible, though currently not necessary
#DOC                 options would be to pass arguements for
#DOC                 different databases and none standard users.
#

sub mysql_logon {
    my ($database, $user, $pass, $host) = @_;
    my ($database_driver, $dbh);
#    $host = $ENV{MYSQL_HOST};
#    $database = $ENV{MYSQL_DB};
    if ( ! defined $user ) {
        $user = $ENV{MYSQL_USER};
        $pass = $ENV{MYSQL_PASS};
    }
    $database = 'seiscomp3' if ( ! defined $database );
    $host = 'localhost' if ( ! defined $host );
    $database_driver = "mysql";
    $dbh = DBI->connect("DBI:$database_driver:$database:$host",
                            $user,
                            $pass,
                            { AutoCommit=>0 }
                           )
            or croak "unable to connect to the database: $DBI::errstr";
    return $dbh;
}

#
#DOC
#DOC    subroutine:  mysql_logoff
#DOC
#DOC                 provides a method of logging off a session.
#DOC                 A value of 1 (failed) or 0 ( succeeded ) is
#DOC                 returned.
#DOC
#DOC                 Requires the database handle to be passed
#DOC                 as an arguement.
#

sub mysql_logoff {
    my $dbh = $_[0];
    my $rc = $dbh->disconnect;
    if ( $rc == 1 ) {
        return 1;
    } else {
        warn "Unable to disconnect from the database: $DBI::errstr";
        return 0;
    }
}

#
#DOC
#DOC    subroutine:  psql_logon
#DOC
#DOC                 provides a method of logging onto the ISC
#DOC                 database. The database handle is returned
#DOC                 or the program will be excited.
#DOC
#DOC                 The database driver and name are hardcoded
#DOC                 as for the forseeable time these are the
#DOC                 only options.
#DOC
#DOC                 Requires that a user has the PSQL_USER
#DOC                 and PSQL_PASS environment variables set.
#DOC
#DOC        Notes:   Possible, though currently not necessary
#DOC                 options would be to pass arguements for
#DOC                 different databases and none standard users.
#

sub psql_logon {
    my ($user,$pass) = @_;
    if ( ! defined $user ) {
        $user = $ENV{PGUSER};
        $pass = $ENV{PGPASSWORD};
    }
    my $db = $ENV{PGDATABASE};
    $db = 'isc' if ( ! defined $db );
    my $dbh = DBI->connect("DBI:Pg:dbname=$db;","$user","$pass",{ AutoCommit=>0 }) or
              croak "unable to connect to the database: " . $dbh->errstr;
    $dbh->{RaiseError} = 0;
    $dbh->{PrintError} = 0;
    $dbh->{AutoCommit} = 0;
    if ( $db eq 'isc' ) {
        if ( $user eq 'capture' ) {
            $dbh->do("SET search_path TO isc") or die "unable to set search_path: " . $dbh->errstr;
        } elsif ( $user eq 'web' ) {
            $dbh->do("SET search_path TO isc") or die "unable to set search_path: " . $dbh->errstr;
        }
    }
    return $dbh;
}

#
#DOC
#DOC    subroutine:  psql_logoff
#DOC
#DOC                 provides a method of logging off a session.
#DOC                 A value of 1 (failed) or 0 ( succeeded ) is
#DOC                 returned.
#DOC
#DOC                 Requires the database handle to be passed
#DOC                 as an arguement.
#

sub psql_logoff {
    my $dbh = $_[0];
    my $rc = $dbh->disconnect;
    if ( $rc == 1 ) {
        return 1;
    } else {
        warn "Unable to disconnect from the database: $DBI::errstr";
        return 0;
    }
}


#
#DOC
#DOC    subroutine:  check_edit2_connections
#DOC
#DOC                 provides a method of checking users log in session
#DOC                 A value of 1 - only one connection for use -> proceede returned.
#DOC                 A value of 0 - more than one connection fro the sameuser connected .
#DOC
#DOC                 Requires the number of loops
#DOC                 as an arguement.
#

sub check_edit2_connections{

    my ($count) = @_;
    my $dbh = psql_logon;
    my @pgusers;
    my $count_users = 0;
    my $next = 0;
    my $sql = "SELECT usename FROM pg_stat_activity";
    my $csr = $dbh->prepare($sql) or die $dbh->errstr;
    $csr->execute() or die "$sql\n".$dbh->errstr;
    while (my ($user) = $csr->fetchrow_array){
        push @pgusers, $user;
    }
    $csr->finish;
    foreach my $user (@pgusers){
         $count_users++ if ($user eq "edit2" and $ENV{PGUSER} eq "edit2" );
    }
    $next = 1 if ($count_users < 2 );
    sleep 1;
    #print "\n$count: @pgusers - $count_users\n";
    $dbh->Connect::psql_logoff or warn $dbh->errstr;
    print STDERR "\t$count $next $count_users $ENV{PGUSER}\n";
    if ($count and !$next) {
        print STDERR "\tSomeone is using edit2 schema ... $count $next $count_users $ENV{PGUSER}";
        print STDERR "\b" x 164;
        print STDERR "    " . "\t\t\t\t\t\tWaiting";
        sleep 1;
        print STDERR "\b" x 164;
        print STDERR "    " . "\t\t\t\t\t\t        ";

    }
    $count++;
    return($count,$next);
}


1;
#
##  End of file Connect.pm
#

