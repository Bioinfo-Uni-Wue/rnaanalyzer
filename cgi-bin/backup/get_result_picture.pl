#!/usr/bin/perl

use strict;
my $query_string = $ENV{'QUERY_STRING'};
$query_string=~/([0-9]+_[0-9abcdef]{32})/;
my $id=$1;
my $INFILENAME = "/var/www/rnaanalyzer/tmp/$id.jpg";
binmode(STDOUT);
 
print "Content-type: image/jpeg\n\n";
$| = 1;
undef $/;
 
open (PICTURE, "<".$INFILENAME);
binmode(PICTURE);
 
my $buff;
while (read(PICTURE, $buff, 8 * 2**10)) {
    print STDOUT $buff;
}

