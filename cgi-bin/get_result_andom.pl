#!/usr/bin/perl -w
use strict;
my $query_string = $ENV{'QUERY_STRING'};
$query_string=~/([0-9]+_[0-9abcdef]{32})/;
my $id=$1;
my $INFILENAME = "/mnt/c/Users/ama55id/Nextcloud/RNA_analyzer/rnaanalyzer/tmp/indexandom$id.html";
binmode(STDOUT);
 
print "Content-type: text/html\n\n";
$| = 1;
undef $/;
 
open (PICTURE, "<".$INFILENAME);
binmode(PICTURE);
 
my $buff;
while (read(PICTURE, $buff, 8 * 2**10)) {
    print STDOUT $buff;
}

