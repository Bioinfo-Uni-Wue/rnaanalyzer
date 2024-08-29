#!/usr/bin/perl -w
use strict;
my $query_string = $ENV{'QUERY_STRING'};
my $seq=$query_string;
my $INFILENAME = "/var/www/rnaanalyzer/htdocs/indexandom.html";

#binmode(STDOUT); 
print "Content-type: text/html\n\n";
$| = 1;
undef $/;
 
open (FILE, "<".$INFILENAME);
#binmode(PICTURE);
my @text=<FILE>;
my $text=join('',@text);
substr($text,index($text,"REMOVEMEPLEASE"),length("REMOVEMEPLEASE"))=$seq; 

print $text;




#my $buff;
#while (read(PICTURE, $buff, 8 * 2**10)) {
#    print STDOUT $buff;
#}

