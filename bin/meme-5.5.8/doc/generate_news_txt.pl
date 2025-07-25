#!/usr/bin/perl

my $state = 0;
while (<stdin>) {
  $state = 1 if ($_ =~ /\QSTART NEWS.TXT\E/);
  $state = 3 if ($_ =~ /\Q<!-- END NEWS.TXT -->\E/);
  if ($state == 1) { 
    $state = 2;
    print "#### MEME SUITE NEWS\n";
    next; 
  }
  if ($state == 2) {
    if ($_ =~ /MEME version (\d\.\d\.\d) \((\w+ (\d+, |)\d\d\d\d)\)/ ) {
      print "<p><hr></p>\n";
      print "* [$2] Release $1 MEME Suite\n";
    } else {
      print;
    }
  }
}
