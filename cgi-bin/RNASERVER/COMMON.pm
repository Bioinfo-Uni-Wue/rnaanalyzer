package Trans::common;

use strict;

sub correctgensq { #takes the sequence and checks it for incorrect bases and returns the corrected sequence
    my $seq=$_[0];
    #print "\n$seq\n++++++++";
    $seq=~s/[^a-zA-Z]//g;
    #print "\n$seq\n++++++++";
    $seq=lc($seq);
    #print "\n$seq\n++++++++";
    $seq=~s/[^agtuc]/n/g; #eliminate all other nts
    #print "\n$seq\n++++++++";
    $seq=~s/t/u/g;
    #print "\n$seq\n++++++++";
    return $seq;
}

return 1;
