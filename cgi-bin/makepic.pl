#!/usr/bin/perl
# written by liang, resposible for paintting a rna structure picture
# using ViennaRNA package
# Image::Magick modules are used to convert the ps pictures into jpg

use warnings;
use Image::Magick;
use CGI;
$co=new CGI;

#print $co->header(-type=>'text/html');
print $co->header(-type=>'image/jpeg');

my $TMPDIR="../tmp/";
my $RNAFoldLocation="/mnt/c/Users/ama55id/Nextcloud/RNA_analyzer/rnaanalyzer/bin/ViennaRNA-2.6.4/src/bin/RNAfold";

if ($co->param()) {
    my $seq=$co->param('seq');
    $image = Image::Magick->new;
    my $TMP=int(rand(65535));
    my $command='echo -e ">'.$TMPDIR.$TMP.'\n'.$seq.'\n" | '.$RNAFoldLocation.' > /dev/null';
    system($command);
    open(IMAGE, $TMPDIR.$TMP."_ss.ps");
    $image->Read(file=>\*IMAGE);
    binmode STDOUT;
    $image->write('jpeg:-');
    close(IMAGE);
    #system("rm ".$TMPDIR.$TMP."_ss.ps");
}
