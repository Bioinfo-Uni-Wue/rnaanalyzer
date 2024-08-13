#!/usr/bin/perl
# written by liang, resposible for paintting a rna structure picture
# using ViennaRNA package
# Image::Magick modules are used to convert the ps pictures into jpg


use Image::Magick;
use CGI;
$co=new CGI;

#print $co->header(-type=>'text/html');
print $co->header(-type=>'image/jpeg');

my $TMPDIR="../tmp/";
my $RNAFoldLocation="/var/www/rnaanalyzer/bin/ViennaRNA-2.4.18/src/binRNAfold";

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
    system("rm ".$TMPDIR.$TMP."_ss.ps");
}
