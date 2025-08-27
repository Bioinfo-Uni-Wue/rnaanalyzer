#!/usr/bin/perl -w -I/storage/srv/bioapps/rnanalyzer/cgi-bin/RNASERVER

use CGI;
use RNASERVER::TRANS2;
use RNASERVER::IRE;
use Bio::Tools::Genscan;
use Cwd;

$rv=CGI::new();

$debug=1;

#Dir-Localisations
$TEMPPICSDIR='/var/www/rnaanalyzer/session/';
$TEMPDIR='/var/www/rnaanalyzer/tmp';
$IPADRESSFORPICTURE='wb2x01.biozentrum.uni-wuerzburg.de';
$GENSCANLIB='/var/www/rnaanalyzer/bin/genscanlinux';
$MAXFOLDINGLEN=1500;
$MAXFOLDINGLENUTR=1500;


print $rv->header();
#print '<body bgcolor="#000080" text="FFFFFF">';
print '<body text="#000000" bgcolor="#C0C0C0">';
print '<font face="monospace">';
print 'Please wait a moment for your results!';
$job=9;
@answergenscan=`/var/www/rnaanalyzer/bin/genscanlinux/genscan $GENSCANLIB/HumanIso.smat $TEMPDIR/$job.genscan >$TEMPDIR/$job.genscanout`;
@testlines=`grep "Intr" $TEMPDIR/$job.genscanout`;
print @testlines;
print "Done! @answergenscan  XXX \n";
&polyAsignal;

sub polyAsignal {
        print "--BeforePolyASignalInGenscan--" if ($debug); #have a look if a polyASignal has been found!!!
       # @polyAlines=`grep "PlyA" $TEMPDIR/80.genscanout`;
        $actualdir=cwd();
        chdir $TEMPDIR;
        print "--Now going to $TEMPDIR--";
        #@polyAlines=`grep "PlyA" $job.genscanout`;
        chdir $actualdir;
        print "--AfterPolyASignalInGenscan--@polyAlines--" if ($debug);
        if (@polyAlines==0) {
                print "<b>Poly-A Signal<sup>1</sup>:</b> none<br>";
        }
        else {
                print "<b>Poly-A Signal<sup>1</sup>:</b> start  -   end<br>";

                foreach $polyAanswer (@polyAlines) {
                        $polyAanswer=~/[0-9. ]+(PlyA) ([+-])[ ]+([0-9]+)/;
                        #print "<br> $polyAanswer !!<br>";
                        if ($1 eq 'PlyA') {
                                #Grep returned an PolyASignal line!
                                printf (" PlyA-Sgl:     %-6d - %6d<br>",$3,$3+5);
                                #print "<br>deb: $1 $2 $3 $4<br>";
                                @polyasignal=(@polyasignal,$3,($3+5));
                        }
                }
        }
}

