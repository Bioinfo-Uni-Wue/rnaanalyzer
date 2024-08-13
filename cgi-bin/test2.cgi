#!/usr/bin/perl -w -I/storage/srv/bioapps/rnanalyzer/cgi-bin/RNASERVER

use CGI;
use RNASERVER::TRANS2;
use RNASERVER::IRE;
use Bio::Tools::Genscan;
use Cwd;

$rv=CGI::new();

$TEMPDIR='/var/www/rnaanalyzer/tmp/';
$GENSCANLIB='/var/www/rnaanalyzer/bin/genscanlinux';
$job = "1";

print $rv->header();

my $cmdline = "/var/www/rnaanalyzer/bin/genscanlinux/genscan $GENSCANLIB/HumanIso.smat $TEMPDIR$job.genscan > $TEMPDIR$job.genscanout";
print $cmdline."\n";
system($cmdline);

system("cat $TEMPDIR$job.genscanout");

