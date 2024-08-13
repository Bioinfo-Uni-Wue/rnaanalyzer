# test-bioperl-modules.pl
#!/usr/bin/perl

use lib '/mnt/c/Users/ama55id/Nextcloud/RNA_analyzer/rnaanalyzer/cgi-bin';
use lib '/mnt/c/Users/ama55id/Nextcloud/RNA_analyzer/rnaanalyzer/bin/ViennaRNA-2.6.4/interfaces/Perl';

use strict;
use warnings;
use Bio::SeqIO;
use Bio::DB::GenBank;
use CGI;
use RNASERVER::TRANS2;
use RNASERVER::IRE;
use Bio::Tools::Genscan;
use Cwd;
use RNA; 



use strict;
use warnings;
use CGI qw(:standard);
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use lib '/mnt/c/Users/ama55id/Nextcloud/RNA_analyzer/rnaanalyzer/cgi-bin/RNASERVER';
use RNASERVER::TRANS2;
use RNASERVER::IRE;
use RNASERVER::RIBOSWITCH;
use Bio::SeqIO;
use Bio::Tools::Genscan;
use Bio::Tools::Run::RemoteBlast;
use LWP::UserAgent;

print header;
print start_html("RNA Analyzer");

# Add your script logic here

print end_html;



print "Bio PMs are installed and working\n";

