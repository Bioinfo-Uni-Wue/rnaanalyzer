<<<<<<< HEAD
#!/usr/bin/perl

use lib ".";
use lib "./RNASERVER/";
=======

>>>>>>> origin/master
use CGI;
use RNASERVER::TRANS2;
use RNASERVER::IRE;
use Bio::Tools::Genscan;
use Cwd;
<<<<<<< HEAD
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use RNASERVER::JobUtil qw(get_next_job_id);
use JSON;
use File::Slurp;
use File::Path qw(make_path);
use Cwd qw(abs_path);
use CGI qw(:standard escapeHTML);
use Bio::Seq;

=======

$rv=CGI::new();
>>>>>>> origin/master

$debug=0;

#Dir-Localisations
<<<<<<< HEAD
$TRNASCANFOLDER=abs_path('../bin/tRNAscan-SE/bin');
$VIENNARNAFOLDDIR=abs_path('../bin/ViennaRNA-2.7.0/bin'); #pointing to the RNAfold dir
$CMSCAN=abs_path('../bin/tRNAscan-SE/infernal-1.1.5/src'); #CMSCAN directory alrady present in transcan
$RFAM=abs_path('../databases/rfam/Rfam.cm'); #RFAM database
$CPC=abs_path('../bin/cpc2/CPC2_standalone-1.0.1/bin'); #Coding potential calculator 2 location
$HMMER=abs_path('../bin/hmmer-3.4/bin'); #hmmer location 
$MIRBASE=abs_path('../databases/mirbase'); #MirBASE database
$MIRANDA=abs_path('../bin/miranda/bin'); #miRanda 
$AUGUSTUS=abs_path('../bin/Augustus/bin/augustus'); #augugtus 
$RBSFINDER=abs_path('../bin/rbs-finder'); #rbs finder
$FIMO=abs_path('../bin/meme-5.5.8/bin'); #fimo rbp scan
$RBPDB=abs_path('../databases/rbpdb/rbpdb.meme');
$RIBOSWDB=abs_path('../databases/riboswitches/riboswitch.cm'); #riboswitch location 
$MAXFOLDINGLEN=5000;
$MAXFOLDINGLENUTR=5000;
$MAXFORNALENGTH=5000;

# aprsing cgi
my $cgi = CGI->new;


$job_id = $ARGV[0] // $cgi->param("job_id");  # passed from batch script
$job = $job_id;


$TEMPDIR = abs_path("../tmp/jobs/job_$job");


unless (-d $TEMPDIR) {
    make_path($TEMPDIR) or die "Can't create job dir: $TEMPDIR";
}
print "<p>Job submitted: $job</p>";


open my $log, ">>", "../tmp/backend_log.txt";
print $log scalar localtime() . " Received job: $job\n";
close $log;

# chcking if backend is alive
# print $cgi->header('text/plain');
# print "Backend is alive. Job ID: $job_id\n";

my $params_json = read_file("$TEMPDIR/params.json");
my $params = decode_json($params_json);

# Now extract them like form params:
my $do_rnamotif         = $params->{RNAmotif};
my $do_augustus         = $params->{run_coding};
my $do_mirna            = $params->{mirna};
my $do_trna             = $params->{trna};
my $do_IRE              = $params->{IRE};
my $do_TRANS            = $params->{TRANS};
my $species             = $params->{species};
my $do_mirnatarget      = $params->{mirna_target};
my $do_rbp              = $params->{RBP};
my $do_ribo             = $params->{ribo};
# my $dnarna              = $params->{dnarna};
my $SEQNAMECHECKED      = $params->{sequence_name};
my $SEQUENCECHECKED     = $params->{sequence_clean};
my $SEQUENCELENGTH      = $params->{sequence_length};


my $html_file = "$TEMPDIR/result.html";

open my $out, ">", $html_file or die "Can't write result page: $!";
select $out;  # Redirect STDOUT to file

# We print into one html which can be accessed again or bookmarked 

# setting up the html headers and links 
print "<!DOCTYPE html>";
print "<html lang='en'>";
print "<html>";
print "<head>";
print "  <title>Batch Results</title>";
print "  <link rel='stylesheet' href='/css/results.css'>";
print "</head><body>";
print "<header>";
print "    <a href='http://localhost'>";    # change after putting on server
print "      <img src='http://localhost/images/logo.png' alt='RNA Analyzer Logo' class='logo' />";
print "    </a>";
print "    <div class='header-text'>";
print "      <h1>RNA Analyzer<sup>3</sup></h1>";
print "      <p>Webserver for RNA Sequence Overview</p>";
print "    </div>";
print "    <div class='header-links'>";
print "      <a href='http://localhost/about.html' target='_blank'>About</a> |";
print "      <a href='http://localhost/contact.html' target='_blank'>Contact</a> |";
print "      <a href='https://www.biozentrum.uni-wuerzburg.de/bioinfo' target='_blank'>Dandekar Lab</a>";
print "    </div>";
print "  </header>";
print "";
print "<main>";
print "<div class='container'>";

   
# Start processing with the single sequence input.
&startproggi;

sub startproggi {
	#For the picture creation do not remove
	open (SEQPIC,">$TEMPDIR/$job.seq"); #don't know if this works!!	
	print SEQPIC ">$SEQNAMECHECKED\n$SEQUENCECHECKED\n"; #shall create a fasta-format sequence file !!!
	close SEQPIC;

    # creating dna sequence file for dna dependent programs 
    (my $dna_sequence = $SEQUENCECHECKED) =~ tr/uU/tT/;
    open(my $fh, '>', "$TEMPDIR/$job.dna") or die "Cannot write input DNA file: $!";
    print $fh ">$job\n$dna_sequence\n";
    close($fh);


    ####### Initializing some variables for the colored output ########
    @exons=();@transsplicing=();@ire=();@smsite=();@aurichregion=();
    @stemggpairs=();@polyasignal=();@utr=();@cds=();
    our @polyatail;
    our @rna_motif;
    our @mirna_loc;
    our @trna_loc;

    print "<div class='info-tip'>";
    print "<b>JOB ID:</b> $job";
    print "<br>";
    print "<b>Sequence Name:</b> ". CGI::escapeHTML($SEQNAMECHECKED);
    print "</div>";
    
    
    if ($SEQUENCELENGTH >= 2500){
            print "<div class='info-warning'>";
            print "<i>\ti.    Structure visualization can take a few seconds to load due to sequence length.</i>\n";
            print "</div>";
        }
    
    print "<br>";
    
        # running analysis
    &analysis;
}

##calls all the new subrotines
sub analysis {
    chdir $TEMPDIR;

    # print "<pre>";

    print "<div class='box'>";
    print "<div class='box-header' onclick='toggleBoxContent(this)'>Structural information</div>";
    print "<div class='box-content-structure'>";
    if (length $SEQUENCECHECKED <= $MAXFOLDINGLEN) {
        &createfolding;
        &checkstems;
        &stemggpairs;
    } else {
        print "<br><b>Length:</b>\t$SEQUENCELENGTH";
        print "<div class='info-error'>";
        print "Some information is only available up to $MAXFOLDINGLEN nt\n" if ($SEQUENCELENGTH > $MAXFOLDINGLEN);
        print "</div>";
    }
    print "</div>";
    print "</div>";

    if ($do_TRANS) {
    print "<div class='box'>";
    &TRANS;
    print "</div>";
    }

    if ($do_IRE) {
    print "<div class='box'>";	
	&IRE;
    print "</div>";
    }

    if ($do_ribo) {
    print "<div class='box'>";	
	&riboswitch;
    print "</div>";
    }

    print "<div class='box'>";
    &ARE;
    print "</div>";

    print "<div class='box'>";
    &smsite;
    print "</div>";

    if ($do_rnamotif) {
    print "<div class='box'>";
    &RNAMOTIF;
    print "</div>";
    }

    if ($do_trna) {
    print "<div class='box'>";
    &tRNA;
    print "</div>";
    }

    if ($do_rbp) {
    print "<div class='box'>"; 
    &rbp;
    print "</div>";
    }

    if ($do_mirna) {
    print "<div class='box'>";
    &microRNA;
    print "</div>";
    }

    # --- Capture transcript models ---
    print "<div class='box'>";
    my %transcripts;
    if ($do_augustus) {
        %transcripts = %{ AUGUSTUS() };   # Returns hashref
    } else {
        %transcripts = %{ CPC2() };       # Also returns hashref
        $cpc = "TRUE"; 
    }
    print "</div>";

    print "<div class='box'>";
    foreach my $tid (keys %transcripts) {
        my $model = $transcripts{$tid};

        my ($utr5_ref, $utr3_ref, $utrprintout_ref, $utr_coords_ref, $polyasignal_ref, $polyatail_ref) = predict_utrs(
            
            
            seq        => $SEQUENCECHECKED,
            cds_start  => $model->{cds_start},
            cds_end    => $model->{cds_end},
            strand     => $model->{strand},
            tss        => $model->{tss},
            tts        => $model->{tts},
            source     => $model->{source},
        );

        if (!@$utr3_ref) {
        print "<i>No 3' UTR predicted. Trying polyA inference...</i><br>";
        my ($polyautr_ref, $signals_ref, $tails_ref) = refineUTRwithPolyA(
            $SEQUENCECHECKED,
            $model->{cds_end},
            $model->{strand},
            length($SEQUENCECHECKED)
        );

        # Instead of storing references, store copies:
        $model->{polyAcoords}     = [@$polyasignal_ref] if $polyasignal_ref && @$polyasignal_ref;
        $model->{polyAtailcoords} = [@$polyatail_ref]   if $polyatail_ref   && @$polyatail_ref;

        $model->{utr3} = $polyautr_ref if $polyautr_ref && @$polyautr_ref;
        $utr3_ref = $polyautr_ref      if $polyautr_ref && @$polyautr_ref;
        }

        # Optionally store results back into the model
        # Always store UTRs
        $model->{utr5}         = $utr5_ref;
        $model->{utr3}         = $utr3_ref;
        $model->{utr_coords}   = $utr_coords_ref;

        #  store polyA info if found by predict_utrs
        # Instead of storing references, store copies:
        $model->{polyAcoords}     = [@$polyasignal_ref] if $polyasignal_ref && @$polyasignal_ref;
        $model->{polyAtailcoords} = [@$polyatail_ref]   if $polyatail_ref   && @$polyatail_ref;   
    }
    print "</div>";

    normalize_transcript_features(
        \%transcripts,
        \@exons,
        \@utr,
        \@polyasignal,
        \@polyatail,
        \@mirnatarget,
        \@structured_regions
    );

    print "<div class='box'>";
    &predprotein;
    print "</div>";

    if ($do_mirnatarget) {
    print "<div class='box'>";
    miRNAtarget(\%transcripts);
    print "</div>";
    } 
    # Optional: enable if scan_rbs supports multi-transcript
    # scan_rbs(\%transcripts);

    print "<div class='box'>";
    &csfce;
    print "</div>";

    print "<div class='box'>";
    &createfoldingpicture;
    print "</div>";

    print "<div class='table-container'>";
    &location_table;
    print "</div>";

    # &drawcoloredsequence;

    # print "</pre>";
}

sub TRANS{
	#this checkes the ciona-consensus !		
			@transcionareturnvalues=RNASERVER::TRANS2::ciona($SEQUENCECHECKED); ## Ist das der CIONA ??? Auf jeden Fall aber SCHISOTSOMA
			print "<div class='box-header' onclick='toggleBoxContent(this)'>Trans-Splicing Analysis:</div>";
			#print "<b> <big>Putative trans-splicing Schistosoma-consensus</b></big> search:";
            print "<div class='box-content'>";
			if (@transcionareturnvalues==1) {
				print "<div class='info-warning'>";
				print "<b>Schistosoma:</b>\tNo Trans-splicing element found.<br>";
                print "</div>";
=======
$TEMPPICSDIROA='/var/www/rnaanalyzer/session/';
$TEMPDIR='/var/www/rnaanalyzer/tmp';
#$IPADRESSFORPICTURE='wb2x01.biozentrum.uni-wuerzburg.de';
$GENSCANLIB='/var/www/rnaanalyzer/bin/genscanlinux';
$GENSCANFILE='/var/www/rnaanalyzer/bin/genscanlinux';
$TRNASCANFOLDER='/var/www/rnaanalyzer/bin/tRNAscan-SE/bin';
$VIENNAFOLDDIR='/var/www/rnaanalyzer/bin/ViennaRNA-2.4.18/src/Utils'; #pointing to the the Program Fold of the Vienna package
$VIENNARNAFOLDDIR='/var/www/rnaanalyzer/bin/ViennaRNA-2.4.18/src/bin'; #pointing to the RNAfold dir
$MAXFOLDINGLEN=1500;
$MAXFOLDINGLENUTR=1500;


print $rv->header();
#print '<body bgcolor="#000080" text="FFFFFF">';
print '<body text="#000000" bgcolor="#C0C0C0">';
print '<font face="monospace">';
print 'Please wait a moment for your results!';

$errorsyetprintout=0; #indicates that an error message has not yet been print out! Will be set to 1 if that has happened!

###### Next lines will provide a way to integrate the use of more than 1 sequence
$BUNCHOFFASTA=$rv->param("FASTA");
$BUNCHOFFASTA=~/^(ON)/;
$BUNCHOFFASTACHECKED=$1;

if ($BUNCHOFFASTA eq 'ON') {
	$DOFASTA=1;
	$xxx='xxx';
	$xxx=~/(xxx)/;

	$SEQUENCESTART=$rv->param("SEQUENCE");
	$SEQUENCESTART=~/(.*)/s;
	$SEQUENCESTARTCHECKED=$1;
	#print "<br><b>$SEQUENCESTARTCHECKED</b><br>";
	$numberoffastasequences=0;
	@FASTASEQ=();
	while ($SEQUENCESTARTCHECKED=~/>/g) {
		$numberoffastasequences++;
	}
	if ($numberoffastasequences<1){
		print '<br><br><b><big>Not ">"  detected! Your Seq is not correct!!<b></big><br>' ;
		print "Go back and correct your sequence or uncheck the SET OF FASTA - box !!!<br>";
		die;
	}
	@FASTASEQ=split (/>/,$SEQUENCESTARTCHECKED);
		#}
		#else {
			#$FASTASEQ[0]=">".$SEQUENCESTARTCHECKED;
			#}
	###debugging ####
	#foreach $fc (@FASTASEQ) {
		#	print "DEBUGGING: $fc<br>";
		#}
	
	#### Okay till now we have checked if there is more than one fasta and we have put each entry into an array
	for ($seqcountfs=1;$seqcountfs<@FASTASEQ;$seqcountfs++){
		$SEQUENCE=$FASTASEQ[$seqcountfs];
		$SEQUENCE='>'.$SEQUENCE;
		
		my $p_ire_switch=$rv->param("IRE");
		my $p_trans_switch=$rv->param("TRANS");
		my $p_fasta_switch=$rv->param("FASTA");
		my $p_origin_switch=$rv->param("ORIGIN");
		my $SEQUENCEX=$SEQUENCE;
		substr($SEQUENCEX,index($SEQUENCEX,"\n"),0)="%0A%0D";

		print '<br>';
		print '<iframe width="920" height="1800" frameborder="0" src="';
		print 'single_seq_rna_analysis.cgi?IRE='.$p_ire_switch;
		print '&ORIGIN='.$p_origin_switch;
		print '&TRANS='.$p_trans_switch.'&FASTA='.$p_fasta_switch.'&SEQUENCE=';
		print "$SEQUENCEX";
		print '"></iframe><hr><br>';


# comment by liang
#		$job=&jobnumber; 
#                #creates a unique number for job! herewith can manage the deal with the picture!!!!
#
#		&startproggi;
	}
}
else {
	    $DOFASTA=0;
	    $SEQUENCESTART=$rv->param("SEQUENCE");
	    $SEQUENCESTART=~/(.*)/s;
	    $SEQUENCESTARTCHECKED=$1;
	    $SEQUENCE=$SEQUENCESTARTCHECKED;
	    $job=&jobnumber; 
             #creates a unique number for job! herewith can manage the deal with the picture!!!!!
	    &startproggi;

}

#&startproggi;

sub startproggi {
	&readandcheckinput;
	####### Initializing some variables for the colored output########
	@exons=();@transsplicing=();@ire=();@smsite=();@aurichregion=();@stemggpairs=();@polyasignal=();@utr=();@promotor=();

	#print $rv->h1("RESULTS:") if  ($passsequence==1 && $passconsensus==1);
	print "<br><br><hr>";
	print $rv->h2("Here are the results for: $SEQNAMECHECKED");
	print $rv->h1("General Information");

#by liang
	&additionalinformation;
	
	print '<font face="monospace">';
	print $rv->h1("Special RNA structure Information");

	&mainrun;

	#&drawcoloredsequence;

	#$DOOPPOSITE=1;
	if ($DOOPPOSITE==1) {
		&oppositestrand;
		print $rv->h3("Now searching the OPPOSITE strand") if  ($passsequence==1 && $passconsensus==1);
		&mainrun;
		&oppositestrand;
	}

&drawcoloredsequence;
print "<hr>References:<br><sup>1</sup>Burge, C. and Karlin, S. (1997) Prediction of complete gene structures in human genomic DNA. J. Mol. Biol. 268, 78-94 ";
print "<br><sup>2</sup>Lowe, T.M. and Eddy, S.R. (1997) tRNAscan-SE: A program for improved detection of transfer RNA genes in genomic sequence, Nucl. Acids Res., 25, 955-964.<br>";

#print "XXXXXXXXXXXXXXXX";
}
#################################################################
#################################################################
#################################################################

sub mainrun {
	print "<pre>";
	if  ($passsequence==1 && $passconsensus==1) {
		if ($DOTRANS==1) {
			#this checkes the ciona-consensus !		
			@transcionareturnvalues=RNASERVER::TRANS2::ciona($SEQUENCECHECKED); ## Ist das der CIONA ??? Auf jeden Fall aber SCHISOTSOMA
			print "<b>Trans-Splicing:</b><br>";
			#print "<b> <big>Putative trans-splicing Schistosoma-consensus</b></big> search:";
			if (@transcionareturnvalues==1) {
				#print "No hit detected<br>";
				print " Schistosoma: none<br>";
>>>>>>> origin/master
			}
			else {
				$hits=pop @transcionareturnvalues;
				for ($count=0;$count<$hits;$count++){
<<<<<<< HEAD
					# print " Schistosoma: <br>";
					# print "  Position:   $transcionareturnvalues[$count*6+0] with 35 bp ustream region<br>";
					# $transcionareturnvalues[$count*6+4]=uc ($transcionareturnvalues[$count*6+4]);
					# $transcionareturnvalues[$count*6+2]=uc ($transcionareturnvalues[$count*6+2]);
					# print "  Stem1:      $transcionareturnvalues[$count*6+4]<br>";
					# print "  Structure:  $transcionareturnvalues[$count*6+5]<br>";
					# print "  Energy:     $transcionareturnvalues[$count*6+3]<br>";
					# print "  Sm-Site:    $transcionareturnvalues[$count*6+2] at pos: $transcionareturnvalues[$count*6+1]<br>";

                    print "<h3>Schistosoma:</h3>\n";
                    print "<table class='table-result'>\n";
                    print "<tr><td>Position</td><td>$transcionareturnvalues[$count*6+0] with 35 bp upstream region</td></tr>\n";

                    $transcionareturnvalues[$count*6+4] = uc($transcionareturnvalues[$count*6+4]);
                    $transcionareturnvalues[$count*6+2] = uc($transcionareturnvalues[$count*6+2]);

                    print "<tr><td>Stem1</td><td>$transcionareturnvalues[$count*6+4]</td></tr>\n";
                    print "<tr><td>Structure</td><td>$transcionareturnvalues[$count*6+5]</td></tr>\n";
                    print "<tr><td>Energy</td><td>$transcionareturnvalues[$count*6+3]</td></tr>\n";
                    print "<tr><td>Sm-Site</td><td>$transcionareturnvalues[$count*6+2] at pos: $transcionareturnvalues[$count*6+1]</td></tr>\n";

                    print "</table>\n";

					@transsplicing=(@transsplicing,$transcionareturnvalues[$count*6+0]-35,($transcionareturnvalues[$count*6+1]+length $transcionareturnvalues[$count*6+2])-1); #for the formatted output of sequence # counting 35 bp upstream why?
=======
					print " Schistosoma: <br>";
					print "  Position:   $transcionareturnvalues[$count*6+0]<br>";
					#print "<br><br>Possible hit detected at nt: $transcionareturnvalues[$count*6+0] (pointing to the ggua of stem1)<br>";
					$transcionareturnvalues[$count*6+4]=uc ($transcionareturnvalues[$count*6+4]);
					$transcionareturnvalues[$count*6+2]=uc ($transcionareturnvalues[$count*6+2]);
					print "  Stem1:      $transcionareturnvalues[$count*6+4]<br>";
					print "  Structure:  $transcionareturnvalues[$count*6+5]<br>";
					print "  Energy:     $transcionareturnvalues[$count*6+3]<br>";
					print "  Sm-Site:    $transcionareturnvalues[$count*6+2] at pos: $transcionareturnvalues[$count*6+1]<br>";
					#print '<table border=1><tr><td align="center"><font face="monospace">Stem1:</td><td align="center"><font face="monospace">'."$transcionareturnvalues[$count*6+4]".'</td></tr>';
					#print '<tr><td align="center"><font face="monospace">Structure:</td><td align="center"><font face="monospace">'."$transcionareturnvalues[$count*6+5]".'</td></tr>';
					#print '<tr><td align="center"><font face="monospace">Energy:</td><td align="center"><font face="monospace">'."$transcionareturnvalues[$count*6+3]".'</td></tr>';
					#print '<tr><td align="center"><font face="monospace">SM-Site:</td><td align="center"><font face="monospace">'."$transcionareturnvalues[$count*6+2] at position: $transcionareturnvalues[$count*6+1]".'</td></tr>';
					#print '</table><br>';
					@transsplicing=(@transsplicing,$transcionareturnvalues[$count*6+0]-35,($transcionareturnvalues[$count*6+1]+length $transcionareturnvalues[$count*6+2])-1); #for the formatted output of sequence
>>>>>>> origin/master
				}
			}
			#this checkes the c. elegans consensus !!
			@transcelegansvalues=RNASERVER::TRANS2::celegans($SEQUENCECHECKED);	
<<<<<<< HEAD
			if (@transcelegansvalues==1){
				#print "No hit detected<br>";
                print "<div class='info-warning'>";
				print "<b>C. elegans:</b>\tNo Trans-splicing element found.";
                print "</div>";
=======
			#print "<b> <big>Putative trans-splicing C.elegans-consensus</b></big> search:";
			if (@transcelegansvalues==1){
				#print "No hit detected<br>";
				print " C. elegans:  none<br>";
>>>>>>> origin/master
			}
			else {
				$hits=pop @transcelegansvalues;
				for ($count=0;$count<$hits;$count++){
					$transcelegansvalues[$count*10+5]=uc($transcelegansvalues[$count*10+5]);
<<<<<<< HEAD
					# print " C.elegans:<br>  Position:   $transcelegansvalues[$count*10+0] with 21 bp ustream region ; pointing to ggua of stem1<br>";
					
					# print "  Stem1:      $transcelegansvalues[$count*10+1]<br>";
					# print "  Structure:  $transcelegansvalues[$count*10+2]<br>";
					# print "  Stem2:      $transcelegansvalues[$count*10+3]<br>";
					# print "  Structure:  $transcelegansvalues[$count*10+4]<br>";
					# print "  Sm-Site:    $transcelegansvalues[$count*10+5]<br>";
					# print "  Stem3:      $transcelegansvalues[$count*10+6]<br>";
					# print "  Structure:  $transcelegansvalues[$count*10+7]<br>";
					# print "  Leader:     $transcelegansvalues[$count*10+8]<br>" if ($transcelegansvalues[$count*10+8] !=0);
					# print "  Leader:     none<br>" if ($transcelegansvalues[$count*10+8] == 0);
                    print "\n<h3>C. elegans:</h3>\n";

                    print "<table class='table-result'>\n";
                    print "<tr><td><b>Position</b></td><td>$transcelegansvalues[$count*10+0] with 21 bp upstream region ; pointing to ggua of stem1</td></tr>\n";
                    print "<tr><td><b>Stem1</b></td><td>$transcelegansvalues[$count*10+1]</td></tr>\n";
                    print "<tr><td><b>Structure (1)</b></td><td>$transcelegansvalues[$count*10+2]</td></tr>\n";
                    print "<tr><td><b>Stem2</b></td><td>$transcelegansvalues[$count*10+3]</td></tr>\n";
                    print "<tr><td><b>Structure (2)</b></td><td>$transcelegansvalues[$count*10+4]</td></tr>\n";
                    print "<tr><td><b>Sm-Site</b></td><td>$transcelegansvalues[$count*10+5]</td></tr>\n";
                    print "<tr><td><b>Stem3</b></td><td>$transcelegansvalues[$count*10+6]</td></tr>\n";
                    print "<tr><td><b>Structure (3)</b></td><td>$transcelegansvalues[$count*10+7]</td></tr>\n";

                    if ($transcelegansvalues[$count*10+8] != 0) {
                        print "<tr><td>Leader</td><td>$transcelegansvalues[$count*10+8]</td></tr>\n";
                    } else {
                        print "<tr><td>Leader</td><td>none</td></tr>\n";
                    }
                    
                    print "</table>";

					@transsplicing=(@transsplicing,$transcelegansvalues[$count*10+0]-21,($transcelegansvalues[$count*10+9]+(length $transcelegansvalues[$count*10+5])+(length $transcelegansvalues[$count*10+6]))); # counting 21 upsteam but why?
				}
			}
    print "</div>";
    # print "DEBUG: @transsplicing";
}

sub IRE {
    @irereturnvalues = RNASERVER::IRE::suboptimalfindire($SEQUENCECHECKED);
    $irelineprintout = 0;
    print "<div class='box-header' onclick='toggleBoxContent(this)'>Iron-resp Element(s):</div>";
    print "<div class='box-content'>";
    
    if (@irereturnvalues > 1) {
        # Group results by position
        my %positions;
        for (my $count = 0; $count <= @irereturnvalues - 1; $count += 7) {
            my $pos = $irereturnvalues[$count + 0];
            push @{$positions{$pos}}, {
                position => $irereturnvalues[$count + 0],
                quality => $irereturnvalues[$count + 1],
                sequence => $irereturnvalues[$count + 2],
                structure => $irereturnvalues[$count + 3],
                energy => $irereturnvalues[$count + 4]
            };
        }
        
        # Display each position with its structures
        foreach my $pos (sort {$a <=> $b} keys %positions) {
            my @structures = @{$positions{$pos}};
            
            print "<table class='table-result'>\n";
            print "<tr><th>Position</th><td>$pos</td></tr>\n";
            print "<tr><th>Sequence</th><td>$structures[0]->{sequence}</td></tr>\n";
            
            @ire = (@ire, $pos - 16, $pos + 22);
            
            # Display all structures for this position
            for (my $i = 0; $i < @structures; $i++) {
                my $struct_num = $i + 1;
                my $struct = $structures[$i];
                
                print "<tr><th>Structure $struct_num</th><td>$struct->{structure}</td></tr>\n";
                printf "<tr><th>Energy $struct_num</th><td>%.2f kcal/mol</td></tr>\n", $struct->{energy};
                
                my $quality_text = ($struct->{quality} == 1) ? "good" : "bad";
                print "<tr><th>Quality $struct_num</th><td>$quality_text</td></tr>\n";
                
                # Add separator line between structures (except for the last one)
                if ($i < @structures - 1) {
                    print "<tr><td colspan='2' style='border-bottom: 1px dashed #ccc; padding: 5px;'></td></tr>\n";
                }
            }
            
            print "</table>";
        }
        
        $irelineprintout = 1;
        
    } else {
        print "<div class='info-warning'>";
        print "<b>Iron-resp Ele.:</b> None";
        print "</div>";
    }
    
    print "</div>";
}

sub csfce {
    
    my $count1=0; #We will mark $count1 as my so that we won't have any problems later
    #$seq='llllllllllllllllluugculllauuuacuglcculllaugcguuccucgucclllllllllllllllll';
    my @seq=split ('',$SEQUENCECHECKED);
    #my @seq=split ('',$seq);
    my @element1=("a","u","g","c","g","u","u","c","c","u","c","g","u","c","c");
    my $putativeCVfound=0;
    #Now we will try to detect those elements Thomas wrote from
    print "<div class='box-header' onclick='toggleBoxContent(this)'>CstF Motif Scan:</div>";
    print "<div class='box-content'>";

    ELEMENT1: for ($count1=0;$count1<@seq-14;$count1++) {
        my $mismatch=0;
        for ($count2=0;$count2<=14;$count2++){
            $mismatch ++ if ($seq[$count1+$count2] ne $element1[$count2]);
            next ELEMENT1 if ($mismatch >=2);
        }
        if ($putativeCVfound == 0) {
            # print "<h4>CstF Motif Hits</h4>\n";
        }

        print "<table class='table-result'>\n";
        print "<tr><th>Element</th><td>Element1</td></tr>\n";
        print "<tr><th>Start Position</th><td>$count1</td></tr>\n";
        print "<tr><th>Mismatches</th><td>$mismatch</td></tr>\n";
        print "</table>\n";

        $putativeCVfound = 1;
    }

    #Next comes Element2, the consensus is: YGUGUYN(0-4)UUYAYUGYGU with 2 mismatches allowed
    ELEMENT2A:  for ($count1=0;$count1<=@seq-20;$count1++) {
        my $ele2a=0;
        $ele2a++ if ($seq[$count1+0] eq 'c' || $seq[$count1+0] eq 'u'); #c or u = y is the first nt of the element2a consensus
        $ele2a++ if ($seq[$count1+1] eq 'g');
        $ele2a++ if ($seq[$count1+2] eq 'u');
        $ele2a++ if ($seq[$count1+3] eq 'g');
        $ele2a++ if ($seq[$count1+4] eq 'u');
        $ele2a++ if ($seq[$count1+5] eq 'c' || $seq[$count1+5] eq 'u');
        my $ele2ab1=$ele2a;
        next ELEMENT2A if ($ele2a<=3); #stoppe wenn nicht mind 4 Punkte vergeben werden
        for (my $count2=0;$count2<=4;$count2++) { #this will allow the 4 possible N(0-4)
            $ele2a=$ele2ab1;
            $ele2a++ if ($seq[$count1+$count2+6] eq 'u');
            $ele2a++ if ($seq[$count1+$count2+7] eq 'u');
            $ele2a++ if ($seq[$count1+$count2+8] eq 'c' || $seq[$count1+$count2+8] eq 'u');
            $ele2a++ if ($seq[$count1+$count2+9] eq 'a');
            $ele2a++ if ($seq[$count1+$count2+10] eq 'c' || $seq[$count1+$count2+10] eq 'u');
            $ele2a++ if ($seq[$count1+$count2+11] eq 'u');
            $ele2a++ if ($seq[$count1+$count2+12] eq 'g');
            $ele2a++ if ($seq[$count1+$count2+13] eq 'c' || $seq[$count1+$count2+10] eq 'u');
            $ele2a++ if ($seq[$count1+$count2+14] eq 'g');
            $ele2a++ if ($seq[$count1+$count2+15] eq 'u');

            # Create table header only once for Element2 results
            if (($ele2a >= 14) && $putativeCVfound == 0) {
                print "<table class='table-result'>\n";
                print "<tr><th><b>CstF:</b></th><th>Start</th><th>Mismatch</th></tr>\n";
                $putativeCVfound = 1;
            }
            printf "<tr><td>Element2a</td><td>%d</td><td>%d</td></tr>\n", $count1, (16 - $ele2a) if ($ele2a >= 14);
        }
    }

    #Next comes Element2b, the consensus is: UUGYUN(0-4)AUUUACU(U/G)N(0-2)YCU with 2 mismatches allowed
    ELEMENT2B: for ($count1=0;$count1<@seq-23;$count1++) {
        my $ele2b=0;
        $ele2b++ if ($seq[$count1+0] eq 'u'); #c or u = y is the first nt of the element2a consensus
        $ele2b++ if ($seq[$count1+1] eq 'u');
        $ele2b++ if ($seq[$count1+2] eq 'g');
        $ele2b++ if ($seq[$count1+3] eq 'c' || $seq[$count1+3] eq 'u');
        $ele2b++ if ($seq[$count1+4] eq 'u');
        my $ele2bb1=$ele2b;
        next ELEMENT2B if ($ele2b<=2); #stoppe wenn nicht mind 3 Punkte vergeben werden
        ELEMENT2BB1: for ($count2=0;$count2<=4;$count2++) { #element 2b (b)reakpunkt 1
            $ele2b=$ele2bb1;
            $ele2b++ if ($seq[$count1+$count2+5] eq 'a');
            $ele2b++ if ($seq[$count1+$count2+6] eq 'u');
            $ele2b++ if ($seq[$count1+$count2+7] eq 'u');
            $ele2b++ if ($seq[$count1+$count2+8] eq 'u');
            $ele2b++ if ($seq[$count1+$count2+9] eq 'a');
            $ele2b++ if ($seq[$count1+$count2+10] eq 'c');
            $ele2b++ if ($seq[$count1+$count2+11] eq 'u');
            $ele2b++ if ($seq[$count1+$count2+12] eq 'u' || $seq[$count1+$count2+12] eq 'g');
            next ELEMENT2BB1 if ($ele2b <=11);
            $ele2bb2=$ele2b;
            for ($count3=0;$count3<=2;$count3++) {
                $ele2b=$ele2bb2;
                $ele2b++ if ($seq[$count1+$count2+$count3+13] eq 'c' || $seq[$count1+$count2+$count3+0] eq 'u');
                $ele2b++ if ($seq[$count1+$count2+$count3+14] eq 'c');
                $ele2b++ if ($seq[$count1+$count2+$count3+15] eq 'u');
                
                # Element2b uses the same table as Element2a (shared $putativeCVfound flag)
                if (($ele2b >= 14) && $putativeCVfound == 0) {
                    print "<table class='table-result'>\n";
                    print "<tr><th><b>CstF:</b></th><th>Start</th><th>Mismatch</th></tr>\n";
                    $putativeCVfound = 1;
                }
                printf "<tr><td>Element2b</td><td>%d</td><td>%d</td></tr>\n", $count1, (16 - $ele2b) if ($ele2b >= 14);
            }
        }
    }

    if ($putativeCVfound == 1) {
        print "</table>\n";
    }

    @seq=();
    print " Those elements are an indication for a processing protein binding motif<br>" if ($putativeCVfound==1);
    print "<div class='info-warning'>";
    print "No CstF Motif found.";
    print "</div>";

    print "</div>";
}

sub stemggpairs {
    print "<div class='box-header' onclick='toggleBoxContent(this)'>Stem GG pairs:</div>";
    print "<div class=box-content>";
=======
					#print "<br><br>Possible hit detected at nt: $transcelegansvalues[$count*10+0] (pointing to the ggua of stem1)<br>";
					print " C.elegans:<br>  Position:   $transcelegansvalues[$count*10+0]  pointing to ggua of stem1<br>";
					
					print "  Stem1:      $transcelegansvalues[$count*10+1]<br>";
					print "  Structure:  $transcelegansvalues[$count*10+2]<br>";
					print "  Stem2:      $transcelegansvalues[$count*10+3]<br>";
					print "  Structure:  $transcelegansvalues[$count*10+4]<br>";
					print "  Sm-Site:    $transcelegansvalues[$count*10+5]<br>";
					print "  Stem3:      $transcelegansvalues[$count*10+6]<br>";
					print "  Structure:  $transcelegansvalues[$count*10+7]<br>";
					print "  Leader:     $transcelegansvalues[$count*10+8]<br>" if ($transcelegansvalues[$count*10+8] !=0);
					print "  Leader:     none<br>" if ($transcelegansvalues[$count*10+8] == 0);
					
					
					
					#print "<table border=1>\n<tr>\n";
					#print '<td align="center"><font face="monospace">Stem1:</td> <td align="center"><font face="monospace">'."$transcelegansvalues[$count*10+1]".'</td></tr>';
					#print '<tr><td align="center"><font face="monospace">Structure:</td><td align="center"><font face="monospace">'."$transcelegansvalues[$count*10+2]".'</td></tr>';
					#print '<tr><td align="center"><font face="monospace">Stem2:</td><td align="center"><font face="monospace">'."$transcelegansvalues[$count*10+3]".'</td></tr><tr><td align="center"><font face="monospace">Structure</td><td align="center"><font face="monospace">'."$transcelegansvalues[$count*10+4]".'</td></tr>';
					#print '<tr><td align="center"><font face="monospace">Sm-Site:</td><td align="center"><font face="monospace">'."$transcelegansvalues[$count*10+5]".'</td></tr>';
					#print '<tr><td align="center"><font face="monospace">Stem3:</td><td align="center"><font face="monospace">'."$transcelegansvalues[$count*10+6]".'</td></tr><tr><td align="center"><font face="monospace">Structure</td><td align="center"><font face="monospace">'."$transcelegansvalues[$count*10+7]".'</td></tr>';
					#print '<tr><td align="center"><font face="monospace">Potential leader:</td><td align="center">'."$transcelegansvalues[$count*10+8]".'</td></tr>' if ($transcelegansvalues[$count*10+8] !=0);
					#print '</table><br>';
					@transsplicing=(@transsplicing,$transcelegansvalues[$count*10+0]-21,($transcelegansvalues[$count*10+9]+(length $transcelegansvalues[$count*10+5])+(length $transcelegansvalues[$count*10+6])));
				}
			}
		}
		
		if ($DOIRE==1){
			#So jetzt machen wir den gleichen Spass wie mit trans
			@irereturnvalues=RNASERVER::IRE::suboptimalfindire($SEQUENCECHECKED);
			#print "<b><big>Iron-responsive element</b></big> search:";
			#print "<b>Iron-resp Ele.:</b> none";
			$irelineprintout=0;
			if (@irereturnvalues>1){
				########new #####
				#$tempp=@irereturnvalues;
				
				$posprintout=0;
				for ($count=0;$count<=@irereturnvalues-1;$count=$count+7){
				    print "<b>Iron-resp Ele.:</b><br>" if ($irelineprintout==0);
				    $irelineprintout=1;
				    
				    if ($posprintout!=$irereturnvalues[$count+0]) {
				    @ire=(@ire,$irereturnvalues[$count+0]-16,$irereturnvalues[$count+0]+22) if ($posprintout!=$irereturnvalues[$count+0]);
				    print " Position:     $irereturnvalues[$count+0]<br>" if ($posprintout!=$irereturnvalues[$count+0]);
				    #print " Quality:      $irereturnvalues[$count+1]<br>";
				   		
				    print " Sequence:     $irereturnvalues[$count+2]<br>" if ($posprintout!=$irereturnvalues[$count+0]); 
			        }
				    print " Structure:    $irereturnvalues[$count+3]";
			    
				    printf ("  %2f kcal/mol ",$irereturnvalues[$count+4]);
				    #print "Quality: good $irereturnvalues[$count+5] $irereturnvalues[$count+6] <br>" if ($irereturnvalues[$count+1] == 1);
				    #print "Quality: bad $irereturnvalues[$count+5] $irereturnvalues[$count+6] <br>" if ($irereturnvalues[$count+1] ==2);
				    print "Quality: good <br>" if ($irereturnvalues[$count+1] == 1);
				    print "Quality: bad <br>" if ($irereturnvalues[$count+1] ==2);

				    $posprintout=$irereturnvalues[$count+0]; 
				 }
					
				########endnew ####
				
			}
			else{
				#print "Nothing found<br>";
			    print "<b>Iron-resp Ele.:</b> none";
			}
			
		}
	}
	else {
		&errors;
	}
}




#########################################################
#########################################################

sub errors {
	print $rv->p("\nSorry, either you did not enter a valid sequence or I could not recognize it!\n") if ($passsequence==0 && $errorsyetprintout==0);
	print $rv->p("\nSorry, you did not choose a consensus pattern!\n") if ($passconsensus==0 && $errorsyetprintout==0);
	print $rv->p("\nPlease go back and correct it!\n") if ($errorsyetprintout==0);;
	$errorsyetprintout=1;
}

#########################################################
#########################################################

sub readandcheckinput {
	$passsequence=1;		#they will be 0 if an error occured here
	$passconsensus=1;
	$xxx='OFF';
	$xxx=~/(OFF)/;  #what the hell is this doing? This shall change the register $1 from ON that it got from IRE to something undefined, here freely choosen off!!!
	#Working with the IRE-check box	
	$IRE=$rv->param("IRE");
	$IRE=~/^(ON)/;
	$IRECHECKED=$1;
	if ($IRECHECKED eq 'ON') {
	    $DOIRE=1;
	}
	else {
	    $DOIRE=0;
	}
	
	#print $rv->p("\nIRE: $IRE\n");
	
	#Working with the TRANS-check box
	$xxx='OFF';
	$xxx=~/(OFF)/;  #what the hell is this doing? This shall change the register $1 from ON that it got from IRE to something undefined, here freely choosen off!!!
	$TRANS=$rv->param("TRANS");
	$TRANS=~/^(ON)/;
	$TRANSCHECKED=$1;
	if ($TRANSCHECKED eq 'ON') {
	    $DOTRANS=1;
	}
	else {
	    $DOTRANS=0;
	}
	print $rv->p("\nTRANS: $TRANS\n") if ($debug);
	
	$passconsensus=0 if ($DOIRE==0 && $DOTRANS==0);
	print $rv->p("\nPASSconsensus: $passconsensus\n") if ($debug);

	#Working with the OPPOSITE checked box
	$xxx='OFF';
	$xxx=~/(OFF)/;  #what the hell is this doing? This shall change the register $1 from ON that it got from IRE to something undefined, here freely choosen off!!!	
	$OPPOSITE=$rv->param("OPPOSITE");
	print $rv->p("\nOPPOSITE: $OPPOSITE\n") if ($debug);
	$OPPOSITE=~/^(YES)/;
	$OPPOSITECHECKED=$1;
	if ($OPPOSITECHECKED eq 'YES') {
	    $DOOPPOSITE=1;
	}
	else {
	    $DOOPPOSITE=0;
	}
	$xxx='OFF';
	$xxx=~/(OFF)/;  #what the hell is this doing? This shall change the register $1 from ON that it got from IRE to something undefined, here freely choosen off!!!		$xxx=~/(OFF)/;  #what the hell is this doing? This shall change the register $1 from ON that it got from IRE to something undefined, here freely choosen off!!!	

	$SEQNAME=$rv->param("SEQNAME");
	$SEQNAME=~/^([ a-zA-Z0-9><;+-]+)/;
	$SEQNAMECHECKED=$1;
	$SEQNAMECHECKED="Your sequence" if ($SEQNAMECHECKED eq "" || $SEQNAMECHECKED eq "ON" | $SEQNAMECHECKED eq "OFF");
	
	print $rv->p("\nSEQNAME: $SEQNAMECHECKED\n") if ($debug);
	$xxx='OFF';
	$xxx=~/(OFF)/;  #what the hell is this doing? This shall change the register $1 from ON that it got from IRE to something undefined, here freely choosen off!!!		$xxx=~/(OFF)/;  #what the hell is this doing? This shall change the register $1 from ON that it got from IRE to something undefined, here freely choosen off!!!	

	#check and correct SEQUENCE
	
	#knocking out the next line due to the use of @FASTASEQ !!!
	#$SEQUENCE=$rv->param("SEQUENCE");
	
	
	# First check if it is in a fasta-format!!!
	#Set $2 to a special value!!
	$temp1='XNOFASTA';
	$temp1=~/^(X)(NOFASTA)/; #So $2 now contains NOFASTA



#################################################################
#	$SEQUENCE=~/^(>[ ]?([A-Za-z0-9: ]+))/;
#	$removefastatag=$1;
#	$SEQNAMECHECKED=$2 if ($2 ne 'NOFASTA');
#	$SEQUENCE=~s/$removefastatag// if ($2 ne 'NOFASTA');
	

## changed into following scripts, by liang
        $SEQUENCE=~/^(>(.*))\n/;
        $removefastatag=$1;
        $SEQNAMECHECKED=$2;
        substr($SEQUENCE,index($SEQUENCE,$removefastatag),length($removefastatag))="";



#Now if it has been fasta, the name has been saved to SEQNAMECHECKED and the tag has been removed	
	$SEQUENCE=~s/[^a-zA-Z]+//g;
	$SEQUENCE=lc($SEQUENCE);
	
	#Read in the origin-field !!
	$ORIGIN=$rv->param("ORIGIN");
	if ($ORIGIN eq 'ON') { #is RNA
		$ORIGINchecked=1;
	}
	else {
		$ORIGINchecked=0; #no RNA (is then DNA)
	}
	
	
	
	
	#Determine if DNA or RNA
	$ucount=0;$tcount=0;$dnarna='';
	while ($SEQUENCE=~m/u/g){
		$ucount++;
	}
	while ($SEQUENCE=~m/t/g){
		$tcount++;
	}
	if ($ucount>$tcount){
		$dnarna='RNA';
	}
	else {
		$dnarna='DNA';
	}
	
	$dnarna='unknown' if ($ucount>20 && $tcount>20); 
	$dnarna='RNA' if ($ORIGINchecked==1);
	#End of DNA RNA detection;
	$SEQUENCE=~s/[^agctu]/n/g;
	$SEQUENCE=~s/t/u/g;
	$SEQUENCE=~/([acgun]+)/; #this ist due to the tainting fuction for the webserver
	$SEQUENCECHECKED=$1;
	$passsequence=0 if ($SEQUENCECHECKED eq '');
	$SEQUENCELENGTH=length $SEQUENCECHECKED;	
        print "<br><b><big>ATTENTION: The length of your sequence is 0, please check your input for errors</b></big><br>" if ($SEQUENCELENGTH==0); 	
	print $rv->p("\nSEQUENCE: $SEQUENCECHECKED\n") if ($debug);
	print $rv->p("\nPASSsequence: $passsequence\n") if ($debug);
}

##################################################
##################################################

sub oppositestrand {
	$SEQUENCECHECKED=reverse ($SEQUENCECHECKED);
	$SEQUENCECHECKED=uc($SEQUENCECHECKED);
	$SEQUENCECHECKED=~s/C/g/g;
	$SEQUENCECHECKED=~s/G/c/g;
	$SEQUENCECHECKED=~s/A/u/g;
	$SEQUENCECHECKED=~s/U/a/g;
	$SEQUENCECHECKED=lc($SEQUENCECHECKED);
}

##################################################
##################################################


sub additionalinformation {
  
	#prints out the length and origin of the pasted sequence!	
	print "<pre>\n";
	$dnarna='unknown' if ($dnarna eq '');
		#print $rv->p("You submitted a $SEQUENCELENGTH nt long sequence of $dnarna origin.");
		
	#For the picture creation	
	open (SEQPIC,">$TEMPDIR/$job.seq"); #don't know if this works!!	
	print SEQPIC ">$job\n$SEQUENCECHECKED\n"; #shall create a fasta-format sequence file !!!
	close SEQPIC;
  
	#Unfortunately, genscan only accepts t's ! It doesn't accept u's ! So, we'll create its own file!

	open (GENSCAN,">$TEMPDIR/$job.genscan");
	$SEQGENSCAN=$SEQUENCECHECKED;
	$SEQGENSCAN=~s/u/t/g;
	print GENSCAN ">$job\n$SEQGENSCAN\n";
	close GENSCAN;
	$SEQGENSCAN='';
	
	#Create the genscan statistics output!
	@answergenscan=`$GENSCANFILE/genscan $GENSCANLIB/HumanIso.smat $TEMPDIR/$job.genscan >$TEMPDIR/$job.genscanout`;
        #    -v -cds	gives more infos for debugging
	
	chdir $TEMPDIR;	

	#Create a picture!!!!
	if (length $SEQUENCECHECKED<=$MAXFOLDINGLEN) { #we should not fold sequences larger than this !!!
		
		##############################
		&createpicture;
		##############################
		
		##############################
		&checkstems; #this has to be run soon after &createpicture (?? has it ??)

		##############################

		#############################################################
		#Now looking for gg pairs in a stem that are
		#inserted between 2 bondings
		#call soubroutine &stemggpairs
		#&stemggpairs;
		
	}
	else {
		print "<br><b>Length:</b>        $SEQUENCELENGTH";
		print "     *some information is only available up to $MAXFOLDINGLEN nt</font>\n" if ($SEQUENCELENGTH >$MAXFOLDINGLEN); 
		print "<br><b>Origin:</b>        $dnarna<br>";
	#print $rv->p("Remember: Your sequence will only be folded up to $MAXFOLDINGLEN nts !!");
	}
	
	&promotor; #search for promotor	
	
	&exons;	  # very time consuming, let s take a look 	
		
	############################################################
	# Trying to calculate potential 5' and 3' UTR
	#
	&calcUTR;


	print "<br><b>3' UTR:</b><br>";

	&polyAsignal; #search for polyAsignal

	
	&ARE; #search the ARE


	print "<br><b>Catalytic RNA:</b><br>";

        #########################################################
	#Looking for a sm-site with the consensus:
	# a/g a/g poly-u a/g a/g , with one other nt allowed in the poly-u stretch
	#

	&smsite;

	
	&tRNA; #search for tRNA

	
	$genscan->close(); #was located in front of the ARE
				
	#########################################################
	#Looking for serveral other regulatory signals
	#Cleavage stimulation factor containing elements
	#element1: AUGCGUUCCUCGUCC
	#element 2a: YGUGUYN(0-4)UUYAYUGYGU
	#element 2b: UUGYUN(0-4)AUUUACU(U/G)N(0-2)YCU
	#allow 2 mismatches each!!! therefore you can't use pattern matching

	&csfce; #Subroutine containing the long search program for those sequences
	#WRONG POSITION!!!!!!!! MUST BE OUTSIDE THE 500 nt barrier!!!!
	
	########################################################
	# Looking for the protein A1 binding site C9E
	&protA1bisite;
	if (length $SEQUENCECHECKED<=$MAXFOLDINGLEN) { 
		&stemggpairs;
		print "<br><b>Here is the folding of your sequence:</b><br>";

                # comment out by liang
		#print '<img src="'."/cgi-bin/get_result_picture.pl?"."$job"."_"."$md5sum".'"><br>';
	        
	}
   print"</pre>";


	$md5in=`md5sum $TEMPDIR/$job.seq`;
	$md5in=~/([0-9abcdef]{32})/;
	$md5sum=$1;
	$actualdir=cwd();
	chdir "$TEMPDIR";
	$answerconvert=system("convert rna.ps -crop 0x0 $job"."_$md5sum.jpg");  # for original vienna package
	print '<img src="'."/cgi-bin/get_result_picture.pl?"."$job"."_"."$md5sum".'"><br>';

}


sub jobnumber {
	#Generiere eine Datei, die die Jobnummern handelt und erhoehe die Jobnummer bei jedem Job um 1 !!!
	my $line=0;
	open (JOBNUMMER,"$TEMPDIR/jobnr.dat") || print $rv->p("Could not open file, will create a new one!");
	$line=<JOBNUMMER>; 
	
	close JOBNUMMER;

	open (JOBNUMMER,">$TEMPDIR/jobnr.dat");
	$line++;
	print JOBNUMMER $line;
	#print $rv->p("New jobnumber $line");
	close JOBNUMMER;

	
	return ($line);
}


sub csfce {
	my $count1=0; #We will mark $count1 as my so that we won't have any problems later
	#$seq='llllllllllllllllluugculllauuuacuglcculllaugcguuccucgucclllllllllllllllll';
	my @seq=split ('',$SEQUENCECHECKED);
	#my @seq=split ('',$seq);
	my @element1=("a","u","g","c","g","u","u","c","c","u","c","g","u","c","c");
	my $putativeCVfound=0;
	#Now we will try to detect those elements Thomas wrote from

	ELEMENT1: for ($count1=0;$count1<@seq-14;$count1++) {
        	my $mismatch=0;
	        for ($count2=0;$count2<=14;$count2++){
        	        $mismatch ++ if ($seq[$count1+$count2] ne $element1[$count2]);
                	next ELEMENT1 if ($mismatch >=2);
	        }
        	print "<b>CstF:</b>          start      mismatch<br>" if ($putativeCVfound==0);
		printf (" Element1:     %-6d     %-2d<br>",$count1,$mismatch);
		$putativeCVfound=1;
	}
	#Next comes Element2, the consensus is: YGUGUYN(0-4)UUYAYUGYGU with 2 mismatches allowed
	ELEMENT2A:  for ($count1=0;$count1<=@seq-20;$count1++) {
        	my $ele2a=0;
	        $ele2a++ if ($seq[$count1+0] eq 'c' || $seq[$count1+0] eq 'u'); #c or u = y is the first nt of the element2a consensus
        	$ele2a++ if ($seq[$count1+1] eq 'g');
	        $ele2a++ if ($seq[$count1+2] eq 'u');
	        $ele2a++ if ($seq[$count1+3] eq 'g');
	        $ele2a++ if ($seq[$count1+4] eq 'u');
        	$ele2a++ if ($seq[$count1+5] eq 'c' || $seq[$count1+5] eq 'u');
	        my $ele2ab1=$ele2a;
        	next ELEMENT2A if ($ele2a<=3); #stoppe wenn nicht mind 4 Punkte vergeben werden
	        for (my $count2=0;$count2<=4;$count2++) { #this will allow the 4 possible N(0-4)
        	        $ele2a=$ele2ab1;
                	$ele2a++ if ($seq[$count1+$count2+6] eq 'u');
	                $ele2a++ if ($seq[$count1+$count2+7] eq 'u');
        	        $ele2a++ if ($seq[$count1+$count2+8] eq 'c' || $seq[$count1+$count2+8] eq 'u');
                	$ele2a++ if ($seq[$count1+$count2+9] eq 'a');
	                $ele2a++ if ($seq[$count1+$count2+10] eq 'c' || $seq[$count1+$count2+10] eq 'u');
        	        $ele2a++ if ($seq[$count1+$count2+11] eq 'u');
                	$ele2a++ if ($seq[$count1+$count2+12] eq 'g');
        	        $ele2a++ if ($seq[$count1+$count2+13] eq 'c' || $seq[$count1+$count2+10] eq 'u');
	                $ele2a++ if ($seq[$count1+$count2+14] eq 'g');
                	$ele2a++ if ($seq[$count1+$count2+15] eq 'u');
	                print "<b>CstF:</b>          start      mismatch<br>" if  (($ele2a>=14)&& $putativeCVfound==0);
			printf (" Element2a     %-6d     %-2d<br>",$count1,(16-$ele2a)) if  ($ele2a>=14);
			$putativeCVfound=1 if ($ele2a>=14);
		}
	}
	#Next comes Element2b, the consensus is: UUGYUN(0-4)AUUUACU(U/G)N(0-2)YCU with 2 mismatches allowed
	ELEMENT2B: for ($count1=0;$count1<@seq-23;$count1++) {
        	my $ele2b=0;
	        $ele2b++ if ($seq[$count1+0] eq 'u'); #c or u = y is the first nt of the element2a consensus
        	$ele2b++ if ($seq[$count1+1] eq 'u');
	        $ele2b++ if ($seq[$count1+2] eq 'g');
        	$ele2b++ if ($seq[$count1+3] eq 'c' || $seq[$count1+3] eq 'u');
	        $ele2b++ if ($seq[$count1+4] eq 'u');
        	my $ele2bb1=$ele2b;
	        next ELEMENT2B if ($ele2b<=2); #stoppe wenn nicht mind 3 Punkte vergeben werden
        	ELEMENT2BB1: for ($count2=0;$count2<=4;$count2++) { #element 2b (b)reakpunkt 1
                	$ele2b=$ele2bb1;
	                $ele2b++ if ($seq[$count1+$count2+5] eq 'a');
        	        $ele2b++ if ($seq[$count1+$count2+6] eq 'u');
                	$ele2b++ if ($seq[$count1+$count2+7] eq 'u');
	                $ele2b++ if ($seq[$count1+$count2+8] eq 'u');
        	        $ele2b++ if ($seq[$count1+$count2+9] eq 'a');
                	$ele2b++ if ($seq[$count1+$count2+10] eq 'c');
	                $ele2b++ if ($seq[$count1+$count2+11] eq 'u');
        	        $ele2b++ if ($seq[$count1+$count2+12] eq 'u' || $seq[$count1+$count2+12] eq 'g');
                	next ELEMENT2BB1 if ($ele2b <=11);
	                $ele2bb2=$ele2b;
        	        for ($count3=0;$count3<=2;$count3++) {
                	        $ele2b=$ele2bb2;
                        	$ele2b++ if ($seq[$count1+$count2+$count3+13] eq 'c' || $seq[$count1+$count2+$count3+0] eq 'u');
	                        $ele2b++ if ($seq[$count1+$count2+$count3+14] eq 'c');
        	                $ele2b++ if ($seq[$count1+$count2+$count3+15] eq 'u');
                	        print "<b>CstF:</b>          start      mismatch<br>" if  (($ele2b>=14)&& $putativeCVfound==0);
				printf (" Element2b:    %-6d     %-2d<br>",$count1,(16-$ele2b)) if  ($ele2b>=14);
	                	$putativeCVfound=1;
			}
        	}
	}
	@seq=();
	print " Those elements are an indication for a processing protein binding motif<br>" if ($putativeCVfound==1);

}






sub stemggpairs {
>>>>>>> origin/master
    my @sequ=split('',$SEQUENCECHECKED);
    my $str=join ('',@structure);
    my $stemggpairfound=0;
    my $stemggleadlineprinted=0;
    @stemggpairs=();
    while ($str=~/\(\(\.\(\(/g) {
        my $anfang=(pos $str);
        if ($sequ[$anfang-3] eq 'g') {    
		#print "Anfang: $anfang <br>";
            while ($str=~/\)\)\.\)\)/g) {
                my $ende=(pos $str);
                if ($sequ[$ende-3] eq 'g') {            
			#print "Ende: $ende <br>";
                    #Now compare the amount of paired!!
                    my $stringtocompare=substr($str,$anfang-5,$ende-$anfang+5);
		    #print "StrtoComp: $stringtocompare <br>";
                    my $seqtocompare=substr($seq,$anfang-5,$ende-$anfang+5);
		    #print "SeqtoComp: $seqtocompare <br>";
                    my $openkl=0;
                    my $closedkl=0;
                    while ($stringtocompare=~/([()])/g ) {
                        $openkl++ if ($1 eq '(');
                        $closedkl++ if ($1 eq ')');
                    }
                    if ($openkl==$closedkl && $ende-$anfang<=50) {
<<<<<<< HEAD
                        # Create table 
                        if ($stemggleadlineprinted == 0) {
                            print "<table class='table-result'>\n";
                            print "<tr><th><b>StemGGpair:</b></th><th>Start</th><th>End</th></tr>\n";
                            $stemggleadlineprinted = 1;
                        }
                        
                        printf "<tr><td>Hit</td><td>%d</td><td>%d</td></tr>\n", $anfang-2, $ende-2;
                        
                        @stemggpairs=(@stemggpairs,$anfang-4,$ende); #these are pointing here
                                                #	((.((    )).))
                        $stemggpairfound=1;		#       ^            ^
                                                #       |            |
                    }
                    
                }
=======
		    	print "<b>StemGGpair:</b>    start     -     end<br>" if ($stemggleadlineprinted==0);
			$stemggleadlineprinted=1;
			
			printf (" Hit:          %-5d     -  %5d<br>",$anfang-2,$ende-2);
                        @stemggpairs=(@stemggpairs,$anfang-4,$ende); #these are pointing here
								#	((.((    )).))
			$stemggpairfound=1;			#       ^            ^
								#       |            |
		    }
		}
>>>>>>> origin/master
                pos($str)=pos($str)-2;
            }        
        }
        pos($str)=$anfang-2;
    }
<<<<<<< HEAD
    if ($stemggleadlineprinted == 1) {
                        print "</table>\n";
                    }

    print "<h3>No Stem GG Pair Found.</h3><br>" if ($stemggpairfound==0); 
    @sequ=();
    $str='';

    print "</div>";
}

sub rbp {
    print "<div class='box-header' onclick='toggleBoxContent(this)'>RNA Binding Protein Motif(s) Scan:</div>";
    print "<div class='box-content'>";

    my $fimo_outfile = "$TEMPDIR/fimo.txt";
    my $input_seq = "$TEMPDIR/$job.seq";  # path to your .meme file

    # Run FIMO
    my $fimo_cmd = "$FIMO/fimo --text --thresh 1e-4 $RBPDB $input_seq > $fimo_outfile";
    system($fimo_cmd);

    my @fimo_results;
    @rbp_locs = ();  # global array to store start/end positions

    # Parse FIMO output
    if (-e $fimo_outfile) {
        open my $fh, '<', $fimo_outfile or die "Cannot open FIMO output: $!";
        while (my $line = <$fh>) {
            next if $. == 1;
            chomp $line;
            my ($motif_id, $alt_id, $seq_name, $start, $stop, $strand, $score, $pvalue, $qvalue, $matched_seq) = split /\t/, $line;

            my ($protein_name, $motif_name);

            if ($alt_id && index($alt_id, '_') != -1) {
                ($protein_name, $motif_name) = split(/_/, $alt_id, 2);
            } elsif ($motif_id && index($motif_id, '_') != -1) {
                ($protein_name, $motif_name) = split(/_/, $motif_id, 2);
            } else {
                $protein_name = $alt_id || $motif_id;
                $motif_name = '';
            }

            push @fimo_results, {
                motif_id     => $motif_id,
                alt_id       => $alt_id || $motif_id,
                protein      => $protein_name,
                motif        => $motif_name,
                seq_name     => $seq_name,
                start        => $start,
                end          => $stop,
                # strand       => $strand,
                score        => $score,
                pvalue       => $pvalue,
                matched_seq  => $matched_seq,
            };
        }
        close $fh;
    }

    # remove reduntant hits for each protein+domain combination
    my %best_hits;
    
    foreach my $hit (@fimo_results) {
        my $combination = $hit->{protein} . "\t" . $hit->{motif};
        
        if (!exists $best_hits{$combination} || $hit->{score} > $best_hits{$combination}->{score}) {
            $best_hits{$combination} = $hit;
        }
    }
    
    my @filtered_results = sort { $a->{start} <=> $b->{start} } values %best_hits;
    
    # return locations for annoation
    @rbp_locs = ();
    foreach my $hit (@filtered_results) {
        push @rbp_locs, $hit->{start}, $hit->{end};
    }


    # print results
    my $total = scalar  @filtered_results;

    if ($total > 0) {
        print "<table class='table-result'>\n";
        print "<tr><th>Protein</th><th>Motif</th><th>Start</th><th>End</th><th>Score</th><th>p-value</th><th>Matched Sequence</th></tr>\n";

        foreach my $hit (@filtered_results) {
            print "<tr>";
            print "<td>$hit->{protein}</td>";
            print "<td>$hit->{motif}</td>";
            print "<td>$hit->{start}</td>";
            print "<td>$hit->{end}</td>";
            # print "<td>$hit->{strand}</td>";
            print "<td>$hit->{score}</td>";
            print "<td>$hit->{pvalue}</td>";
            print "<td><tt>$hit->{matched_seq}</tt></td>";
            print "</tr>";
        }

        print "</table>";
    } else {
        print "<b>No Protein Binding Motif found above threshold.</b>\n";
    }

    print "</div>";
}

sub ARE {
    $arepresent=0;
    $are_pos=1;
    my $are_table_opened = 0;
    #Check for so called ARE = Au-rich regions; consensus (AUUUA)n of ~50 bases
    print "<div class='box-header' onclick='toggleBoxContent(this)'>Au-rich regions:</div>";

    print "<div class='box-content'>";
    while ($SEQUENCECHECKED=~/([ag]uuu[ag](uuu[ag])+)/g) {
        $are_len=length($1);
        $are_pos=pos($SEQUENCECHECKED)-$are_len;
        if ($are_len >=9){
            # Open table only once when first ARE is found
            if ($are_table_opened == 0) {
                print "<table class='table-result'>\n";
                print "<tr><th><b>ARE</b></th><th>Start</th><th>End</th><th>Sequence</th><th>Mismatches</th></tr>\n";
                $are_table_opened = 1;
            }
            
            $arepresent=1;
            $mismatchinare=0;
            @aretemp=split('',$1);
            for ($arecount=0;$arecount<@aretemp;$arecount++){
                $mismatchinare++ if ($aretemp[$arecount] eq 'g');
            }
            
            printf "<tr><td>Hit</td><td>%d</td><td>%d</td><td>%s</td><td>%d</td></tr>\n",
                   $are_pos, $are_pos+$are_len-1, $1, $mismatchinare;
            @aurichregion=(@aurichregion,$are_pos+1,$are_pos+$are_len);
        }
    }

    # Close table if it was opened
    if ($are_table_opened == 1) {
        print "</table>\n";
    }

    if ($arepresent==0) {
        print "<div class='info-warning'>";
        print "No Au-rich region found.";   #     *(AU-rich region of at least 30 nt)<br>";
        print "</div>";
    }
    print "</div>";
}

sub tRNA {
    # Looking for tRNAs using tRNAscan-SE

    print "<div class='box-header' onclick='toggleBoxContent(this)'>tRNA Scan:</div>";
    print "<div class='box-content'>";
    my $trnascan_output = "$TEMPDIR/$job.trnascanout";
    # my $trnascan_out    = "$TEMPDIR/$job.trnascanlog";
    
    my $trnascan = "$TRNASCANFOLDER/tRNAscan-SE -Q -y -f $trnascan_output $TEMPDIR/$job.seq";
    
    system($trnascan);
        
    my @results;
    @trna_loc = ();
    
    open my $fh_trna, '<', $trnascan_output or die "Cannot open tRNA result: $!";
    while (my $line = <$fh_trna>) {
        chomp $line;
        next if $line =~ /^$/;  # Skip empty lines
        
        # Parse the first line with location and length info
        if ($line =~ /^(\S+)\s+\((\d+)-(\d+)\)\s+Length:\s+(\d+)\s+bp/) {
            my ($name, $from, $to, $length) = ($1, $2, $3, $4);
            
            push @trna_loc, $from, $to;
            
            # Read the next line for Type, Anticodon, and Score
            my $info_line = <$fh_trna>;
            chomp $info_line if $info_line;
            
            my ($type, $anticodon, $anticodon_pos, $score) = ('', '', '', '');
            if ($info_line && $info_line =~ /Type:\s+(\S+)\s+Anticodon:\s+(\S+)\s+at\s+(\d+-\d+)\s+.*Score:\s+([\d.]+)/) {
                ($type, $anticodon, $anticodon_pos, $score) = ($1, $2, $3, $4);
            }
            
            push @results, {
                name          => $name,
                from          => $from,
                to            => $to,
                length        => $length,
                type          => $type,
                anticodon     => $anticodon,
                anticodon_pos => $anticodon_pos,
                score         => $score + 0,  # force numeric
            };
            
            # Skip the sequence and structure lines
            <$fh_trna>;  # Skip separator line
            <$fh_trna>;  # Skip sequence line
            <$fh_trna>;  # Skip structure line
        }
    }
    close $fh_trna;
    
    my $total = scalar @results;
    
    if ($total > 0) {
        print "<table class='table-result'>\n";
        print  "<tr><th>Name</th><th>Start</th><th>End</th><th>Length</th><th>Type</th><th>Anticodon</th><th>Anticodon Pos</th><th>Score</th></tr>\n";
        # printf $format, "Name", "Position", "Length", "Type", "Anticodon", "Anticodon Pos", "Score";
        # print "-" x 85, "\n";
        
        # Sort by score descending (higher scores are better)
        @results = sort { $b->{score} <=> $a->{score} } @results;
        
        foreach my $hit (@results) {
            # my $position = "$hit->{from}-$hit->{to}";
                print "<tr>"; 
                print "<td>$hit->{name}</td>"; 
                print "<td>$hit->{from}</td>"; 
                print "<td>$hit->{to}</td>";
                print "<td>$hit->{length}</td>"; 
                print "<td>$hit->{type}</td>"; 
                print "<td>$hit->{anticodon}</td>"; 
                print "<td>$hit->{anticodon_pos}</td>"; 
                print "<td>$hit->{score}</td>";
                print "<tr>";
        }

        print "</table>";
    } else {
        print "<div class='info-warning'>";
        print "No region matching a tRNA was found.";
        print "</div>";
    }
    
    print "</div>";
    # print "DEBUG: @trna_loc" if @trna_loc;
}

sub smsite {
    print "<div class='box-header' onclick='toggleBoxContent(this)'>Catalytic RNA site(s) scan:</div>";
    print "<div class='box-content'>";
    $smlength=-1;
    $smpos=-1;
    $leadlineprinted=0; #the first line not yet printed
    my $table_opened = 0;
    
    while ($SEQUENCECHECKED=~/([ag][ag](u+([agc]?)u+)[ag][ag])/g){
        $smlength=length $1;
        $smpos=pos($SEQUENCECHECKED)-$smlength+1;
        
        # Open table only once when first motif is found
        if ($table_opened == 0) {
            print "<table class='table-result'>\n";
            print "<tr><th><b>snRNP-motifs:</b></th><th>Start</th><th>Sequence</th><th>Quality</th></tr>\n";
            $table_opened = 1;
        }

        if (length $3 == 1 && length $2 >=4) {
            printf "<tr><td>snRNP-motif</td><td>%d</td><td>%s</td><td>+</td></tr>\n", $smpos, $1;
            @smsite=(@smsite,$smpos,$smpos+((length $1)-1));
            $leadlineprinted=1;
        }
        if (length $3 == 0 && length $2>=4) {
            printf "<tr><td>Put. sm-site</td><td>%d</td><td>%s</td><td>++</td></tr>\n", $smpos, $1;
            @smsite=(@smsite,$smpos,$smpos+((length $1)-1));
            $leadlineprinted=1;
        }
    }
    
    # Close table if it was opened
    if ($table_opened == 1) {
        print "</table>\n";
    }
    
    print "<div class='info-warning'><b>No snRNP-motifs found.</b></div>" if ($leadlineprinted==0);
    
    ##### OUTPUT if Seq is RNA has no cds but smsite --> structured perhaps catalytic RNA !

    if ($grepanswer=~/NO EXONS/ && $ORIGINchecked==1 && @smsite>0){
        print "<br>As I could not detect a coding sequence on this RNA, but there are 1 or more sn-RNP motifs (sm-sites),<br>it might be possible that this is a catalytic RNA!!<br>";
    }
    print "</div>";
}

sub createfolding {
    # comment out by liang
    # new ViennaRNA does not have the old FOLD program but rather incorporated in the same program RNAfold AA
		my $infile = "$TEMPDIR/$job.seq";
        my $outfile = "$TEMPDIR/$job.foldout";

        die "Infile not found: $infile" unless -e $infile;
        die "RNAfold binary not found: $VIENNARNAFOLDDIR/RNAfold" unless -x "$VIENNARNAFOLDDIR/RNAfold";
        my $cmd = "cat $infile | $VIENNARNAFOLDDIR/RNAfold > $outfile 2>&1";
        system("$cmd") == 0 or die "RNAfold failed";
        # write_file("$TEMPDIR/debug_command.sh", "#!/bin/bash\n$cmd\n");
        # chmod 0755, "$TEMPDIR/debug_command.sh";
        open(my $fh, '<', $outfile) or die "Can't open output";
        <$fh>;  # skip FASTA header
        my $seq = <$fh>;
        my $struct_line = <$fh>;

        $struct_line =~ /([().]+)\s+\(([-\d.]+)\)/;
        @structure = split('', $1);
        $structure = $1;
        $energy = $2;


        my $svg_file = "$TEMPDIR/${SEQNAMECHECKED}_ss.svg";

        # Make sure RNAplot output is ready
        system("$VIENNARNAFOLDDIR/RNAplot --infile=$TEMPDIR/$job.foldout -f svg --filename-full");

        # Read SVG file content
        open(my $svgfh, '<', $svg_file) or die "Cannot open SVG file: $!";
        my $svg_content = do { local $/; <$svgfh> };
        close($svgfh);

        # Add ID to <svg> tag if not present
        $svg_content =~ s/<svg /<svg id="static_rna_ss" width="auto" height="600" style="background:white;" /;

        # Output
        print "<h3>RNA Structure (Annotated Structure Below):</h3>\n";
        print "$svg_content\n";
        print "<p style='font-size: 0.9em; color: gray; text-align: center;'> Drag to pan, scroll to zoom</p>\n";

        # Add svg-pan-zoom script
        print "<script src='/js/svg-pan-zoom.min.js'></script>\n";
        print "<script>\n";
        print "  svgPanZoom('#static_rna_ss', {\n";
        print "    zoomEnabled: true,\n";
        print "    controlIconsEnabled: true,\n";
        print "    fit: true,\n";
        print "    center: true\n";
        print "  });\n";
        print "</script>\n";
		
}

sub checkstems {

    # new version now checks 3 bp pairing 
    # added checks for hairpins (only terminal?)
    # also added stacks to not miss or recount the stems and/or hairpins
    # finally also added bulge incorporation without breaking stem pairs
    # need external input before final publishing

    my @pairing = ();
    my %visited;
    my @stack;
    my $stem_count = 0;
    my $hairpin_count = 0;
    my $comment;
    
    # Build pairing map
    for (my $i = 0; $i < @structure; $i++) {
        if ($structure[$i] eq '(') {
            push @stack, $i;
        } elsif ($structure[$i] eq ')') {
            if (@stack) {
                my $j = pop @stack;
                $pairing[$i] = $j;
                $pairing[$j] = $i;
            }
        }
    }
    
    # Find stems
    my @stems = ();
    
    for (my $i = 0; $i < @structure; $i++) {
        next if $visited{$i};
        next unless defined $pairing[$i];
        next if $pairing[$i] < $i;  # Process each pair only once
        
        my $left = $i;
        my $right = $pairing[$i];
        
        # Minimum loop length check
        my $loop_length = $right - $left - 1;
        next if $loop_length < 3;
        
        my $count = 1;
        my @positions = ($left, $right);
        
        # Extend stem with bulge tolerance
        my ($l, $r) = ($left, $right);
        while (1) {
            my $found = 0;
            
            # Try consecutive pairing first
            my $nl = $l + 1;
            my $nr = $r - 1;
            
            if ($nl < $nr && defined $pairing[$nl] && $pairing[$nl] == $nr && 
                !$visited{$nl} && !$visited{$nr}) {
                $count++;
                $l = $nl;
                $r = $nr;
                push @positions, $nl, $nr;
                $found = 1;
            }
            # Try with 1-nucleotide bulge
            elsif (!$found) {
                # Try bulge on left side
                $nl = $l + 2;
                $nr = $r - 1;
                if ($nl < $nr && defined $pairing[$nl] && $pairing[$nl] == $nr && 
                    !$visited{$nl} && !$visited{$nr}) {
                    $count++;
                    $l = $nl;
                    $r = $nr;
                    push @positions, $nl, $nr;
                    $found = 1;
                }
                # Try bulge on right side
                elsif (!$found) {
                    $nl = $l + 1;
                    $nr = $r - 2;
                    if ($nl < $nr && defined $pairing[$nl] && $pairing[$nl] == $nr && 
                        !$visited{$nl} && !$visited{$nr}) {
                        $count++;
                        $l = $nl;
                        $r = $nr;
                        push @positions, $nl, $nr;
                        $found = 1;
                    }
                }
            }
            
            last unless $found;
        }
        
        # Check minimum stem length >= 3
        if ($count >= 3) {
            $stem_count++;
            
            # Store stem information
            push @stems, {
                left => $left,
                right => $right,
                length => $count,
                positions => [@positions]
            };
            
            # Mark positions as visited
            $visited{$_} = 1 for @positions;
        }
    }
    
    # Better hairpin detection
    for my $stem (@stems) {
        my $left = $stem->{left};
        my $right = $stem->{right};
        
        # Check if this stem forms a terminal hairpin
        my $has_internal_stems = 0;
        
        for my $other_stem (@stems) {
            next if $other_stem == $stem;
            
            # Check if other stem is completely inside this stem's loop
            if ($other_stem->{left} > $left && $other_stem->{right} < $right) {
                $has_internal_stems = 1;
                last;
            }
        }
        
        # Additional check: ensure the loop region is reasonable for a hairpin
        my $loop_size = $right - $left - 1;
        if (!$has_internal_stems && $loop_size >= 3 && $loop_size <= 30) {
            $hairpin_count++;
        }
    }
    
    # Output results
    print "<table class='table-result'>";
    print "<tr><th>Length</th><td>$SEQUENCELENGTH</td></tr>\n";
    print "<tr><th>Energy</th><td>$energy kcal/mol</td></tr>\n";
    print "<tr><th>Stems</th><td>$stem_count stem structure(s)</td></tr>\n";
    print "<tr><th>Hairpins</th><td>$hairpin_count hairpin(s)</td></tr>\n";
    
    print "</table>";
    # Prediction logic
    my $structure_length = scalar @structure;
    my $avg_spacing = $structure_length / ($stem_count || 1);
    
    if ($stem_count >= 15 && $avg_spacing < 100) {
        $comment = "Highly structured RNA can be likely rRNA, tRNA, or any other regulatory RNA\n";
    } elsif ($hairpin_count >= 1 && $stem_count <= 3) {
        $comment = "Simple structured RNA, possible miRNA, siRNA, or regulatory element\n";
    } elsif ($stem_count >= 1) {
        $comment = "Some secondary structure detected, may have biological significance\n";
    } else {
        $comment = "Minimal secondary structure\n";
    }
    print "<div class='info-info'>$comment</div>";
    
    if ($hairpin_count > 0 && $stem_count > 0 && $hairpin_count / $stem_count > 0.7) {
        print "\t\tHigh hairpin content — characteristic of miRNA precursors\n";
    }
    
    
    print "<br>";
    
    return ($stem_count, $hairpin_count, $energy);
}

sub checkstemsonly {
    # revamped checkstemsonly based on the logic of the new checkstems
    my $inseq = $_[0];
    my $inwhattodo = $_[1];   # check if structure given, 0 = strcuture, 1 = sequence which will be folded
    my $struct;
    my @structure;
    my $energy;
    
    # Handle RNA folding or structure input
    if ($inwhattodo == 1 && length $inseq <= $MAXFOLDINGLENUTR) {
        @struct = `echo $inseq | $VIENNARNAFOLDDIR/RNAfold`; 
        if ($? != 0) { ##added by AA, error handling
            return ("RNAfold failed", 0);
        }
        $struct[1] =~ /([().]+) \(([+-. 0-9]+)\)/;
        @structure = split('', $1);
        $energy = $2;
    }
    if ($inwhattodo == 1 && length $inseq > $MAXFOLDINGLENUTR) {
        return ("Too long for detection", 1); 
    }
    if ($inwhattodo == 0) {
        @structure = split ('', $inseq);
        $energy = 0;
    }

    my @pairing = ();
    my %visited;
    my @stack;
    my $fldstemsauf = 0;
    my $fldstemszu = 0;
    
    # pairng map creation
    for (my $i = 0; $i < @structure; $i++) {
        if ($structure[$i] eq '(') {
            push @stack, $i;
        } elsif ($structure[$i] eq ')') {
            if (@stack) {
                my $j = pop @stack;
                $pairing[$i] = $j;
                $pairing[$j] = $i;
            }
        }
    }
    
    # Scanning for biologically valid stems
    for (my $i = 0; $i < @structure; $i++) {
        next if $visited{$i};
        next unless defined $pairing[$i];
        next if $pairing[$i] < $i;  # avoid recounting the same pair 
        
        my $left = $i;
        my $right = $pairing[$i];
        
        # BIOLOGICAL VALIDATION: Check minimum loop length
        my $loop_length = $right - $left - 1;
        next if $loop_length < 3;  # Skip biologicall imposiible loops
        
        my $count = 1;
        my @positions = ($left, $right);
        
        # Extend stem with improved bulge tolerance
        my ($l, $r) = ($left, $right);
        while (1) {
            my $found = 0;
            
            # Try consecutive pairs first
            my $nl = $l + 1;
            my $nr = $r - 1;
            
            if ($nl < $nr && defined $pairing[$nl] && $pairing[$nl] == $nr && 
                !$visited{$nl} && !$visited{$nr}) {
                $count++;
                $l = $nl;
                $r = $nr;
                push @positions, $nl, $nr;
                $found = 1;
            }
            # Try 1-nt symmetric bulges
            elsif (!$found) {
                $nl = $l + 2;
                $nr = $r - 2;
                
                if ($nl < $nr && defined $pairing[$nl] && $pairing[$nl] == $nr && 
                    !$visited{$nl} && !$visited{$nr}) {
                    $count++;
                    $l = $nl;
                    $r = $nr;
                    push @positions, $nl, $nr;
                    $found = 1;
                }
            }
            
            last unless $found;
        }
        
        # Check minimum stem length (reduced from 5 to 3 for biological relevance)
        if ($count >= 3) {
            $fldstemsauf++;
            $visited{$_} = 1 for @positions;
        }
    }
    
    # Set fldstemszu to same value (maintaining original interface behavior)
    # Your original code had logic issues with separate auf/zu counting
    $fldstemszu = $fldstemsauf;
    
    # Return results in same format as original
    return ($fldstemsauf, $fldstemszu, $energy);
}

sub predprotein {
	# $predprotforAnDom=0;
	print "<div class='box-header' onclick='toggleBoxContent(this)'>Predicted Protein Sequence:</div>";
    print "<div class='box-content'>";
    if (defined $cpc){
        print "<i>Protein is predicted from CPC2 output.</i><br>";
    }
	if (@predprot>1){
		# print "<br>";	
		for ($count1=1;$count1<=@predprot;$count1++) {
			print $predprot[$count1-1];
			# print "<br>" if ($count1%120==0);
		}
	print "<br>";
	}
	else {
        # print "<br>";
        print "<div class='info-warning'>";
		print "No protein sequence was predicted. Check coding potential of the sequence above.";
        print "</div>";
	}

    print "</div>";
}

# microrna search should be full length and should show potential micrornas with a warning. 
sub RNAMOTIF {
    my $tblout_file = "$TEMPDIR/$job.tblout";  # Table format output
	my $output_file = "$TEMPDIR/$job.out";     # Full verbose output

	my $cmd = "$CMSCAN/cmscan -E 0.001 --tblout $tblout_file -o $output_file $RFAM $TEMPDIR/$job.seq > /dev/null 2>&1";
	system($cmd);

	# my $format = "%-12s %-12s %-12s %-6s %-8s %-10s %-20s\n";

	my $found = 0;
    @rna_motif =();
    my @results;
    my $found = 0;
    my $i = 0;

    print "<link rel='stylesheet' href='/css/fornac.css'>";
    print "<script src='/js/d3.v3.min.js'></script>";
    print "<script src='/js/fornac.js'></script>";

	print "<div class='box-header' onclick='toggleBoxContent(this)'>RNA motif(s) Scan:</div>";
    print "<div class='box-content'>";
	
	open my $fh_tbl, '<', $tblout_file or die "Cannot open RNAmotif scan file: $!";
	while (my $line = <$fh_tbl>) {
		next if $line =~ /^#/;  # Skip comments
		chomp $line;
		my @columns = split(/\s+/, $line, 18);
		next unless @columns >= 16;

		my ($match, $family, $from, $to, $score, $e_value, $description) = ($columns[0], $columns[1], $columns[7], $columns[8], $columns[14], $columns[15], $columns[17]);

        push @rna_motif, $from, $to;

        ### Extract motif subsequence from full $seq
        my $start = $from - 1;
        my $length = $to - $from + 1;
        my $motif_seq = substr($SEQUENCECHECKED, $start, $length);

        # print ("Motif Seq: $motif_seq\n");

        ### Fold subsequence using RNAfold
        my $cmd = qq{echo "$motif_seq" | $VIENNARNAFOLDDIR/RNAfold --noPS 2>/dev/null};
        my $foldout = `$cmd`;  # run it!
        my @fold_lines = split("\n", $foldout);

        my $dot_bracket = '';
        $dot_bracket = $1 if $fold_lines[1] && $fold_lines[1] =~ /([().]+)\s+\([^)]+\)/;
        # print "Structure: $dot_bracket\n";
        # now save the results 
		my $family_link = "<a href=\"https://rfam.org/family/$family\" target=\"_blank\">$family</a>";

		# my $row = sprintf($format, $match, "$family_link     ", $from, $to, $score, $e_value, $description);
        # adding forna visual
       
        my $div_id = "rna_ss_$i";
        my $forna_html = "";

        $forna_html .= "<button onclick=\"toggleStructure$i()\">View Structure</button>\n";
        $forna_html .= "<div id='$div_id' style='width: 250px; height: 250px; display: none; margin: 10px; overflow: hidden;'></div>\n";
        $forna_html .= "<script>\n";
        $forna_html .= "  function toggleStructure$i() {\n";
        $forna_html .= "    var el = document.getElementById('$div_id');\n";
        $forna_html .= "    if (el.style.display === 'none') {\n";
        $forna_html .= "      el.style.display = 'block';\n";
        $forna_html .= "      var container = new fornac.FornaContainer('#$div_id', {\n";
        $forna_html .= "        animation: false,\n";
        $forna_html .= "        applyForce: false,\n";
        $forna_html .= "        labelInterval: 0,\n";
        $forna_html .= "        allowPanningAndZooming: true,\n";
        $forna_html .= "        structurePadding: 0,\n";
        $forna_html .= "        drawBackground: false,\n";
        $forna_html .= "        layout: 'naview'\n";
        $forna_html .= "      });\n";
        $forna_html .= "      var options = {\n";
        $forna_html .= "        structure: '$dot_bracket',\n";
        $forna_html .= "        sequence: '$motif_seq'\n";
        $forna_html .= "      };\n";
        $forna_html .= "      container.addRNA(options.structure, options);\n";
        $forna_html .= "    } else {\n";
        $forna_html .= "      el.innerHTML = '';\n";
        $forna_html .= "      el.style.display = 'none';\n";
        $forna_html .= "    }\n";
        $forna_html .= "  }\n";
        $forna_html .= "</script>\n";

        # push @results, $row . $forna_html;
        push @results, {
            match       => $match,
            family      => $family_link,
            from        => $from,
            to          => $to,
            score       => $score,
            evalue      => $e_value,
            description => $description,
            structure   => $forna_html,
        };

		$found = 1;
        $i++;
	}
	
	close $fh_tbl;

	
	if ($found) {
		# # Print header only if results exist
		# printf $format, "Match", "Family", "From", "To", "Score", "E-Value", "Description";
		# print "-" x 80, "\n";  # Simple separator

		# # Print stored results
		# print @results;
        print "<table class='table-result'>";
        print "<tr><th>Match</th><th>Family</th><th>From</th><th>To</th><th>Score</th><th>E-value</th><th>Description</th><th>Structure</th></tr>\n";

        foreach my $hit (@results) {
            print "<tr>";
            print "<td>$hit->{match}</td>";
            print "<td>$hit->{family}</td>";
            print "<td>$hit->{from}</td>";
            print "<td>$hit->{to}</td>";
            print "<td>$hit->{score}</td>";
            print "<td>$hit->{evalue}</td>";
            print "<td>$hit->{description}</td>";
            print "<td>$hit->{structure}</td>";
            print "</tr>";
        }
        
        print "</table>";

	} else {
        print "<div class='info-warning'>";
		print "No RNA motif recognized\n";
        print "</div>";
    }
    print "</div>";
}

sub CPC2 {

	my $cpc_input = "$TEMPDIR/$job.dna";
	my $cpc_output = "$TEMPDIR/$job.cpc2";

	my $RUN_CPC="python3 $CPC/CPC2.py --ORF -r -i $cpc_input -o $cpc_output";
	my $exit_code = system($RUN_CPC);
	if ($exit_code != 0) {
    	print "CPC2 execution failed with exit code: $exit_code\n";
	}

	print "<div class='box-header' onclick='toggleBoxContent(this)'>Coding potential Analysis:</div>";
    print "<div class='box-content'>";

	open(my $fh_cpc2, "<", "$cpc_output.txt") or die "Cannot open CPC2 result $cpc_output: $!";
	my @results;
	my $found = 0;

	# my $format = "%-10s %-18s %-15s %-10s %-10s %-10s %-15s %-10s %-10s\n";

    my ($orf_start, $peptide_length, $label);

	# Read file line by line
	while (my $line = <$fh_cpc2>) {
		next if $line =~ /^#/;  # Skip comment/header lines
		chomp $line;

		my @columns = split(/\t/, $line);  # Split by tab
		next unless @columns >= 9;  # Ensure enough columns exist

        $peptide_length = int($columns[2]);
        $orf_start      = int($columns[6]);
        $label          = $columns[9];
        $strand         = $columns[7];
		# Extract relevant fields
        my ($id, $transcript_length, $peptide_length, $fickett_score, $pI, $orf_integrity, $orf_start, $strand, $coding_probability, $label) = @columns[0..9];

		# Format numeric values with 2 decimal places
		$fickett_score = sprintf("%.2f", $fickett_score);
		$pI = sprintf("%.2f", $pI);
		$coding_probability = sprintf("%.6f", $coding_probability);  # Keep precision for probability

		# Store formatted row
		# push @results, sprintf($format, $id, $transcript_length, $peptide_length, $fickett_score, $pI, $orf_integrity, $orf_start, $coding_probability, $label);

        push @results, {
            id          => $id,
            tlength     => $transcript_length,
            plength     => $peptide_length,
            fscore      => $fickett_score,
            pI          => $pI,
            orf_int     => $orf_integrity,
            orf_start   => $orf_start,
            coding_prob => $coding_probability,
            label       => $label,

        };

        $found = 1;
    }
    close $fh_cpc2;

    if ($found) {
        # printf $format, "ID", "Transcript Length", "Peptide Length", "Fickett", "pI", "ORF", "ORF Start", "Coding Prob.", "Label";
        # print "-" x 110, "\n";
        # print @results;

        # making tbaular output
        print "<table class='table-result'>";
        print "<tr><th>ID</th><th>Transcript Length</th><th>Peptide Length</th><th>Fickett</th><th>pI</th><th>ORF</th><th>ORF Start</th><th>Coding Prob.</th><th>Label</th></tr>";

        foreach my $hit (@results) {
            print "<tr>";
            print "<td>$hit->{id}</td>";
            print "<td>$hit->{tlength}</td>";
            print "<td>$hit->{plength}</td>";
            print "<td>$hit->{fscore}</td>";
            print "<td>$hit->{pI}</td>";
            print "<td>$hit->{orf_int}</td>";
            print "<td>$hit->{orf_start}</td>";
            print "<td>$hit->{coding_prob}</td>";
            print "<td>$hit->{label}</td>";
            print "</tr>";
        }
        print "</table>"

    } else {
        print "<div class='info-warning'>";
        print "No result in CPC2 output\n";
        print "</div>";
    }
    
    if (defined $label && $label eq 'coding' && $orf_start > 0 && $peptide_length > 0) {

    my $orf_end = $orf_start + ($peptide_length * 3) - 1;

    # Extract forward-oriented sequence regardless of strand
    my $cds_nt = substr($SEQUENCECHECKED, $orf_start - 1, $peptide_length * 3);

    # Reverse complement if necessary
    if ($strand eq '-') {
        
        $cds_nt =~ tr/ACGTacgt/TGCAtgca/;
        $cds_nt = reverse($cds_nt);
    }

    # Translate
    my $cds_obj = Bio::Seq->new(-seq => $cds_nt, -alphabet => 'dna');
    my $aa_obj  = $cds_obj->translate(-stop_symbol => '');
    my $protein_seq = $aa_obj->seq;

    $protein_seq =~ s/\s//g;
    @predprot = split('', $protein_seq);

    $transcripts{"t1"} = {
        cds_start => $orf_start,
        cds_end   => $orf_end,
        strand    => $strand,
        source    => 'cpc2',
        protein   => $protein_seq,
    };
    } else {
        print "<div class='info-info'>";
        print "<i>Transcript predicted to be noncoding. Running structural region scan instead.</i>\n";
        print "</div>";
        print "<div class='box'>";
        &scan_structured_regions;
        print "</div>";
    }
    print "</div>";
    return \%transcripts;
    
}

sub microRNA {

		my $mirbase_output = "$TEMPDIR/$job.mirtbl";
		my $mirbase_out    = "$TEMPDIR/$job.mir";

		my $mirna_search = "$HMMER/nhmmer --rna --watson -Z 3.73 -E 1 --tblout $mirbase_output -o $mirbase_out $TEMPDIR/$job.seq $MIRBASE/hairpin.fa";

		system($mirna_search);

        print "<div class='box-header' onclick='toggleBoxContent(this)'>miRNA(s) scan:</div>";
        print "<div class='box-content'>";

		my @results;
        @mirna_loc = ();

		open my $fh_tbl, '<', $mirbase_output or die "Cannot open miRNA result: $!";
		while (my $line = <$fh_tbl>) {
			next if $line =~ /^#/;
			chomp $line;
			my @columns = split(/\s+/, $line);
			next unless @columns >= 17;
			my $desc_full = join(" ", @columns[15..$#columns]);

			my ($match, $from, $to, $e_value, $score) =
    			($columns[0], $columns[7], $columns[8], $columns[12], $columns[13]);

            push @mirna_loc, $from, $to; 

			# Extract accession and description
			my ($accession, $desc_text) = $desc_full =~ /(MI\w+\d+)\s+(.*)/;
			unless ($accession) {
				$accession = '-';
				$desc_text = $desc_full;
			}

			push @results, {
				match       => $match,
				from        => $from,
				to          => $to,
				score       => $score + 0,    # force numeric
				e_value     => $e_value + 0,  # force numeric
				accession   => $accession,
				description => $desc_text
			};
		}
		close $fh_tbl;

		my $total = scalar @results;

		if ($total) {
			# print "<b>miRNA search:</b><br>\n";
			# printf $format, "Match", "From", "To", "E-Value", "Score", "Accession", "Description";
			# print "-" x 120, "\n";
            #making tabular
            print "<table class='table-result'>";
            print "<tr><th>Match</th><th>From</th><th>To</th><th>E-Value</th><th>Score</th><th>Accession</th><th>Description</th></tr>";

			# Prioritize human miRNA
			my (@human_hits, @other_hits);
			foreach my $hit (@results) {
				if ($hit->{match} =~ /^hsa-/i || $hit->{description} =~ /Homo sapiens/i) {
					push @human_hits, $hit;
				} else {
					push @other_hits, $hit;
				}
			}

			# Sort both groups by E-value ascending
			@human_hits = sort { $a->{e_value} <=> $b->{e_value} } @human_hits;
			@other_hits = sort { $a->{e_value} <=> $b->{e_value} } @other_hits;

			my @top_hits = (@human_hits, @other_hits);
			@top_hits = @top_hits[0..2] if @top_hits > 3;

			foreach my $hit (@top_hits) {
				my $link = "<a href=\"https://www.mirbase.org/hairpin/$hit->{accession}\" target=\"_blank\">$hit->{accession}</a>";
				# printf $format, $hit->{match}, $hit->{from}, $hit->{to}, $hit->{e_value}, $hit->{score}, $link, $hit->{description};
                print "<tr>";
                print "<td>$hit->{match}</td>";
                print "<td>$hit->{from}</td>";
                print "<td>$hit->{to}</td>";
                print "<td>$hit->{e_value}</td>";
                print "<td>$hit->{score}</td>";
                print "<td>$link</td>";
                print "<td>$hit->{description}</td>";
                print "</tr>";           
            }
            print "</table>";
			print "<div class='info-warning'>Total microRNA hits found:</b> $total</div>";
		} else {
            print "<div class='info-warning'>";
			print "No regions matching a mircroRNA was found.\n";
            print "</div>";
        }

		print "</div>";


}

sub riboswitch {
    my $ribosw_tblout = "$TEMPDIR/$job.ribo.tblout";  # Table format output
	my $ribosw_out = "$TEMPDIR/$job.ribo";     # Full verbose output

	my $cmd = "$CMSCAN/cmscan -E 1e-5 --tblout $ribosw_tblout -o $ribosw_out $RIBOSWDB $TEMPDIR/$job.seq > /dev/null 2>&1";
	system($cmd);

    @ribosw =();
    my @results;
    my $found = 0;
    my $i = 0;

	print "<div class='box-header' onclick='toggleBoxContent(this)'>Riboswitch Scan:</div>";
    print "<div class='box-content'>";
	
	open my $fh_tbl, '<', $ribosw_tblout or die "Cannot open Riboswitch out file: $!";
	while (my $line = <$fh_tbl>) {
		next if $line =~ /^#/;  # Skip comments
		chomp $line;
		my @columns = split(/\s+/, $line, 18);
		next unless @columns >= 16;

		my ($match, $family, $from, $to, $score, $e_value, $description) = ($columns[0], $columns[1], $columns[7], $columns[8], $columns[14], $columns[15], $columns[17]);

        push @ribosw, $from, $to;

       	my $family_link = "<a href=\"https://rfam.org/family/$family\" target=\"_blank\">$family</a>";

        ### Extract motif subsequence from full $seq
        my $start = $from - 1;
        my $length = $to - $from + 1;
        my $ribosw_seq = substr($SEQUENCECHECKED, $start, $length);

        #adding stem check
        my @check_stem = checkstemsonly($ribosw_seq, 1);

        my $stemsinfo;
        if ($check_stem[0] != $check_stem[1]) {
            $stemsinfo = "$check_stem[0]-$check_stem[1]";
        } else {
            $stemsinfo = "$check_stem[0]";
        }

        my $energy = ($check_stem[2] != 1) ? $check_stem[2] : "";

        # adding biological context
        my $bio_context = 'Metabolite sensing; gene regulation';  # default
        if ($match =~ /TPP/i)      { $bio_context = 'Vitamin B1 metabolism; thiamine biosynthesis'; }
        elsif ($match =~ /FMN/i)   { $bio_context = 'Riboflavin (B2) biosynthesis; flavin transport'; }
        elsif ($match =~ /SAM/i)   { $bio_context = 'Methionine metabolism; methylation reactions'; }
        elsif ($match =~ /AdoCbl/i){ $bio_context = 'Vitamin B12 metabolism; cobalamin transport'; }
        elsif ($match =~ /Cobalamin/i){ $bio_context = 'Vitamin B12 metabolism; cobalamin transport'; }
        elsif ($match =~ /Glycine/i){ $bio_context = 'Amino acid metabolism; glycine cleavage system'; }
        elsif ($match =~ /Lysine/i) { $bio_context = 'Amino acid biosynthesis; lysine transport'; }
        elsif ($match =~ /glmS/i)   { $bio_context = 'Cell wall synthesis; self-cleaving riboswitch'; }
        elsif ($match =~ /purine/i) { $bio_context = 'Purine metabolism; adenine/guanine transport'; }
        elsif ($match =~ /preQ1/i)  { $bio_context = 'Queuosine precursor; tRNA modification'; }
        elsif ($match =~ /ykoK/i)   { $bio_context = 'Magnesium sensing; Mg2+ homeostasis'; }
        elsif ($match =~ /yybP/i)   { $bio_context = 'Manganese homeostasis; metal ion transport'; }
        elsif ($match =~ /c-di-GMP/i){ $bio_context = 'Cyclic diguanylate signaling; biofilm regulation'; }

        # push @results, $row . $forna_html;
        push @results, {
            match       => $match,
            family      => $family_link,
            from        => $from,
            to          => $to,
            score       => $score,
            evalue      => $e_value,
            description => $description,
            stems       => $stemsinfo,
            energy      => $energy,
            context     => $bio_context,
        };

		$found = 1;
        $i++;
	}
	
	close $fh_tbl;

	
	if ($found) {
		# # Print header only if results exist
		# printf $format, "Match", "Family", "From", "To", "Score", "E-Value", "Description";
		# print "-" x 80, "\n";  # Simple separator

		# # Print stored results
		# print @results;
        print "<table class='table-result'>";
        print "<tr><th>Match</th><th>Family</th><th>From</th><th>To</th><th>Score</th><th>E-value</th><th>Stems</th><th>Energy</th><th>Description</th><th>Function</th></tr>\n";

        foreach my $hit (@results) {
            print "<tr>";
            print "<td>$hit->{match}</td>";
            print "<td>$hit->{family}</td>";
            print "<td>$hit->{from}</td>";
            print "<td>$hit->{to}</td>";
            print "<td>$hit->{score}</td>";
            print "<td>$hit->{evalue}</td>";
            print "<td>$hit->{stems}</td>";
            print "<td>$hit->{energy}</td>";
            print "<td>$hit->{description}</td>";
            print "<td>$hit->{context}</td>";
            print "</tr>";
        }
        
        print "</table>";

	} else {
        print "<div class='info-warning'>";
		print "No Riboswitch sequence recognized\n";
        print "</div>";
    }
    print "</div>";
}

# microRNA target scan using miRanda
# slow, takes over 2 minutes to scan
# trying to implement only 3' UTR scan so it less intesive but still biologically relevant with fallback.
sub miRNAtarget {
    my ($transcripts_ref) = @_;
    my $mirna_db     = "$MIRBASE/mature.fa";
    my $raw_out      = "$TEMPDIR/$job.miranda.out";
    my $miranda_out  = "$TEMPDIR/$job.miranda.tsv";
    my $utr_fasta    = "$TEMPDIR/$job.utr3.fa";

    my $seq = $SEQUENCECHECKED;
    $seq =~ tr/uU/tT/;


    print "<div class='box-header' onclick='toggleBoxContent(this)'>miRNA target prediction(s):</div>\n";
    print "<div class='box-content'>";

    # 1. Create UTR3 FASTA file (or fallback)
    my $found_valid = 0;
    open my $fa_out, '>', $utr_fasta or die "Cannot write $utr_fasta: $!";

    if (%$transcripts_ref) {
        foreach my $tid (keys %$transcripts_ref) {
            next unless $transcripts_ref->{$tid}{utr3} && @{ $transcripts_ref->{$tid}{utr3} };  # Ensure utr3 exists and has elements

            my @coords = @{ $transcripts_ref->{$tid}{utr3} };
            for (my $i = 0; $i < @coords; $i += 2) {
                my ($start, $end) = @coords[$i, $i+1];
                next unless defined $start && defined $end && $end >= $start;

                next if $start < 1 || $end > length($seq) || $end < $start;
                my $utr_seq = substr($seq, $start - 1, $end - $start + 1);
                next unless $utr_seq =~ /[ACGT]/i;

                $found_valid++;
                print $fa_out ">$tid|utr3_${start}_${end}\n$utr_seq\n";
            }
        }
    }


    print "<div class='info-info'>Extracted $found_valid valid 3' UTR regions for target prediction.</i></div>";

    # no valid UTRs 
    if (!$found_valid) {
        print "<div class='info-warning'>";
        print "<i>No valid 3' UTRs found. Scanning full sequence instead.</i></div>";
        print $fa_out ">full_sequence\n$seq\n";
    }

    close $fa_out;

    my $cmd = "python3 $MIRANDA/miranda_wrapper.py --parsed_out $miranda_out --miranda_bin $MIRANDA/miranda --tmpdir $TEMPDIR --job_id $job  $mirna_db $utr_fasta $raw_out";
    my $exit_code = system($cmd);
    if ($exit_code != 0) {
        print "miRanda wrapper execution failed with exit code: $exit_code\n";
        return ();
    }

    open(my $fh, "<", $miranda_out) or die "Can't open $miranda_out: $!";
    my $header = <$fh>;  # Skip header
    my @lines = <$fh>;
    chomp @lines;
    close($fh);

    my @sorted = sort {
        (split /\t/, $a)[4] <=> (split /\t/, $b)[4] ||
        (split /\t/, $a)[5] <=> (split /\t/, $b)[5] ||
        (split /\t/, $a)[3] <=> (split /\t/, $b)[3]
    } @lines;

    my @filtered;
    my @current_group;

    foreach my $line (@sorted) {
        my ($query, $mirna, $score, $energy, $start, $end) = (split /\t/, $line)[0,1,2,3,4,5];

        if (!@current_group) {
            push @current_group, $line;
            next;
        }

        my ($prev_start, $prev_end) = (split /\t/, $current_group[-1])[4,5];

        if ($start <= $prev_end) {
            push @current_group, $line;
        } else {
            my $best = (sort { (split /\t/, $a)[3] <=> (split /\t/, $b)[3] } @current_group)[0];
            push @filtered, $best;
            @current_group = ($line);
        }
    }

    if (@current_group) {
        my $best = (sort { (split /\t/, $a)[3] <=> (split /\t/, $b)[3] } @current_group)[0];
        push @filtered, $best;
    }


    print "<table class='table-result'>";
    print "<tr><th>miRNA</th><th>From</th><th>To</th><th>Energy</th><th>Score</th></tr>";

    my @regions;
    foreach my $line (@filtered) {
        my (undef, $mirna, $score, $energy, $start, $end) = split /\t/, $line;
        print "<tr><td>$mirna</td><td>$start</td><td>$end</td><td>$energy</td><td>$score</td></tr>";
        push @regions, [$start, $end];
    }
    print "</table>";
    print "</div>";
    return @regions;
}

########################
# augustus replacing the old genscan; need to check the whole sub for errors
# Flag inferred UTRs as “low confidence” when:
# ORF is very short
# ORF start is <15 nt from sequence start (no room for 5′ UTR)
# Sequence ends shortly after ORF (truncated transcript)
# Provide scoring or confidence per UTR:
# E.g., "Predicted 5′ UTR: 37 bp (contains weak Shine-Dalgarno motif, ΔG = -3.2 kcal/mol)"
# Refactored AUGUSTUS + UTR + PolyA logic to handle multiple transcripts

sub AUGUSTUS {

    print "<div class='box-header' onclick='toggleBoxContent(this)'>Gene Prediction Analysis:</div>";
    print "<div class='box-content'>";
    my ($species) = @_;
    $species ||= "human";
    my $utr_flag = ($species =~ /^(human|fly|zebrafish)$/i) ? "--UTR=on" : "";

    my $output_gff = "$TEMPDIR/$job.augustus";
    my $input_dna  = "$TEMPDIR/$job.dna";

    my $augustus_cmd = "$AUGUSTUS --softmasking=0 --protein=on $utr_flag --species=$species $input_dna > $output_gff 2>&1";
    system($augustus_cmd) == 0 or die "AUGUSTUS run failed: $!";

    open(my $GFF, '<', $output_gff) or die "Can't open AUGUSTUS output: $!";

    my %transcripts;
    my $current_tid;
    my @protein_lines;
    my $capturing_protein = 0;

    while (my $line = <$GFF>) {
        chomp $line;

        if ($line =~ /^# protein sequence = \[(.*)$/) {
            $capturing_protein = 1;
            push @protein_lines, $1;
            next;
        }
        if ($capturing_protein) {
            $line =~ s/^#\s?//;
            if ($line =~ /^(.*)\]$/) {
                push @protein_lines, $1;
                $capturing_protein = 0;
                $transcripts{$current_tid}{protein} = join('', @protein_lines);
                @protein_lines = ();
                next;
            }
            push @protein_lines, $line;
            next;
        }

        next if $line =~ /^#/;
        my @fields = split("\t", $line);
        next unless @fields >= 9;

        my ($type, $start, $end, $strand, $attr) = @fields[2,3,4,6,8];
        my ($tid) = $attr =~ /transcript_id \"(.*?)\"/;
        $current_tid = $tid if defined $tid;
        next unless defined $tid;

        $transcripts{$tid}{strand} = $strand if $type eq 'gene';
        $transcripts{$tid}{tss} = $start if $type eq 'tss';
        $transcripts{$tid}{tts} = $end   if $type eq 'tts';
        if ($type eq 'CDS') {
            push @{ $transcripts{$tid}{cds} }, [$start, $end];
            $transcripts{$tid}{cds_start} //= $start;
            $transcripts{$tid}{cds_end} = $end;
        }
        if ($type eq 'exon') {
            push @{ $transcripts{$tid}{exons} }, [$start, $end];
        }
        $transcripts{$tid}{source} = 'augustus';
    }
    close $GFF;

    my $i = 1;
    my $j = 1;
    if (keys %transcripts == 0) {
        print "<div class='info-warning'>";
        print "No gene predictions were found.</div>";
    } else {
    foreach my $tid (sort keys %transcripts) {
        my $model = $transcripts{$tid};
        
        my ($gene_id) = $tid =~ /^(g\d+)\./;
        print "<h3>Transcript $i (ID: $tid";
        print ", gene: $gene_id" if $gene_id;
        print ")</h3><br>";
        $i++;

        print "<table class='table-result'>";
        print "<tr><th>Exon</th><th>Start</th><th>End</th><th>Type</th></tr>";

        foreach my $exon (@{ $model->{exons} }) {
            my ($start, $end) = @$exon;
            my @segments = split_exon_by_cds($start, $end, @{ $model->{cds} });
            foreach my $seg (@segments) {
                my ($seg_start, $seg_end, $type) = @$seg;
                print "<tr><td>Exon $j</td><td>$seg_start</td><td>$seg_end</td><td>$type</td></tr>";
                $j++;
            }
        }
        print "</table>";
        $model->{protein} =~ s/\s//g if exists $model->{protein};

        print "<br>";
        }
    }

    @predprot = ();              
    foreach my $tid (sort keys %transcripts) {
        if (my $protein_seq = $transcripts{$tid}{protein}) {
            $protein_seq =~ s/\s//g;
            my @residues = split('', $protein_seq);
            push @predprot, @residues;
            $predprotforAnDom .= $protein_seq;
        }
    }
    print "Parsed type=$type, tid=$tid, start=$start, end=$end<br>" if defined $tid;
    print "</div>";

    return \%transcripts;
}

sub split_exon_by_cds {
    my ($exon_start, $exon_end, @cds_ranges) = @_;
    my @segments;
    my $cursor = $exon_start;

    foreach my $cds (@cds_ranges) {
        my ($cds_start, $cds_end) = @$cds;
        next if $cds_end < $exon_start || $cds_start > $exon_end;

        if ($cds_start > $cursor) {
            push @segments, [$cursor, $cds_start - 1, 'noncoding'];
        }

        my $coding_start = ($cds_start > $cursor) ? $cds_start : $cursor;
        my $coding_end = ($cds_end < $exon_end) ? $cds_end : $exon_end;
        push @segments, [$coding_start, $coding_end, 'coding'];
        $cursor = $coding_end + 1;
    }

    push @segments, [$cursor, $exon_end, 'noncoding'] if $cursor <= $exon_end;
    return @segments;
}

sub predict_utrs {
    my (%args) = @_;

    my $seq        = $args{seq};
    my $cds_start  = $args{cds_start};
    my $cds_end    = $args{cds_end};
    my $strand     = $args{strand} || '+';
    my $source     = $args{source} || 'unknown';
    my $tss_pos    = $args{tss};
    my $tts_pos    = $args{tts};

    $seq =~ tr/uU/tT/;
    my $seq_length = length($seq);

    my @new5primeutr = ();
    my @new3primeutr = ();
    my @utrprintout  = ();
    my @utr          = ();
    my @results      = ();

    print "<div class='box-header' onclick='toggleBoxContent(this)'>UTR(s) Prediction:</div>";
    print "<div class='box-content'>";

    print "<div class='info-info'>UTR prediction source: $source</i></div>";

    if ($source eq 'augustus' && (defined $tss_pos || defined $tts_pos)) {
        print "<i>Using AUGUSTUS-predicted TSS/TTS boundaries.</i><br>";

        if (defined $tss_pos && $tss_pos < $cds_start) {
            push @new5primeutr, $tss_pos, $cds_start - 1 if $tss_pos < $cds_start;
            push @utrprintout, 5, $tss_pos, $cds_start - 1;
        }
        if (defined $tts_pos && $cds_end < $tts_pos) {
            push @new3primeutr, $cds_end + 1, $tts_pos if $cds_end + 1 <= $tts_pos;
            push @utrprintout, 3, $cds_end + 1, $tts_pos;
        }
    } else {
        print "<i>Inferring UTRs from CDS boundaries...</i><br>";

        if ($strand eq '+') {
            if ($cds_start > 1) {
                my $utr_start = 1;
                my $utr_end   = $cds_start - 1;
                push @new5primeutr, $utr_start, $utr_end if $utr_start <= $utr_end;
            }

            if ($cds_end < $seq_length) {
                my $utr_start = $cds_end + 1;
                my $utr_end   = $seq_length;
                push @new3primeutr, $utr_start, $utr_end if $utr_start <= $utr_end;
            }
        } else {
            if ($cds_end < $seq_length) {
                my $utr_start = $cds_end + 1;
                my $utr_end   = $seq_length;
                push @new5primeutr, $utr_start, $utr_end if $utr_start <= $utr_end;
            }

            if ($cds_start > 1) {
                my $utr_start = 1;
                my $utr_end   = $cds_start - 1;
                push @new3primeutr, $utr_start, $utr_end if $utr_start <= $utr_end;
            }
        }

        @utrprintout = (
            5, @new5primeutr,
            3, @new3primeutr
        );
    }

    @utr = sort { $a <=> $b } (@new5primeutr, @new3primeutr);

    # print "<b>UTR:</b>           start  -   end   -  stems - energy<br>";
    for (my $i = 0; $i < @utrprintout; $i += 3) {
        my ($type, $start, $end) = @utrprintout[$i, $i+1, $i+2];

        next if $end < $start || $start < 1 || $end > $seq_length;

        my $utr_seq = substr($seq, $start - 1, $end - $start + 1);
        next unless length($utr_seq) > 0;

        my @returnout = checkstemsonly($utr_seq, 1);

        # printf(" %d'            %-6d - %6d", $type, $start, $end);
        # print "       $returnout[0]-$returnout[1]" if $returnout[0] != $returnout[1];
        # print "       $returnout[0]" if $returnout[0] == $returnout[1];
        # print "     $returnout[2]<br>" if $returnout[2] != 1;

        # collecting info for tabular output

        my $stems_info;
        if ($returnout[0] != $returnout[1]) {
            $stems_info = "$returnout[0]-$returnout[1]";
        } else {
            $stems_info = "$returnout[0]";
        }

        my $energy = ($returnout[2] != 1) ? $returnout[2] : "";
        my @motifs = ();

        if ($type == 5) {
            if ($utr_seq =~ /AGGAGG/i) {
                my $pos = $-[0] + $start;
                push @motifs, "SD motif AGGAGG at $pos<br>";
            }
            if ($utr_seq =~ /gcc[AG]ccATGG/i) {
                my $pos = $-[0] + $start;
                my $matched_seq = $&;
                push @motifs, "Kozak motif $matched_seq at $pos<br>";
            }
        }

        if ($type == 3) {
            foreach my $motif (qw(AATAAA ATTAAA TATAAA AAGAAA AGTAAA AATATA)) {
                if ($utr_seq =~ /$motif/i) {
                    my $pos = $-[0] + $start;
                    push @motifs, "PolyA signal $motif at $pos";
                    my $signal_end = $pos + length($motif) - 1;
                    push @polyasignal, $pos, $signal_end;
                    last;
                }
            }
            if ($utr_seq =~ /A{10,}/i) {
                my $tail_pos = $-[0] + $start;
                push @motifs, "PolyA tail near $tail_pos<br>";
                my $tail_end = $+[0] - 1 + $start;
                push @polyatail, $tail_pos, $tail_end;
            }
        }

        push @results, {
            type        => "${type}'",
            start       => $start,
            end         => $end,
            length      => $end - $start + 1,
            stems       => $stems_info,
            energy      => $energy,
            motifs      => join(";<br>", @motifs) || "None"
        };

    }

     if (@results) {
        print "<table class='table-result'>";
        print "<tr><th>UTR Type</th><th>Start</th><th>End</th><th>Length</th><th>Stems</th><th>Energy</th><th>Motifs</th></tr>";

        foreach my $hit (@results) {
            print "<tr>";
            print "<td>$hit->{type}</td>";
            print "<td>$hit->{start}</td>";
            print "<td>$hit->{end}</td>";
            print "<td>$hit->{length}</td>";
            print "<td>$hit->{stems}</td>";
            print "<td>$hit->{energy}</td>";
            print "<td>$hit->{motifs}</td>";
            print "</tr>";
        }
        print "</table>";
    } else {
        print "<div class='info-warning'>";
        print "No valid UTRs inferred.</div>" unless @utr;
    }
    print "</div>";
    return (\@new5primeutr, \@new3primeutr, \@utrprintout, \@utr, \@polyasignal, \@polyatail);
}

sub refineUTRwithPolyA {
    my ($sequence, $cds_end, $strand, $seq_length) = @_;
    $sequence =~ tr/uU/tT/;

    my @new3primeutr = ();
    my @polyasignal  = ();
    my @polyatail    = ();
    my @utrprintout  = ();
    my @utr          = ();

    my @motifs = qw(AATAAA ATTAAA TATAAA AAGAAA AATATA AGTAAA);
    my $window_start = $cds_end + 1;
    my $window_end   = $seq_length;
    my $utr_seq      = substr($sequence, $window_start - 1, $window_end - $window_start + 1);

    my ($signal_pos, $signal_motif, $tail_pos) = (-1, '', -1);

    foreach my $motif (@motifs) {
        if ($utr_seq =~ /$motif/i) {
            $signal_pos = $-[0] + $window_start;
            $signal_motif = $motif;
            my $signal_end = $signal_pos + length($motif) - 1;
            push @polyasignal, $signal_pos, $signal_end;
            print "<div class='info-info'>";
            print "PolyA signal ($motif) detected at $signal_pos<br>";
            print "</div>";
            last;
        }
    }

    while ($sequence =~ /A{10,}/g) {
        my $match_start = $-[0];
        my $match_end   = $+[0] - 1;

        if ($match_start >= $cds_end && $match_start - $cds_end < 200) {
            $tail_pos  = $match_start + 1;
            my $tail_end = $match_end + 1;
            push @polyatail, $tail_pos, $tail_end;
            print "<div class='info-info'>";
            print "PolyA tail detected near $tail_pos<br>";
            print "</div>";
            last;
        }
    }

    if ($signal_pos > 0 || $tail_pos > 0) {
        my $utr_start = $cds_end + 1;
        my $utr_end   = $tail_pos > 0 ? $tail_pos : ($signal_pos + 20);

        push @new3primeutr, $utr_start, $utr_end;
        push @utrprintout, 3, $utr_start, $utr_end;
        push @utr, $utr_start, $utr_end;

        print "<div class='info-info'>";
        print "Inferred 3' UTR based on polyA: $utr_start - $utr_end<br>";
        print "</div>";
    } else {
        print "<div class='info-warning'>";
        print "No strong polyA signal/tail detected in 3' region.<br>";
        print "</div>";
    }

    return (\@new3primeutr, \@polyasignal, \@polyatail);
}

sub scan_rbs {
    my ($transcripts_ref) = @_;
    my $seq_file = "$TEMPDIR/$job.dna.fa";
    my $coords_file = "$TEMPDIR/$job.rbs.coords";
    my $output_file = "$TEMPDIR/$job.rbs.out";

    # 1. Write CDS coordinates for all transcripts
    open my $coord_fh, '>', $coords_file or die "Can't write RBS coords: $!";
    foreach my $tid (sort keys %$transcripts_ref) {
        my $model = $transcripts_ref->{$tid};
        next unless defined $model->{cds_start} and defined $model->{cds_end};
        print $coord_fh "$tid\t$model->{cds_start}\t$model->{cds_end}\n";
    }
    close $coord_fh;
    my $utr5_len = 0;
    if (ref $utr5_ref eq 'ARRAY' && @$utr5_ref == 2) {
        $utr5_len = abs($utr5_ref->[1] - $utr5_ref->[0]) + 1;
        $utr5_len = 20 if $utr5_len < 20;  # fallback
    }

    # 2. Run the RBS finder
    my $cmd = "perl $RBSFINDER/rbs_finder.pl $seq_file $coords_file $output_file $utr5_len AGGAGG";
    system($cmd) == 0 or die "RBS finder failed: $?";

    # 3. Parse summary output minimally
    open my $rbs_out, '<', $output_file or die "Can't open RBS output: $!";
    my ($with_rbs, $without_rbs, $total);

    while (my $line = <$rbs_out>) {
        if ($line =~ /have RBS.*?= (\d+)/) {
            $with_rbs += $1;
        }
        if ($line =~ /have no RBS= (\d+)/) {
            $without_rbs = $1;
        }
        if ($line =~ /Total.*?= (\d+)/) {
            $total = $1;
        }
    }
    close $rbs_out;

    # 4. Display minimal summary
    print "<b>RBS Scan:</b> ";
    if ($with_rbs) {
        my $percent = sprintf("%.1f", 100 * $with_rbs / $total);
        print "$with_rbs of $total transcript(s) have a Shine-Dalgarno motif ($percent%)<br>";
    } else {
        print "No Shine-Dalgarno motifs detected in 5' UTR regions.<br>";
    }
}

sub scan_structured_regions {
    my $count = 0;
    my $len = length $SEQUENCECHECKED;
    my $printout = 0;
    my @high_struct_regions;

    while ($count <= $len) {
        my $query_len = ($count + 150 < $len) ? 150 : $len - $count;
        my $query = substr($SEQUENCECHECKED, $count, $query_len);
        my @answer = &checkstemsonly($query, 1);  # (stem_start, stem_end, energy)

        if ($answer[0] >= 3 && $answer[1] >= 3) {
            push @high_struct_regions, [
                $count + 1, $count + $query_len,
                ($answer[0] == $answer[1]) ? "$answer[0]" : "$answer[0]-$answer[1]", $answer[2]
            ];
            $printout = 1;
        }

        last if $count + 150 >= $len;
        $count += 150;
    }
    print "<div class='box-header' onclick='toggleBoxContent(this)'>Structured region(s) scan:</div>";
    print "<div class='box-content'>";
    if ($printout) {
        my $i = 1;
        # print "<b>Structured regions detected:</b><br>";
        print "<table class='table-result'>";
        print "<tr><th>Region</th><th>From</th><th>To</th><th>Stems</th><th>Energy</th></tr>";
        # printf "%-12s %-6s %-6s %-10s %-10s<br>", "Region", "From", "To", "Stems", "Energy";
        # print "-" x 60 . "<br>";
        foreach my $r (@high_struct_regions) {
            print "<tr>";
            print "<td>$i</td>";
            print "<td>$r->[0]</td>";
            print "<td>$r->[1]</td>";
            print "<td>$r->[2]</td>";
            print "<td>$r->[3]</td>";
            $i++;
        }
        print "</table>";
        print "<div class='info-info'>Highly structured regions found. Consider tRNA, rRNA, or ncRNA elements.</div>";
    } else {
        print "<div class='info-warning'>";
        print "No regions with significant RNA structure detected.";
        print "</div>";
    }
    print "<br>";
    print "</div>";
}

sub normalize_transcript_features {
    my ($transcripts_ref, $exons_ref, $utr_ref, $polyasignal_ref, $polyatail_ref, $mirnatarget_ref, $structured_ref) = @_;

    @$exons_ref         = ();
    @$utr_ref           = ();
    @$polyasignal_ref   = ();
    @$polyatail_ref     = (); 
    @$mirnatarget_ref   = ();
    @$structured_ref    = ();

    foreach my $tid (keys %$transcripts_ref) {
        my $model = $transcripts_ref->{$tid};
        if (defined $model->{cds_start} && defined $model->{cds_end}) {
            push @$exons_ref, $model->{cds_start}, $model->{cds_end};
        }

        push @$utr_ref, @{ $model->{utr5} }               if ref $model->{utr5} eq 'ARRAY';
        push @$utr_ref, @{ $model->{utr3} }               if ref $model->{utr3} eq 'ARRAY';
        push @$polyasignal_ref, @{ $model->{polyAcoords} }     if ref $model->{polyAcoords} eq 'ARRAY';
        push @$polyatail_ref,   @{ $model->{polyAtailcoords} } if ref $model->{polyAtailcoords} eq 'ARRAY';
        push @$mirnatarget_ref, @{ $model->{mirnatargets} }     if ref $model->{mirnatargets} eq 'ARRAY';
        push @$structured_ref,  @{ $model->{structured_regions} } if ref $model->{structured_regions} eq 'ARRAY';

        

    }
    # print "DEBUG (normalize): final polyasignal = @$polyasignal_ref\n";
    # print "DEBUG (normalize): final polyatail   = @$polyatail_ref\n";
}

sub createfoldingpicture {
    my $seq_file       = "$TEMPDIR/$job.seq";
    my $foldout_file   = "$TEMPDIR/$job.foldout";

    # Print HTML content
    print "<div class='box-header' onclick='toggleBoxContent(this)'>RNA Structure Analysis:</div>";
    print "<div class='box-content-structure'>";

        

	# if ($SEQUENCELENGTH <= $MAXFOLDINGLEN && $SEQUENCELENGTH > $MAXFORNALENGTH) {
    #     my $svg_file = "$TEMPDIR/${SEQNAMECHECKED}_ss.svg";

    #     # Make sure RNAplot output is ready
    #     system("$VIENNARNAFOLDDIR/RNAplot --infile=$TEMPDIR/$job.foldout -f svg --filename-full");

    #     # Read SVG file content
    #     open(my $svgfh, '<', $svg_file) or die "Cannot open SVG file: $!";
    #     my $svg_content = do { local $/; <$svgfh> };
    #     close($svgfh);

    #     # Add ID to <svg> tag if not present
    #     $svg_content =~ s/<svg /<svg id="rna_ss" width="650" height="650" /;

    #     # Output
    #     print "<h3>RNA Structure Visualization:</h3>\n";
    #     print "$svg_content\n";
    #     print "<p style='font-size: 0.9em; color: gray;'> Drag to pan, scroll to zoom</p>\n";

    #     # Add svg-pan-zoom script
    #     print "<script src='/js/svg-pan-zoom.min.js'></script>\n";
    #     print "<script>\n";
    #     print "  svgPanZoom('#rna_ss', {\n";
    #     print "    zoomEnabled: true,\n";
    #     print "    controlIconsEnabled: true,\n";
    #     print "    fit: true,\n";
    #     print "    center: true\n";
    #     print "  });\n";
    #     print "</script>\n";
    # }



    if ($SEQUENCELENGTH <= $MAXFOLDINGLEN) {
        
        my $svg_file = "$TEMPDIR/${SEQNAMECHECKED}_ss.svg";
        my $ps_url  = "/tmp/jobs/job_$job/${SEQNAMECHECKED}_ss.ps";
        my $svg_url = "/tmp/jobs/job_$job/${SEQNAMECHECKED}_ss.svg";

        # # Make sure RNAplot output is ready
        # system("$VIENNARNAFOLDDIR/RNAplot --infile=$TEMPDIR/$job.foldout -f svg --filename-full");

        # Read sequence and structure from RNAfold output
        open(my $fh, '<', $foldout_file) or die "Cannot open foldout file: $!";
        my $header = <$fh>;  # Skip header line (e.g. >job123)
        my $sequence = <$fh>;
        chomp($sequence);
        my $structure_line = <$fh>;
        chomp($structure_line);
        $structure_line =~ /^([().]+)\s+/;
        my $structure = $1;
        close($fh);
          
        print "<div id='rna_ss'></div>";
        print "<p style='font-size: 0.9em; color: gray; text-align: center;'> Drag to pan, scroll to zoom</p>";
        print "<b>Download Folding As: </b>\n";
        print "<a href='$svg_url' target='_blank'><button>SVG File</button></a>";
        print "<a href='$ps_url' target='_blank'><button>PS File</button></a>\n";

        # Include required scripts
        print "<link rel='stylesheet' href='/css/fornac.css'>\n";
        print "<script src='/js/d3.v3.min.js'></script>\n";
        print "<script src='/js/fornac.js'></script>";

        my $color_text = "";

        # Process UTR 
        
        for (my $j = 0; $j < @utr; $j += 2) {
            my ($from, $to) = @utr[$j, $j + 1];
            for (my $i = $from; $i <= $to; $i++) {
                $color_text .= "$i:lightblue ";
            }
        }

        # Process exons 
        for (my $j = 0; $j < @exons; $j += 2) {
            my ($from, $to) = @exons[$j, $j + 1];
            for (my $i = $from; $i <= $to; $i++) {
                $color_text .= "$i:green ";
            }
        }

        #polya signals
        for (my $j = 0; $j < @polyasignal; $j += 2) {
            my ($from, $to) = @polyasignal[$j, $j + 1];

            ($from, $to) = ($to, $from) if $from > $to;

            for (my $i = $from; $i <= $to; $i++) {
                $color_text .= "$i:lime ";
            }
        }

        for (my $j = 0; $j < @polyatail; $j += 2) {
            my ($from, $to) = @polyatail[$j, $j + 1];

            ($from, $to) = ($to, $from) if $from > $to;

            for (my $i = $from; $i <= $to; $i++) {
                $color_text .= "$i:lime ";
            }
        }


        # process motif
        for (my $j = 0; $j < @rna_motif; $j += 2) {
            my ($from, $to) = @rna_motif[$j, $j + 1];
            for (my $i = $from; $i <= $to; $i++) {
                $color_text .= "$i:red ";
            }
        }

        # process mirna locations
        for (my $j = 0; $j < @mirna_loc; $j += 2) {
            my ($from, $to) = @mirna_loc[$j, $j + 1];

            # Ensure correct order
            ($from, $to) = ($to, $from) if $from > $to;

            for (my $i = $from; $i <= $to; $i++) {
                $color_text .= "$i:yellow ";
            }
        }

        #trna
        for (my $j = 0; $j < @trna_loc; $j += 2) {
            my ($from, $to) = @trna_loc[$j, $j + 1];

            ($from, $to) = ($to, $from) if $from > $to;

            for (my $i = $from; $i <= $to; $i++) {
                $color_text .= "$i:orange ";
            }
        }

        # smsites
        for (my $j = 0; $j < @smsite; $j += 2) {
            my ($from, $to) = @smsite[$j, $j + 1];

            ($from, $to) = ($to, $from) if $from > $to;

            for (my $i = $from; $i <= $to; $i++) {
                $color_text .= "$i:magenta ";
            }
        }

        #riboswitches
        for (my $j = 0; $j < @ribosw; $j += 2) {
            my ($from, $to) = @ribosw[$j, $j + 1];

            ($from, $to) = ($to, $from) if $from > $to;

            for (my $i = $from; $i <= $to; $i++) {
                $color_text .= "$i:royalblue ";
            }
        }

        #trans-splciing
        for (my $j = 0; $j < @transsplicing; $j += 2) {
            my ($from, $to) = @transsplicing[$j, $j + 1];

            ($from, $to) = ($to, $from) if $from > $to;

            for (my $i = $from; $i <= $to; $i++) {
                $color_text .= "$i:cyan ";
            }
        }

        #protein binding motifs
        for (my $j = 0; $j < @rbp_locs; $j += 2) {
            my ($from, $to) = @rbp_locs[$j, $j + 1];

            ($from, $to) = ($to, $from) if $from > $to;

            for (my $i = $from; $i <= $to; $i++) {
                $color_text .= "$i:pink ";
            }
        }


        $color_text =~ s/\s+$//;  # Trim trailing space

        # print "DEBUG: $color_text";
        # Inject JavaScript block to visualize RNA
        print "<script>\n";
        print "  var container;\n";

        print "  window.addEventListener('load', function () {\n";
        print "    var scrollY = window.scrollY;\n";  # prevent scrolling lock to forna
        print "    \n";
        print "    container = new fornac.FornaContainer(\"#rna_ss\", {\n";
        print "      animation: false,\n";
        print "      labelInterval: 50,\n";
        print "      allowPanningAndZooming: true,\n";
        print "      drawBackground: false,\n";
        print "      layout: 'naview'\n";
        print "    });\n";

        print "    var options = {\n";
        print "      structure: '$structure',\n";
        print "      sequence: '$sequence'\n";
        print "    };\n";
        print "    container.addRNA(options.structure, options);\n";
        print "    var colorText = \"$color_text\";\n";
        print "    container.addCustomColorsText(colorText);\n";
        print "    window.scrollTo(0,scrollY);\n";   # prevent scrolling lock to forna
        print "  });\n";
        print "</script>\n";

        # creating legend for visual interpretation 
        print "<div class='legend'>";
        print "<div class='legend-title'>Legend:</div>";
        print "<div class='legend-items'>";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: red; margin-left: 0px;'></span> Motifs\t";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: lightblue; margin-left: 10px;'></span> UTRs\t";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: green; margin-left: 10px;'></span> Exons\t";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: orange; margin-left: 10px;'></span> TRNA\t";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: yellow; margin-left: 10px;'></span> MiRNA\t";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: purple; margin-left: 10px;'></span> SM-site/snRNP-motif\t";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: royalblue; margin-left: 10px;'></span> Riboswitches(s)";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: cyan; margin-left: 10px;'></span> TRANS-splicing\t";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: lime; margin-left: 10px;'></span> PolyA Signal/PolyA Tail\t";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: pink; margin-left: 10px;'></span> Protein Binding Site(s)";

        print "</div>";
        print "</div>";

    } else {
        print "<div class='info-error'>";
        print "<br><b>Maximum folding limit reached</b><br>";
        print "</div>";
    }

    print "</div>";
}

sub location_table {

    sub format_flat_ranges {
        my @data = @_;
        my @formatted;
        for (my $i = 0; $i < @data; $i += 2) {
            last if $i + 1 > $#data;
            push @formatted, "$data[$i] to $data[$i+1]";
        }
        return join(", ", @formatted);  # Comma-separated
    }


    # Begin styling and table

    print "<h2>Feature table</h2>\n";
    print "<table class='table-loc'>";

    print "<tr><th>Structure</th><td>Location(s)</td></tr>\n";

    if (@rna_motif) {
        my $motif_str = format_flat_ranges(@rna_motif);
        print "<tr><th>RNA Motifs</th><td>$motif_str</td></tr>\n";
    }

    if (@exons) {
        my $exons_str = format_flat_ranges(@exons);
        print "<tr><th>Coding Sequence</th><td>$exons_str</td></tr>\n";
    }

    if (@utr) {
        my $utr_str = format_flat_ranges(@utr);
        print "<tr><th>UTRs</th><td>$utr_str</td></tr>\n";
    }

    if (@trna_loc) {
        my $trna_str = format_flat_ranges(@trna_loc);
        print "<tr><th>tRNA(s)</th><td>$trna_str</td></tr>\n";
    }

    if (@mirna_loc) {
        my $mirna_str = format_flat_ranges(@mirna_loc);
        print "<tr><th>miRNA(s)</th><td>$mirna_str</td></tr>\n";
    }

    if (@ribosw) {
        my $ribosw_str = format_flat_ranges(@ribosw);
        print "<tr><th>Riboswitch(s)</th><td>$ribosw_str</td></tr>\n";
    }


    if (@transsplicing) {
        my $trans_str = format_flat_ranges(@transsplicing);
        print "<tr><th>tRNA(s)</th><td>$trans_str</td></tr>\n";
    }

    if (@polyasignal) {
        my $polyasignal_str = format_flat_ranges(@polyasignal);
        print "<tr><th>PolyA motif(s)</th><td>$polyasignal_str</td></tr>\n";
    }

    if (@polyatail) {
        my $polyatail_str = format_flat_ranges(@polyatail);
        print "<tr><th>PolyA motif(s)</th><td>$polyatail_str</td></tr>\n";
    }

    if (@rbp_locs) {
        my $rbp_str = format_flat_ranges(@rbp_locs);
        print "<tr><th>Protein Binding Site(s)</th><td>$rbp_str</td></tr>\n";
    }


    print "</table>\n";
}

print "<script>";
print "  function toggleBoxContent(header) {";
print "    const content = header.nextElementSibling;";
print "    if (content.style.display === 'none') {";
print "      content.style.display = 'block';";
print "    } else {";
print "      content.style.display = 'none';";
print "    }";
print "  }";
print "</script>";




write_file("$TEMPDIR/result.txt", "done\n");

print "</div></main></body></html>";
close $out;
=======
    print "<b>StemGGpair:</b>    none<br>" if ($stemggpairfound==0); 
    @sequ=();
    $str='';
}

sub drawcoloredsequence {
	
	#Here we will create a colored output of the sequence
	#this subroutine requires arrays containing the required information!
	#e.g. @exons for the exons in the format (exonstart,exonend,exonstart,exonend.........)
	# further @transsplicing, @ire, @smsite, @aurichregion, @stemggpairs and further more to be added, see below

	#What to do to invent a new color or feature!
	#1. create an array outside this sub with start and endpoint
	#2. choose in which category it is put. eg FONTCOLOR, UNDERLINED, BOLD/ITALIC
	#3. take care that the start and endpoint are marked in the hash!
	#4. Create an entry below when the nts are counted up, that the feature is writte to the seq-array
	#5. Mark it in the section for the sequence formatting!


	
        my @seq=split ('',$SEQUENCECHECKED);

	#my @exons=(70,80,90,100,160,170);
	#my @transsplicing=(35,58,90,120,125,139, 150,210);
	#	my @ire=(20,25);
	#my @smsite=(34,55,66,77);
	#my @aurichregion=(50,60,88,99);
	#my @stemggpairs=(140,160);

	#print ("<br>ire contains: @ire <br>");
    my %change_underline=();      #the changes will be written into these hashes first!
    my %change_fontcolor=();
    my %change_bold=();
    my %change_italic=();
    my @nooutputpossible=();

    my $c=0;
    #Creating arrays to store, where the layout can't be changed further
    my @possible_underline=@seq; #initialize them to be as big as the seq
    my @possible_fontcolor=@seq;
    #my @possible_bold_italic=@seq;
    #my @possible_italic_only=@seq;
    for ($c=0;$c<=@seq-1;$c++) {
        $possible_underline[$c]=0; #and set them to 0 indicating that this field can be altered
        $possible_fontcolor[$c]=0;     #a number above 0 codes for a color or bold or underlined
	#$possible_bold_italic[$c]=0;
    }


    print '<font face="monospace">';

    ########### ATTENTION !!!!! ####################
    ### Old problem the index of array starts at 0 !!! #####
    ### So for the next lines to be correct we will add #####
    ### an x to the beginning of the array !!!!!!!!!! #########
    ### this won't be printed out, but then the ##########
    ### nomenclature is okay !!!!#####################
    @seq=('x',@seq);
    ####

    ###########################################
    ### Codes for the features in the @possible_ arrays ###
    ### transsplicing = 2   #########################
    ### ire = 3                   #########################
    ### smsite=4               #########################
    ### aurichregion=5     #########################
    ### exonboundary=6   #########################
    ### polyasignal=7

    ###################################################################
    ############ Underlined##############################################
    ###################################################################
    my $mark_in_seq_possible=0;
    my $d=0;
    #print "DEBUG: @transsplicing DEGUEND";
    #Processing the trans-splicing hits
    for ($c=0;$c<=@transsplicing-1;$c=$c+2){
        $mark_in_seq_possible=1;
        for ($d=$transsplicing[$c];$d<=$transsplicing[$c+1];$d++) {     #Okay, this checks whether a part of the seq is marked yet
            $mark_in_seq_possible=0 if ($possible_underline[$d]!=0);
        }
        if ($mark_in_seq_possible==1) {
            #1. mark all fields in the array as nomorepossible
            for ($d=$transsplicing[$c];$d<=$transsplicing[$c+1];$d++) {
                $possible_underline[$d] = 1;
            }
            #Add the lines to the hash
            $change_underline{$transsplicing[$c]}='transstart';
            $change_underline{$transsplicing[$c+1]}='transend';
        }
        else { #If it is already marked
            push (@nooutputpossible,"Trans-splicing hit from $transsplicing[$c] to $transsplicing[$c+1]");
        }
    }

    #Processing the IRE hits
    #print "<p>DEBURG IRE: @ire</p>";
    #print "PETERDEBUGGING";
    for ($c=0;$c<=@ire-1;$c=$c+2){
        $mark_in_seq_possible=1;
        for ($d=$ire[$c];$d<=$ire[$c+1];$d++) {     #Okay, this checks whether a part of the seq is marked yet
            $mark_in_seq_possible=0 if ($possible_underline[$d]!=0);
        }
        if ($mark_in_seq_possible==1) {
            #1. mark all fields in the array as nomorepossible
            for ($d=$ire[$c];$d<=$ire[$c+1];$d++) {
                $possible_underline[$d] = 1;
            }
            #Add the lines to the hash
            $change_underline{$ire[$c]}='irestart';
            $change_underline{$ire[$c+1]}='ireend';
        }
        else { #If it is already marked
            push (@nooutputpossible,"IRE hit from $ire[$c] to $ire[$c+1]");
        }
    }

    ####################################################################################
    ########################## Font Color #################################################
    ####################################################################################
   
    #Processing Promotors
    for ($c=0;$c<=@promotor-1;$c=$c+2){
            $mark_in_seq_possible=1;
            for ($d=$promotor[$c];$d<=$promotor[$c+1];$d++) {     #Okay, this checks whether a part of the seq is marked yet
	             $mark_in_seq_possible=0 if ($possible_fontcolor[$d]!=0);
 	    }
	    if ($mark_in_seq_possible==1) {
        	#1. mark all fields in the array as nomorepossible
                for ($d=$promotor[$c];$d<=$promotor[$c+1];$d++) {
     	           $possible_fontcolor[$d] = 8; #8 means purple
                }
                $change_fontcolor{$promotor[$c]}='promotorstart';
                $change_fontcolor{$promotor[$c+1]}='promotorend';
            }
            else { #If it is already marked
                push (@nooutputpossible,"Promotor from $promotor[$c] to $promotor[$c+1]");
            }
    }

    
    #Processing the PolyASignal
    for ($c=0;$c<=@polyasignal-1;$c=$c+2){
        $mark_in_seq_possible=1;
        for ($d=$polyasignal[$c];$d<=$polyasignal[$c+1];$d++) {     #Okay, this checks whether a part of the seq is marked yet
            $mark_in_seq_possible=0 if ($possible_fontcolor[$d]!=0);
        }
        if ($mark_in_seq_possible==1) {
            #1. mark all fields in the array as nomorepossible
            for ($d=$polyasignal[$c];$d<=$polyasignal[$c+1];$d++) {
                $possible_fontcolor[$d] = 7;
            }
            $change_fontcolor{$polyasignal[$c]}='polyasignalstart';
            $change_fontcolor{$polyasignal[$c+1]}='polyasignalend';
        }
        else { #If it is already marked
            push (@nooutputpossible,"PolyA-Signal from $polyasignal[$c] to $polyasignal[$c+1]");
        }
    }
    #Processing the smsites
    for ($c=0;$c<=@smsite-1;$c=$c+2){
        $mark_in_seq_possible=1;
        for ($d=$smsite[$c];$d<=$smsite[$c+1];$d++) {     #Okay, this checks whether a part of the seq is marked yet
            $mark_in_seq_possible=0 if ($possible_fontcolor[$d]!=0);
        }
        if ($mark_in_seq_possible==1) {
            #1. mark all fields in the array as nomorepossible
            for ($d=$smsite[$c];$d<=$smsite[$c+1];$d++) {
                $possible_fontcolor[$d] = 5;
            }
            $change_fontcolor{$smsite[$c]}='smsitestart';
            $change_fontcolor{$smsite[$c+1]}='smsiteend';
        }
        else { #If it is already marked
            push (@nooutputpossible,"Sm-Site from $smsite[$c] to $smsite[$c+1]");
        }
    }

    #Processing the au-rich-regions
    for ($c=0;$c<=@aurichregion-1;$c=$c+2){
        $mark_in_seq_possible=1;
        for ($d=$aurichregion[$c];$d<=$aurichregion[$c+1];$d++) {     #Okay, this checks whether a part of the seq is marked yet
            $mark_in_seq_possible=0 if ($possible_fontcolor[$d]!=0);
        }
        if ($mark_in_seq_possible==1) {
            #1. mark all fields in the array as nomorepossible
            for ($d=$aurichregion[$c];$d<=$aurichregion[$c+1];$d++) {
                $possible_fontcolor[$d] = 5;
            }
            #Add the lines to the hash
            $change_fontcolor{$aurichregion[$c]}='aurichregionstart';
            $change_fontcolor{$aurichregion[$c+1]}='aurichregionend';
        }
        else { #If it is already marked
            push (@nooutputpossible,"AU-rich region from $aurichregion[$c] to $aurichregion[$c+1]");
        }
    }
    #Processing stemggpairs
    for ($c=0;$c<=@stemggpairs-1;$c=$c+2){
        $mark_in_seq_possible=1;
        for ($d=$stemggpairs[$c];$d<=$stemggpairs[$c+1];$d++) {     #Okay, this checks whether a part of the seq is marked yet
            $mark_in_seq_possible=0 if ($possible_fontcolor[$d]!=0);
        }
        if ($mark_in_seq_possible==1) {
            #1. mark all fields in the array as nomorepossible
            for ($d=$stemggpairs[$c];$d<=$stemggpairs[$c+1];$d++) {
                $possible_fontcolor[$d] = 4;
            }
            #Add the lines to the hash
            $change_fontcolor{$stemggpairs[$c]}='stemggpairsstemstart';
            $change_fontcolor{$stemggpairs[$c]+1}='stemggpairsstemend';
            $change_fontcolor{$stemggpairs[$c]+2}='stemggpairsbulge';
            $change_fontcolor{$stemggpairs[$c]+3}='stemggpairsstemstart';
            
            $change_fontcolor{$stemggpairs[$c+1]-3}='stemggpairsstemend';
            $change_fontcolor{$stemggpairs[$c+1]-2}='stemggpairsbulge';
            $change_fontcolor{$stemggpairs[$c+1]-1}='stemggpairsstemstart';
            $change_fontcolor{$stemggpairs[$c+1]}='stemggpairsstemend';
        }
        else { #If it is already marked
            push (@nooutputpossible,"Stem-GG-Pair from $stemggpairs[$c] to $stemggpairs[$c+1]");
        }
    }
     

    ############################################################################
    ########################### Bold ###########################################
    ############################################################################
    #Processing the exons
    for ($c=0;$c<=@exons-1;$c=$c+2){
        #Add the lines to the hash
        $change_bold{$exons[$c]}='exonstart';
        $change_bold{$exons[$c+1]}='exonend';
    }

    ############################################################################
    ########################### Italic #########################################
    ############################################################################
    #Processing the untranslated regions UTR
    for ($c=0;$c<=@utr-1;$c=$c+2){
        #Add the lines to the hash
        $change_italic{$utr[$c]}='utrstart';
        $change_italic{$utr[$c+1]}='utrend';
    }

    #print $rv->h3("");
    print $rv->h1("Summary");
    
    #Print out the predicted protein! reformat for a nicer looking on the page	
    &predprotein;


    #Now create the html output!!!!!!!!!!!!

    print "<br><b>Colored Sequence:</b><br>";

    my $ntausgabe=1;

    my $boldonout=0;
    my $fontcolor=0; #0=black 1=red 2=blue 3=navy 4=fuchsia
    my $underlined=0; #0=no 1=yes

    for ($c=1;$c<=@seq-1;$c++){
        ###########################    
        ### OUTPUT UNDERLINED ######
        ###########################
        if (defined $change_underline{$c}) {    
            if ($change_underline{$c} eq 'transstart' || $change_underline{$c} eq 'irestart')  {
                $seq[$c]='<u>'.$seq[$c];
                $underlined=1;
            }
            if ($change_underline{$c} eq 'transend' || $change_underline{$c} eq 'ireend') {
                $seq[$c]=$seq[$c].'</u>';
                $underlined=0;
            }
        }       
        
        ###########################
        ###### OUTPUT BOLD #########
        ###########################
        if (defined $change_bold{$c}) {
            if ($change_bold{$c} eq 'exonstart') {
                $seq[$c]='<b>'.$seq[$c];
                $boldonout=1;
            }
            if ($change_bold{$c} eq 'exonend') {
                $seq[$c]=$seq[$c].'</b>';
                $boldonout=0;
            }
        }
        
       
	###########################
	###### OUTPUT ITALIC ######
	###########################
	if (defined $change_italic{$c}) {
             if ($change_italic{$c} eq 'utrstart') {
                     $seq[$c]='<i>'.$seq[$c];
                     $italiconout=1;
             }
             if ($change_italic{$c} eq 'utrend') {
                     $seq[$c]=$seq[$c].'</i>';
                     $italiconout=0;
             }
        }
	
        ############################
        ####### OUTPUT COLOR ########
        ############################
        if (defined $change_fontcolor{$c}) {
		

       	   ######## Promotor ##########
           if ($change_fontcolor{$c} eq 'promotorstart') {
                 $seq[$c]='<font color=purple>'.$seq[$c];
                 $fontcolor=8;
           }
           if ($change_fontcolor{$c} eq 'promotorend') {
                $seq[$c]=$seq[$c].'</font>';
                $fontcolor=0;
           }
           ######## PolyA-Signal ##########
           if ($change_fontcolor{$c} eq 'polyasignalstart') {
               $seq[$c]='<font color=lime>'.$seq[$c];
               $fontcolor=5;
           }
           if ($change_fontcolor{$c} eq 'polyasignalend') {
                 $seq[$c]=$seq[$c].'</font>';
                 $fontcolor=0;
           }
	   ######## Sm-Site ##########
           if ($change_fontcolor{$c} eq 'smsitestart') {
                $seq[$c]='<font color=red>'.$seq[$c];
                $fontcolor=1;
            }
            if ($change_fontcolor{$c} eq 'smsiteend') {
                $seq[$c]=$seq[$c].'</font>';
                $fontcolor=0;
            }
            ####### AU-rich region #######
            if ($change_fontcolor{$c} eq 'aurichregionstart') {
                $seq[$c]='<font color=blue>'.$seq[$c];
                $fontcolor=2;
            }
            if ($change_fontcolor{$c} eq 'aurichregionend') {
                $seq[$c]=$seq[$c].'</font>';
                $fontcolor=0;
            }
            ####### Stem-GG-Pairs #######
            if ($change_fontcolor{$c} eq 'stemggpairsstemstart') {
                $seq[$c]='<font color=fuchsia>'.$seq[$c];
                $fontcolor=4;  
            }
            if ($change_fontcolor{$c} eq 'stemggpairsstemend') {
                $seq[$c]=$seq[$c].'</font>';
                $fontcolor=0;
            }
            if ($change_fontcolor{$c} eq 'stemggpairsbulge') {
                $seq[$c]='<font color=green>'.$seq[$c].'</font>';
            }
        }
        
     
        #############################
        ###Formatting the Output ########
        ### into suitable lines ###########
        #############################
        
        
        if ($c % 50 == 0) {   
            if ($boldonout==1) {
                $seq[$c]=$seq[$c].'</b>';
                $seq[$c+1]='<b>'.$seq[$c+1];
            }
	    if ($italiconout==1) {
		$seq[$c]=$seq[$c].'</i>';
		$seq[$c+1]='<i>'.$seq[$c+1]; 
	    }
            if ($underlined==1) {
                $seq[$c]=$seq[$c].'</u>';
                $seq[$c+1]='<u>'.$seq[$c+1];
            }
            if ($fontcolor != 0) {
                $seq[$c]=$seq[$c].'</font>';
                $seq[$c+1]='<font color = red>'.$seq[$c+1] if ($fontcolor==1);
                $seq[$c+1]='<font color = blue>'.$seq[$c+1] if ($fontcolor==2);
                $seq[$c+1]='<font color = fuchsia>'.$seq[$c+1] if ($fontcolor==4);
		$seq[$c+1]='<font color = lime>'.$seq[$c+1] if ($fontcolor==5);
		$seq[$c+1]='<font color = purple>'.$seq[$c+1] if ($fontcolor==8);
            }
        }
    }

    #for debugging only
    #$seq[30]='<span style="underline-color: yellow">'.$seq[30];
    #$seq[45]=$seq[45].'</span>';
    ################


    for ($c=1;$c<=@seq-1;$c++){
        print "$seq[$c]";
        if ($c % 50 == 0) {
            print "     ".$c."<br>";
        }
        if ($c %10 == 0) {
            print "";
        }
    }
    #print;
    @seq=();
    print "<br>Output not possible:<br>@nooutputpossible<br>" if (@nooutputpossible>0);
    
    if (length $SEQUENCECHECKED<=$MAXFOLDINGLEN) {
	&drawcoloredstructure;
    }

    
    print "<br><br>Legend:<br><b>BOLD</b> marks EXONS<br><u>UNDERLINED</u> marks IRE or TRANSSPLICNG hits<br>";
    print "<i>ITALIC</i> marks putative UTRs<br>";
    print "<font color=red>RED </font> marks SMSITES or snRNP-binding motifs<br><font color=blue>BLUE</font> marks AU-rich regions<br>";
    print "<font color=fuchsia>FUCHSIA</font> with <font color=green>GREEN</font> marks Stemm-GG-Pairs<br>";
    print "<font color=lime>LIME</font> marks PolyA-signal<br><font color=purple>PURPLE</font> marks a Promotor<br>";
    print "<br>Pr.A1.bin.site = Protein A1 binding site<br>";

	
    
    #print "<br><br>%change_fontcolor<br><br>";
}


sub calcUTR {
	#Trying to calculate 3' and 5' regions! But this is really very basic done
	
        #	$dnara='DNA|RNA';
	#print "<br>Exons: @exons<br>";
	@utrprintout=();
	$numberofexons=@exons/2;
	#if ($dnarna eq 'RNA' && $numberofexons==1) { #wenn es mehr exons sind, dann stimmt etwas nicht! und es erfolgt keine Ausgabe!
        
		#then write the 3' and 5' UTR
	    #@utr=(@utr,1,$exons[0]-1,$exons[1]+1,(length $SEQUENCECHECKED)); #stimmt das auch ??
	    #print "DEBUGGING2";

    if ($numberofexons==1){
	###### RNA 5' and 3' Detection ###########
	############### S T A R T ################



	###New: checking for coding sequence!!! Perhaps this should better be done in UTR, but nevermind
	@new3primeutr=();
	@new5primeutr=();
	@singleexongenscan=();
	@singleexonboundaries=();
	@singleexongenscan=`grep "Sngl +" $TEMPDIR/$job.genscanout`;
	#okay we have obtained line if a single exon is present!!!!
	if (@singleexongenscan>0) {
		foreach $testline (@singleexongenscan) {
			$testline=~m/Sngl [+][ ]+([0-9]+)[ ]+([0-9]+)/;
			push @singleexonboundaries,$1,$2;

		}
	}
	@polyagenscan=();
	@polyagenscanstart=();
	@polyagenscan=`grep "PlyA +" $TEMPDIR/$job.genscanout`;
	#okay we have obtained line if a single exon is present!!!!
	if (@polyagenscan>0) {
		foreach $testline (@polyagenscan) {
			$testline=~m/PlyA [+][ ]+([0-9]+)[ ]+([0-9]+)/;
			#@polyagenscanstart=@polyagenscanstart,$1;
			push @polyagenscanstart,$1;
			
		}
		
	}	# if it works we have the cds boundaries and the polasignal
	for ($count=0;$count<=@singleexonboundaries-1;$count=$count+2){
		#$singleexonboundaries[$count+1] 3'
		$count2=0;
		#		print "D1";
		$leaveendless=0;
		while ($leaveendless<5000 && defined $polyagenscanstart[$count2] && defined $singleexonboundaries[$count] &&($polyagenscanstart[$count2]<$singleexonboundaries[$count])) {
			$rightpolyaindex=$count2+1;
			$count2++;
			$endless++;
			#	print "D2";
		}
		#jetzt gilt: 3' UTR von $singleexonboundaries[$count] bis $polyagenscanstart[$rightpolyindex]
		push @new3primeutr,$singleexonboundaries[$count+1]+1,$polyagenscanstart[$rightpolyaindex]+5 if ($polyagenscanstart[$rightpolyaindex]>1);	
		#oder man sagt, die 3' region geht dann bis zum Schluss, koennte man auch sagen, das soll aber mal okay sein!
	}
	#now 5' UTR works only with the first sngl exon because here we can located the 5'UTR Start to 1 bzw 0, otherwise we can't, because we
	#do not know where ist starts, is it the PolyASignal, PolyATail or elsewhere ??
	if ($singleexonboundaries[0]==$exons[0]){ #the it might be rna and the 5'UTR is 1 to $exons[0]-1
		push @new5primeutr,1,$exons[0]-1;
	}
	
	#print "<br>New3' @new3primeutr<br>PlyAstasrt @polyagenscanstart Bound: @singleexonboundaries<br>";
	#print "Exons: @exons Exonsend";
	#print "<br>5' UTR: @new5primeutr End 5' UTR! <br>";




	#################################################
	######################## E N D ##################
	#################################################


     }
	#}
	#print "<br>Debugging: $numberofexons DR: $dnarna<br>";
        #######################################
	#if (($dnarna eq 'DNA') || ($dnarna eq 'unknown')) { #is it really DNA ??
        if (($dnarna eq 'DNA') && ($numberofexons>1)){    
	   for (my $c=0;$c<=@polyasignal-1;$c=$c+2){
		#Wenn wir also ein PolyASignal haben
                #Ausprobieren ob das geht
                #pos($SEQUENCECHECKED)=$polyasignalstart;
		#print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
		$newdnarna='';        
                while ($SEQUENCECHECKED=~/aaaaaaaaaaa/g) { #Wir suchen jetzt nach einem moeglichen PolyA-Tail!
                    $polyatailstart=pos($SEQUENCECHECKED); #wenn wir sie jetzt wirklich haben
		    #print "<br> PAT: $polyatailstart";
		    if ($polyatailstart-$polyasignal[$c]<=40 && $polyatailstart-$polyasignal[$c]>5) {#nimmt also nur eine aaaaa Sequ die nahe am polyasignal dran ist.
                        #dann ist es wohl RNA und wir koennen die 3' UTR markieren
                        $newdnarna='RNA';
			#print "<br> DEBUGGING: Polyasignal mit PolyATail!!! <br>";
			#Now we will look for the last exons where the 3' UTR starts!!
			foreach $exonfind (@exons){
				$actual3prime=$exonfind if ($exonfind<$polyasignal[$c]);
			}
			@utr=(@utr,$actual3prime+1,$polyasignal[$c]+5);
			@utrprintout=(@utrprintout,3,$actual3prime+1,$polyasignal[$c]+5);
			@new3primeutr=(@new3primeutr,$actual3prime+1,$polyasignal[$c]+5);
			last;
                    }
                }
                $newdnarna='DNA' if ($newdnarna eq '');
            }
            if (@promotor>0) {
                for ($c=0;$c<=@promotor-1;$c=$c+2){
			#Dann muss es DNA sein und Promotor bis InitialExon = 5'UTR
			$newdnarna='DNA';
			#find next initial exon
			foreach $exonfind (@exons){
	                	$actual5prime=$exonfind if ($exonfind<$polyasignal[$c]);
				if ($exonfind>$promotor[$c+1]) {
					$actual5prime=$exonfind;
					last;
				}
                	}
			@utr=(@utr,$promotor[$c],$actual5prime-1);
			@utrprintout=(@utrprintout,5,$promotor[$c],$actual5prime-1);
			@new5primeutr=(@new5primeutr,$promotor[$c],$actual5prime-1);
		}          
       	   }
       }
	#### OLD
	#if (@utrprintout>0) {
		#	print "<b>UTR:</b>           start  -   end<br>";	
		#		for ($c=0;$c<=@utrprintout-1;$c=$c+3) {
			#	printf (" %1d'            %-6d - %6d<br>",$utrprintout[$c+0],$utrprintout[$c+1],$utrprintout[$c+2]); 
			#		}
			#	}
			#	else {
				#	print "<b>UTR:</b>           none detected<br>";
				#	}

	#### NEW
	
	print "<b>UTR:</b>           start  -   end   -  stems - energy<br>";	
	if (@new5primeutr>0){
		#print "<b>5' UTR:</b>        start  -   end<br>";
		for ($c=0;$c<=@new5primeutr-1;$c=$c+2) {
			printf (" 5'            %-6d - %6d",$new5primeutr[$c+0],$new5primeutr[$c+1]);
			#Testing the UTR Folding!!!
			my $temp=substr($SEQUENCECHECKED,$new5primeutr[$c+0],$new5primeutr[$c+1]-$new5primeutr[$c+0]);
			my @returnout=&checkstemsonly($temp,1);
			#print "<br>Attention: @returnout ENDATTENTION<br>";
			print "       $returnout[0]" if ($returnout[0] == $returnout[1]);
			print "       $returnout[0]-$returnout[1]" if ($returnout[0] != $returnout[1]);
			print "     $returnout[2]<br>" if ($returnout[2] != 1);	
			
		}
	}
	else {
		print "<b>5' UTR:</b>        none detected<br>";
	}
	
	if (@new3primeutr>0){
		#print "<b>3' UTR:</b>        start  -   end<br>";
		for ($c=0;$c<=@new3primeutr-1;$c=$c+2) {
			printf (" 3'            %-6d - %6d",$new3primeutr[$c+0],$new3primeutr[$c+1]);
			#Testing the UTR Folding!!!
			my $temp=substr($SEQUENCECHECKED,$new3primeutr[$c+0],$new3primeutr[$c+1]-$new3primeutr[$c+0]);
			my @returnout=&checkstemsonly($temp,1);
			#print "<br>Attention: @returnout ENDATTENTION<br>";
			#print "Posted: $temp<br>";
			print "       $returnout[0]" if ($returnout[0] == $returnout[1]);
			print "       $returnout[0]-$returnout[1]" if ($returnout[0] != $returnout[1]);
			print "     $returnout[2]<br>" if ($returnout[2] != 1);	
		}
		if ($returnout[0]+$returnout[1]>=7) { #means there are more than 3,5 stem loops!!!!!
			print "<br>Potential stability elements might be located in this 3' UTR !!! <br>";
		}
	}
	else {
		print "<b>3' UTR:</b>        none detected<br>";
	}
	
	#Now we shall create array for the UTR for the colored sequence
	@utr=(@utr,@new5primeutr,@new3primeutr);
	@utr=sort {$a <=> $b} @utr;

		


	
}

sub protA1bisite  {

	#finding the protein A1 binding motif
	#allowing 5 mismatches but no insertions/deltions
	#
	my $line='cuggauuauucaacugaaugccucacucagagaaugaa';
	my @protA1bisi=split('',$line);
	#my $se='nnnnnnnnnnnncuggauuauucaacugaaugccucacucagagaaugaannnnnnnnnnnnnnncuggauuauucaacugaaugccucnacucagagaaugaannnnnnnnnnnnnnnnnnnnnnnnnnn';
	my @sequence=split('',$SEQUENCECHECKED);
	my $prota1startingline=0;
	#allow 5 mismatches out of 38!
	PROA1: for (my $wh=0;$wh<=@sequence-38;$wh++){
	        my $pa1mismatches=0;
	        for (my$c=0;$c<=37;$c++) {
	                $pa1mismatches++ if ($sequence[$wh+$c] ne $protA1bisi[$c]);
	                next PROA1 if ($pa1mismatches>5);
	        }
	        if ($pa1mismatches<=5) {#this is a hit
                	print "<b>Pr.A1 bin.site:</b>start  -  mismatch  - seq<br>" if ($prota1startingline==0);
			$prota1startingline=1;
			my $pa1seq=substr($SEQUENCECHECKED,$wh,38);;			
			printf ("               %-5d        %-2d       %s<br>",$wh,$pa1mismatches,$pa1seq);
			
        	}
	}
	print "<b>Pr.A1 bin.site:</b>none<br>" if ($prota1startingline==0);
}

sub promotor {
	#Have a look at a promotor
	@temparray=`grep "Prom" $TEMPDIR/$job.genscanout`;
	if (@temparray==0) {
		 print "<b>Promotor:</b>      none<br>";
	}
	else {
		print "<b>Promotor:</b>      start  -   end<br>";
		foreach $templine (@temparray) {
			$templine=~/[0-9. ]+(Prom) ([+-])[ ]+([0-9]+)[ ]+([0-9]+)/;
			if ($1 eq 'Prom' && $2 eq '+') {
				#A promotor is present
				printf (" Prom:         %-6d - %6d<br>",$3,$4);#  found on $2 strand from $3 to $4"); #we are now loosing information! $2 contains the strand
				@promotor=(@promotor,$3,$4);
			}
		}
	}


}

sub polyAsignal {
	#have a look if a polyASignal has been found!!!
	@polyAlines=`grep "PlyA" $TEMPDIR/$job.genscanout`;
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

sub ARE {
    $arepresent=0;
    $are_pos=1;
    #Check for so called ARE = Au-rich regions; consensus (AUUUA)n of ~50 bases

    while ($SEQUENCECHECKED=~/([ag]uuu[ag](uuu[ag])+)/g) {
        $are_len=length($1);
        $are_pos=pos($SEQUENCECHECKED)-$are_len;
        if ($are_len >=9){
	    printf ("<b>ARE</b>:          %-6d - %6d   possi. match:   %s<br>",
		    $are_pos,$are_pos+$are_len-1,
		    $1);
	    $arepresent=1;
            $mismatchinare=0;
            @aretemp=split('',$1);
            for ($arecount=0;$arecount<@aretemp;$arecount++){
                $mismatchinare++ if ($aretemp[$arecount] eq 'g');
            }
	    if ($mismatchinare) {
		printf ("<b>ARE</b>:          %-6d - %6d   mismatch:  %2d<br>",
			$are_pos,($are_pos+$are_len-1),$mismatchinare);
	    }
	    @aurichregion=(@aurichregion,$are_pos+1,$are_pos+$are_len);
        }
    }

    if ($arepresent==0) {
        print "<b>ARE:</b>           None       *(AU-rich region of at least 30 nt)<br>";
    }
}
############ARE_old is peters script#############
sub ARE_old {
	$arepresent=0;
	$arepos=1;
	#Check for so called ARE = Au-rich regions; consensus (AUUUA)n of ~50 bases
	$leadlineprinted=0;
	while ($SEQUENCECHECKED=~m/([ag]uuu[ag](uuu[ag])+)/g) {
		$arelen=length $1;	
		$arepos=pos($SEQUENCECHECKED)-$arelen+1; #HIIIIIEEEEELLLLLFFFFFEE ! Funktioniert nicht, keine Ahnung warum!
		if ($arelen >=9){ 
			print "<b>ARE:</b>           none       *(AU-rich region of at least 10 nt)<br>" if ($leadlineprinted==0);
			$leadlineprinted=1;
			$mismatchinare=0;
			@aretemp=split('',$1);
			for ($arecount=0;$arecount<@aretemp;$arecount++){
				$mismatchinare++ if ($aretemp[$arecount] eq 'g');
			}
			printf (" ARE:          %-6d - %6d   mismatch:  %2d<br>",$arepos,($arepos+$arelen-1),$mismatchinare);
			$arepresent=1;
			@aurichregion=(@aurichregion,$arepos,$arepos+$arelen-1);
		}
	}
	if ($arepresent==0) {
		print "<b>ARE:</b>           none       *(AU-rich region of at least 30 nt)<br>";
	}
}


sub tRNA {
	######################################################################
	#Looking for tRNAs using tRNAscan-SE

	$answertrnascan=`$TRNASCANFOLDER/tRNAscan-SE -Q -y -f $TEMPDIR/$job.trnascanout $TEMPDIR/$job.genscan`;
	open (TRNA,"$TEMPDIR/$job.trnascanout");
	$line=<TRNA>;
	if ($line=~/Length/){
		print "<b>tRNA<sup>2</sup>:</b><br>";
		print "$line<br>";
		while ($line=<TRNA>){
			print "$line<br>";
		}
	}
	else {
		print "<b>tRNA<sup>2</sup>:</b>         none<br>";
	}
}

sub smsite {
	$smlength=-1;
	$smpos=-1;
        $leadlineprinted=0; #the first line not yet printed
	while ($SEQUENCECHECKED=~/([ag][ag](u+([agc]?)u+)[ag][ag])/g){
		$smlength=length $1;
		$smpos=pos($SEQUENCECHECKED)-$smlength+1;
		
		if (length $3 == 1 && length $2 >=4) {
			#print "Potential snRNP binding motif, similar to a sm-site, at position $smpos with the sequence $1 <br>";
			print "<b>snRNP-motifs:</b>  start      sequence                quality<br>" if ($leadlineprinted==0);
			$leadlineprinted=1;
			printf (" snRNP-motif:  %-6d     %-20s          +<br>",$smpos,$1);
			@smsite=(@smsite,$smpos,$smpos+((length $1)-1));
		}
		if (length $3 == 0 && length $2>=4) {
			print "<b>snRNP-motifs:</b>  start      sequence                quality<br>" if ($leadlineprinted==0);
			$leadlineprinted=1;
			printf (" Put. sm-site: %-6d     %-20s          ++<br>",$smpos,$1);
			@smsite=(@smsite,$smpos,$smpos+((length $1)-1));
		}
	}
	print "<b>snRNP-motifs:</b>  none<br>" if ($leadlineprinted==0);
	
	##### OUTPUT if Seq is RNA has no cds but smsite --> structured perhaps catalytic RNA !

	if ($grepanswer=~/NO EXONS/ && $ORIGINchecked==1 && @smsite>0){
		print "<br>As I could not detect a coding sequence on this RNA, but there are 1 or more sn-RNP motifs (sm-sites),<br>it might be possible that this is a catalytic RNA!!<br>";
	}
}

sub exons {
    
	#print $rv->h3("Results of the exon/gene search");	
	#print "Exons:&nbsp;&nbsp;";
	$grepanswer=`grep "NO EXONS" $TEMPDIR/$job.genscanout`;

	#print $rv->p("No gene/exon detected") if ($grepanswer=~/NO EXONS/);
	#print "none<br>" if ($grepanswer=~/NO EXONS/);
	if ($grepanswer=~/NO EXONS/) {
		print "<b>Exons<sup>1</sup>:</b>        none<br>";
	}
	else {
		print "<b>Exons<sup>1</sup>:</b>        start  -   end       type\n";
	}
	print $grepanswer;
	$genscan = Bio::Tools::Genscan->new(-file => "$TEMPDIR/$job.genscanout");
	while(my $gene = $genscan->next_prediction()) {
	    my @exonspresent=$gene->exons();	    
	    foreach my $exon ($gene->exons()) {
		my $loc = $exon->location;
		#print "Coding exon: " if ($exon->is_coding()==1);
		#print "Noncoding exon: " if ($exon->is_coding()==0);
		printf (" Exon:         %-6d - %6d",$loc->start, $loc->end);
		print "     coding<br>" if ($exon->is_coding()==1);
		print "     noncoding<br>" if ($exon->is_coding()==0);
		@exons=(@exons,$loc->start,$loc->end); #adding to the drawseq array!		
	    }

	    foreach my $utr ($gene->utrs()) {#What the hell does this print out ????
		printf ("UTR:           %-6d - %6d     ",$loc->start,$loc->end);
		print $utr->primary_tag, " <br>";
	    }
	    #print "---------------------------------\n";
		my $prot = $gene->predicted_protein;
	   	@predprot=split ('',$prot->seq);		
	}

# comment out by liang 

if ($grepanswer=~/NO EXONS/ && $ORIGINchecked==1) {

		#then we have an RNA without a coding region
		#were we will scan for highly structured regions doing it that way:
		#we start at the beginning an take 150nt, fold them, calc the stems, if above a certain level, we will print out a message!
		my $count=0;
		my $len=length $SEQUENCECHECKED;
		my $query='';
		my @answer=();
		my $printout=0;
		my $last=0;
		print "<b>Region:        From     To    Stems   Energy Remark</b><br>";
		for ($count=0;$count<=$len;$count=$count+100){
			if ($count+150<$len){
				$query=substr($SEQUENCECHECKED,$count,150);
				printf ("<br> Pos:       %6d - %6d",$count,$count+150);
			}
			else {
				$query=substr($SEQUENCECHECKED,$count,$len-$count);
	           		printf ("<br> Pos:       %6d - %6d",$count,$len);	
				$last=1;
			}
			@answer=&checkstemsonly($query,1);
			#print "<br>ANswer: @answer";
			print "      $answer[0]" if ($answer[0] == $answer[1]);
			print "    $answer[0]-$answer[1]" if ($answer[0] != $answer[1]);
			print "    $answer[2]";
			if ($answer[0]>=3 && $answer[1]>=3){ #wenn mehr oder genau 3 stems / 150 nt
				print "   **<br>";
				$printout=1;
			}
			else {
				print "<br>";
			}
			last if ($last==1);
		}
		if ($printout==1) {
			print "** Three or more stem loops in this region! This is a highly structured region. <br>Please check whether tRNA, rRNA or another highly structured RNA is encoded here!<br>";
		}
		else {
			print "<br>";
		}
	}
}

sub createpicture {
# comment out by liang
		@answerout=`$VIENNAFOLDDIR/Fold $job.seq`;
		#print $rv->p("@answerout");
		#print $rv->p("Answerd:<br>@answerout");  #3_ss.ps		
		
		#Evaluate the number of stems in the structure
		
		$joinedanswer=join('',@answerout);
		$joinedanswer=~s/\n//g;
		$joinedanswer=~/([.(][.()]+) (\([ -]*[0-9]+\.[0-9]+\))/;
		@structure=split('',$1);
		$energy=$2;
		$energy=~s/[()]//g;
		
		##### Creating a unique id using md5sum
# comment out by liang
		$md5in=`md5sum $TEMPDIR/$job.seq`;
		$md5in=~/([0-9abcdef]{32})/;
		$md5sum=$1;
		
		$actualdir=cwd();
		chdir "$TEMPDIR";
		
		#Create the picture !!

# comment out by liang
# because I will use a seperated cgi to complete it
# for changed vienna package (our version)
#		$answerconvert=system("convert "."$job".'_ss.ps'.' -crop 0x0'." $job"."_"."$md5sum".".jpg");

#		$answerconvert=system("convert rna.ps -crop 0x0 $job"."_$md5sum.jpg");  # for original vienna package
		#print $rv->p("Answerconvert: $answerconvert");  #3_ss.ps
#                $answermove=system("mv $job.jpg $TEMPPICSDIR/$job.jpg");
		#print $rv->p("Answermove: $answermove"); 
		
		#next line will be commented out, so that it can be used at an other postiion in the proggi
		#print '<img src="'."../results/$job.jpg".'">';
		
		#Print this information right after the picture !!!
		print "<br><b>Length:</b>        $SEQUENCELENGTH";
	        print "     *some information is only available up to $MAXFOLDINGLEN nt</font>\n" if ($SEQUENCELENGTH >$MAXFOLDINGLEN); 
	        print "<br><b>Origin:</b>        $dnarna<br>";
}

sub checkstems {
		#evaluate the structure
		$fldklauf=0; #Klammern auf pro stem
		$fldstemsauf=0;
		$fldklzu=0;
		$fldstemszu=0;

		for ($i=0;$i<=@structure-1;$i++) {
		    $fldklauf++  if ($structure[$i] eq '(');
		    $fldstemsauf++ if ($structure[$i] eq ')' && $fldklauf>=5);
		    $fldklauf=0 if ($structure[$i] eq ')');
		    
		    $fldklzu++  if ($structure[$i] eq ')');
		    $fldstemszu++ if ($structure[$i] eq '(' && $fldklzu>=5);
		    $fldklzu=0 if ($structure[$i] eq '(') 
		} 
		#Now a correction for Stemzu!
		$fldstemszu++ if ($fldklzu>=5);
		print "<b>Energy:</b>        $energy kcal/mol<br>";
		print "<b>Stems:</b>         $fldstemsauf Stem-Structure/s<br>" if ($fldstemsauf==$fldstemszu);
		print "<b>Stems:</b>         $fldstemsauf-$fldstemszu Stem-Structure/s<br>" if ($fldstemsauf<$fldstemszu);
		print "<b>Stems:</b>         $fldstemszu-$fldstemsauf Stem-Structure/s<br>" if ($fldstemsauf>$fldstemszu);
		if ($fldstemsauf != 0 && $fldstemszu != 0) {
			if ((2*@structure)/($fldstemsauf+$fldstemszu)<60) { #The last number determines when the structure is interesting! 60 means at least on Stem in each 60 nt!
				print "<br>The sequence seems to contain a lot of secondary structure. If the RNA structure search below<br>does not find a result, it might be interesting to have a closer look at the structures.<br>";
				print 'You might find it useful to look in the book <br><br>"RNA Motifs and Regulatory Elements"<br>Thomas Dandekar (Ed.)<br>Published by Springer<br>ISBN 3-540-41701<br>';
			}
		}
		#An message if it is highly structured that it could be an rRNA
		if (($fldstemauf+$fldstemzu)/2>=10) {
			print "<br>          Highly structured RNA, could this be a ribosomal RNA ?<br>";
		}
		#Okay, that's it for the stem detection!
}

sub checkstemsonly {
		#evaluate the structure
		my $inseq=$_[0];
		my $inwhattodo=$_[1];   # --> 0 means is structure, evaluate, 1--> sequence, please fold first!
		my $struct;
		my @structure;
		my $energy;
		if ($inwhattodo==1 && length $inseq<=$MAXFOLDINGLENUTR){
			@struct=`echo $inseq | $VIENNARNAFOLDDIR/RNAfold`;
			$struct[1]=~/([().]+) \(([+-. 0-9]+)\)/;
			@structure=split('',$1);
			$energy=$2;
			#print "DEBUGMARK1: $struct[1]";
		}
		if ($inwhattodo==1 && length $inseq>$MAXFOLDINGLENUTR){
			return ("Too long for detection",1); #the 1 shows that we did not produce a result (not negative due to negative energy)
		}
		if ($inwhattodo==0){
			@structure=split ('',$inseq);
			$energy=0;
		}
		

		#print "Debug: Here is thew struct: @structure END";	
		my $fldklauf=0; #Klammern auf pro stem
		my $fldstemsauf=0;
		my $fldklzu=0;
		my $fldstemszu=0;

		for (my $i=0;$i<=@structure-1;$i++) {
		    $fldklauf++  if ($structure[$i] eq '(');
		    $fldstemsauf++ if ($structure[$i] eq ')' && $fldklauf>=5);
		    $fldklauf=0 if ($structure[$i] eq ')');
		    
		    $fldklzu++  if ($structure[$i] eq ')');
		    $fldstemszu++ if ($structure[$i] eq '(' && $fldklzu>=5);
		    $fldklzu=0 if ($structure[$i] eq '(') 
		} 
		#Now a correction for Stemzu!
		$fldstemszu++ if ($fldklzu>=5);

		return ($fldstemsauf,$fldstemszu,$energy) if ($fldstemsauf==$fldstemszu);
		return ($fldstemsauf,$fldstemszu,$energy) if ($fldstemsauf<$fldstemszu);
		return ($fldstemszu,$fldstemsauf,$energy) if ($fldstemsauf>$fldstemszu);

		print "ERROR ERROR!!!!  This line should N O T be reached!!!";
	
		#######ATTENTION THESES LINES BELOW ARE NOT PRINTED!!!! BECAUSE WE RETURN ABOVE !!!!!!#####
		###########################################################################################
		###########################################################################################


		print "<b>Energy:</b>        $energy kcal/mol<br>";
		print "<b>Stems:</b>         $fldstemsauf Stem-Structure/s<br>" if ($fldstemsauf==$fldstemszu);
		print "<b>Stems:</b>         $fldstemsauf-$fldstemszu Stem-Structure/s<br>" if ($fldstemsauf<$fldstemszu);
		print "<b>Stems:</b>         $fldstemszu-$fldstemsauf Stem-Structure/s<br>" if ($fldstemsauf>$fldstemszu);
		if ($fldstemsauf != 0 && $fldstemszu != 0) {
			if ((2*@structure)/($fldstemsauf+$fldstemszu)<60) { #The last number determines when the structure is interesting! 60 means at least on Stem in each 60 nt!
				print "<br>The structure seems to contain a lot of secondary structure. If the RNA structure search below<br>does not find a result, it might be interesting to have a closer look at the structures.<br>";
				print 'You might find it useful to look in the book <br><br>"RNA Motifs and Regulatory Elements"<br>Thomas Dandekar (Ed.)<br>Published by Springer<br>ISBN 3-540-41701<br>';
			}
		}
		#An message if it is highly structured that it could be an rRNA
		if (($fldstemauf+$fldstemzu)/2>=10) {
			print "<br>          Highly structured RNA, could this be a ribosomal RNA ?<br>";
		}
		#Okay, that's it for the stem detection!

}

sub predprotein {
	$predprotforAnDom=0;
	print "<b>Coding Sequence:</b><br><b>Pred. Protein<sup>1</sup>:</b>";
	if (@predprot>1){
		print "<br>";	
		for ($count1=1;$count1<=@predprot;$count1++) {
			print $predprot[$count1-1];
			print "<br>" if ($count1%60==0);
		}
	print "<br>";
	$predprotforAnDom=join('',@predprot);
	&runAnDom;
	}
	else {
		print " none";
	}
}

sub runAnDom {
         ##### Creating a unique id using md5sum
                $md5inandom=`md5sum $TEMPDIR/$job.seq`;
                $md5inandom=~/([0-9abcdef]{32})/;
                $md5andom=$1;
        $andomtemp="indexandom$job"."_"."$md5andom".".html";
	`cp /var/www/rnaanalyzer/htdocs/indexandom.html /var/www/rnaanalyzer/tmp/$andomtemp`;
	#`cp /var/www/html/RNAAnalyzer/html/andomstart.html /var/www/html/RNAAnalyzer/results/andomstart$job.html`;
	
	`perl -p -i -e 's/REMOVEMEPLEASE/$predprotforAnDom/g' /var/www/rnaanalyzer/tmp/$andomtemp`;
	#`perl -p -i -e 's/REMOVEME/Main$job/g' /var/www/html/RNAAnalyzer/results/andomstart$job.html`;

	#print '<br><a target="_blank" href="../tmp/'."indexandom$job.html".'">Analyze the Predicted Protein with AnDom</a><br>';
	print '<a target="_blank" href="'."/cgi-bin/get_result_andom.pl?"."$andomtemp".'">Directly Analyze the predicted Protein using AnDom</a><br>';
        #print "<br>Further analyzation with AnDom is temporaryly not available<br>\n";	
}

#okay, now creating the colored RNAStructure for rev2
sub drawcoloredstructure {
	#Here we will create a colored output of the sequence
	#this subroutine requires arrays containing the required information!
	#e.g. @exons for the exons in the format (exonstart,exonend,exonstart,exonend.........)
	# further @transsplicing, @ire, @smsite, @aurichregion, @stemggpairs and further more to be added, see below

	#What to do to invent a new color or feature!
	#1. create an array outside this sub with start and endpoint
	#2. choose in which category it is put. eg FONTCOLOR, UNDERLINED, BOLD/ITALIC
	#3. take care that the start and endpoint are marked in the hash!
	#4. Create an entry below when the nts are counted up, that the feature is writte to the seq-array
	#5. Mark it in the section for the sequence formatting!


		
	#my @seq=split ('',$SEQUENCECHECKED);
	my @seq=@structure;
	#my @exons=(70,80,90,100,160,170);
	#my @transsplicing=(35,58,90,120,125,139, 150,210);
	#my @ire=(20,25);
	#my @smsite=(34,55,66,77);
	#my @aurichregion=(50,60,88,99);
	#my @stemggpairs=(140,160);

	#print ("<br>ire contains: @ire <br>");
    my %change_underline=();      #the changes will be written into these hashes first!
    my %change_fontcolor=();
    my %change_bold=();
    my %change_italic=();
    my @nooutputpossible=();

    my $c=0;
    #Creating arrays to store, where the layout can't be changed further
    my @possible_underline=@seq; #initialize them to be as big as the seq
    my @possible_fontcolor=@seq;
    #my @possible_bold_italic=@seq;
    #my @possible_italic_only=@seq;
    for ($c=0;$c<=@seq-1;$c++) {
        $possible_underline[$c]=0; #and set them to 0 indicating that this field can be altered
        $possible_fontcolor[$c]=0;     #a number above 0 codes for a color or bold or underlined
	#$possible_bold_italic[$c]=0;
    }


    print '<font face="monospace">';

    ########### ATTENTION !!!!! ####################
    ### Old problem the index of array starts at 0 !!! #####
    ### So for the next lines to be correct we will add #####
    ### an x to the beginning of the array !!!!!!!!!! #########
    ### this won't be printed out, but then the ##########
    ### nomenclature is okay !!!!#####################
    @seq=('x',@seq);
    ####

    ###########################################
    ### Codes for the features in the @possible_ arrays ###
    ### transsplicing = 2   #########################
    ### ire = 3                   #########################
    ### smsite=4               #########################
    ### aurichregion=5     #########################
    ### exonboundary=6   #########################
    ### polyasignal=7

    ###################################################################
    ############ Underlined##############################################
    ###################################################################
    my $mark_in_seq_possible=0;
    my $d=0;

    #Processing the trans-splicing hits
    for ($c=0;$c<=@transsplicing-1;$c=$c+2){
        $mark_in_seq_possible=1;
        for ($d=$transsplicing[$c];$d<=$transsplicing[$c+1];$d++) {     #Okay, this checks whether a part of the seq is marked yet
            $mark_in_seq_possible=0 if ($possible_underline[$d]!=0);
        }
        if ($mark_in_seq_possible==1) {
            #1. mark all fields in the array as nomorepossible
            for ($d=$transsplicing[$c];$d<=$transsplicing[$c+1];$d++) {
                $possible_underline[$d] = 1;
            }
            #Add the lines to the hash
            $change_underline{$transsplicing[$c]}='transstart';
            $change_underline{$transsplicing[$c+1]}='transend';
        }
        else { #If it is already marked
            push (@nooutputpossible,"Trans-splicing hit from $transsplicing[$c] to $transsplicing[$c+1]");
        }
    }

    #Processing the IRE hits
    for ($c=0;$c<=@ire-1;$c=$c+2){
        $mark_in_seq_possible=1;
        for ($d=$ire[$c];$d<=$ire[$c+1];$d++) {     #Okay, this checks whether a part of the seq is marked yet
            $mark_in_seq_possible=0 if ($possible_underline[$d]!=0);
        }
        if ($mark_in_seq_possible==1) {
            #1. mark all fields in the array as nomorepossible
            for ($d=$ire[$c];$d<=$ire[$c+1];$d++) {
                $possible_underline[$d] = 1;
            }
            #Add the lines to the hash
            $change_underline{$ire[$c]}='irestart';
            $change_underline{$ire[$c+1]}='ireend';
        }
        else { #If it is already marked
            push (@nooutputpossible,"IRE hit from $ire[$c] to $ire[$c+1]");
        }
    }

    ####################################################################################
    ########################## Font Color #################################################
    ####################################################################################
   
    #Processing Promotors
    for ($c=0;$c<=@promotor-1;$c=$c+2){
            $mark_in_seq_possible=1;
            for ($d=$promotor[$c];$d<=$promotor[$c+1];$d++) {     #Okay, this checks whether a part of the seq is marked yet
	             $mark_in_seq_possible=0 if ($possible_fontcolor[$d]!=0);
 	    }
	    if ($mark_in_seq_possible==1) {
        	#1. mark all fields in the array as nomorepossible
                for ($d=$promotor[$c];$d<=$promotor[$c+1];$d++) {
     	           $possible_fontcolor[$d] = 8; #8 means purple
                }
                $change_fontcolor{$promotor[$c]}='promotorstart';
                $change_fontcolor{$promotor[$c+1]}='promotorend';
            }
            else { #If it is already marked
                push (@nooutputpossible,"Promotor from $promotor[$c] to $promotor[$c+1]");
            }
    }

    
    #Processing the PolyASignal
    for ($c=0;$c<=@polyasignal-1;$c=$c+2){
        $mark_in_seq_possible=1;
        for ($d=$polyasignal[$c];$d<=$polyasignal[$c+1];$d++) {     #Okay, this checks whether a part of the seq is marked yet
            $mark_in_seq_possible=0 if ($possible_fontcolor[$d]!=0);
        }
        if ($mark_in_seq_possible==1) {
            #1. mark all fields in the array as nomorepossible
            for ($d=$polyasignal[$c];$d<=$polyasignal[$c+1];$d++) {
                $possible_fontcolor[$d] = 7;
            }
            $change_fontcolor{$polyasignal[$c]}='polyasignalstart';
            $change_fontcolor{$polyasignal[$c+1]}='polyasignalend';
        }
        else { #If it is already marked
            push (@nooutputpossible,"PolyA-Signal from $polyasignal[$c] to $polyasignal[$c+1]");
        }
    }
    #Processing the smsites
    for ($c=0;$c<=@smsite-1;$c=$c+2){
        $mark_in_seq_possible=1;
        for ($d=$smsite[$c];$d<=$smsite[$c+1];$d++) {     #Okay, this checks whether a part of the seq is marked yet
            $mark_in_seq_possible=0 if ($possible_fontcolor[$d]!=0);
        }
        if ($mark_in_seq_possible==1) {
            #1. mark all fields in the array as nomorepossible
            for ($d=$smsite[$c];$d<=$smsite[$c+1];$d++) {
                $possible_fontcolor[$d] = 5;
            }
            $change_fontcolor{$smsite[$c]}='smsitestart';
            $change_fontcolor{$smsite[$c+1]}='smsiteend';
        }
        else { #If it is already marked
            push (@nooutputpossible,"Sm-Site from $smsite[$c] to $smsite[$c+1]");
        }
    }

    #Processing the au-rich-regions
    for ($c=0;$c<=@aurichregion-1;$c=$c+2){
        $mark_in_seq_possible=1;
        for ($d=$aurichregion[$c];$d<=$aurichregion[$c+1];$d++) {     #Okay, this checks whether a part of the seq is marked yet
            $mark_in_seq_possible=0 if ($possible_fontcolor[$d]!=0);
        }
        if ($mark_in_seq_possible==1) {
            #1. mark all fields in the array as nomorepossible
            for ($d=$aurichregion[$c];$d<=$aurichregion[$c+1];$d++) {
                $possible_fontcolor[$d] = 5;
            }
            #Add the lines to the hash
            $change_fontcolor{$aurichregion[$c]}='aurichregionstart';
            $change_fontcolor{$aurichregion[$c+1]}='aurichregionend';
        }
        else { #If it is already marked
            push (@nooutputpossible,"AU-rich region from $aurichregion[$c] to $aurichregion[$c+1]");
        }
    }
    #Processing stemggpairs
    for ($c=0;$c<=@stemggpairs-1;$c=$c+2){
        $mark_in_seq_possible=1;
        for ($d=$stemggpairs[$c];$d<=$stemggpairs[$c+1];$d++) {     #Okay, this checks whether a part of the seq is marked yet
            $mark_in_seq_possible=0 if ($possible_fontcolor[$d]!=0);
        }
        if ($mark_in_seq_possible==1) {
            #1. mark all fields in the array as nomorepossible
            for ($d=$stemggpairs[$c];$d<=$stemggpairs[$c+1];$d++) {
                $possible_fontcolor[$d] = 4;
            }
            #Add the lines to the hash
            $change_fontcolor{$stemggpairs[$c]}='stemggpairsstemstart';
            $change_fontcolor{$stemggpairs[$c]+1}='stemggpairsstemend';
            $change_fontcolor{$stemggpairs[$c]+2}='stemggpairsbulge';
            $change_fontcolor{$stemggpairs[$c]+3}='stemggpairsstemstart';
            
            $change_fontcolor{$stemggpairs[$c+1]-3}='stemggpairsstemend';
            $change_fontcolor{$stemggpairs[$c+1]-2}='stemggpairsbulge';
            $change_fontcolor{$stemggpairs[$c+1]-1}='stemggpairsstemstart';
            $change_fontcolor{$stemggpairs[$c+1]}='stemggpairsstemend';
        }
        else { #If it is already marked
            push (@nooutputpossible,"Stem-GG-Pair from $stemggpairs[$c] to $stemggpairs[$c+1]");
        }
    }
     

    ############################################################################
    ########################### Bold ###########################################
    ############################################################################
    #Processing the exons
    for ($c=0;$c<=@exons-1;$c=$c+2){
        #Add the lines to the hash
        $change_bold{$exons[$c]}='exonstart';
        $change_bold{$exons[$c+1]}='exonend';
    }

    ############################################################################
    ########################### Italic #########################################
    ############################################################################
    #Processing the untranslated regions UTR
    for ($c=0;$c<=@utr-1;$c=$c+2){
        #Add the lines to the hash
        $change_italic{$utr[$c]}='utrstart';
        $change_italic{$utr[$c+1]}='utrend';
    }

    #print $rv->h3("");
    #print $rv->h1("Summary");
    
    #Print out the predicted protein! reformat for a nicer looking on the page	
    #&predprotein;


    #Now create the html output!!!!!!!!!!!!

    print "<br><b>Colored Structure:</b><br>";

    my $ntausgabe=1;

    my $boldonout=0;
    my $fontcolor=0; #0=black 1=red 2=blue 3=navy 4=fuchsia
    my $underlined=0; #0=no 1=yes

    for ($c=1;$c<=@seq-1;$c++){
        ###########################    
        ### OUTPUT UNDERLINED ######
        ###########################
        if (defined $change_underline{$c}) {    
            if ($change_underline{$c} eq 'transstart' || $change_underline{$c} eq 'irestart')  {
                $seq[$c]='<u>'.$seq[$c];
                $underlined=1;
            }
            if ($change_underline{$c} eq 'transend' || $change_underline{$c} eq 'ireend') {
                $seq[$c]=$seq[$c].'</u>';
                $underlined=0;
            }
        }       
        
        ###########################
        ###### OUTPUT BOLD #########
        ###########################
        if (defined $change_bold{$c}) {
            if ($change_bold{$c} eq 'exonstart') {
                $seq[$c]='<b>'.$seq[$c];
                $boldonout=1;
            }
            if ($change_bold{$c} eq 'exonend') {
                $seq[$c]=$seq[$c].'</b>';
                $boldonout=0;
            }
        }
        
       
	###########################
	###### OUTPUT ITALIC ######
	###########################
	if (defined $change_italic{$c}) {
             if ($change_italic{$c} eq 'utrstart') {
                     $seq[$c]='<i>'.$seq[$c];
                     $italiconout=1;
             }
             if ($change_italic{$c} eq 'utrend') {
                     $seq[$c]=$seq[$c].'</i>';
                     $italiconout=0;
             }
        }
	
        ############################
        ####### OUTPUT COLOR ########
        ############################
        if (defined $change_fontcolor{$c}) {
		

       	   ######## Promotor ##########
           if ($change_fontcolor{$c} eq 'promotorstart') {
                 $seq[$c]='<font color=purple>'.$seq[$c];
                 $fontcolor=8;
           }
           if ($change_fontcolor{$c} eq 'promotorend') {
                $seq[$c]=$seq[$c].'</font>';
                $fontcolor=0;
           }
           ######## PolyA-Signal ##########
           if ($change_fontcolor{$c} eq 'polyasignalstart') {
               $seq[$c]='<font color=lime>'.$seq[$c];
               $fontcolor=5;
           }
           if ($change_fontcolor{$c} eq 'polyasignalend') {
                 $seq[$c]=$seq[$c].'</font>';
                 $fontcolor=0;
           }
	   ######## Sm-Site ##########
           if ($change_fontcolor{$c} eq 'smsitestart') {
                $seq[$c]='<font color=red>'.$seq[$c];
                $fontcolor=1;
            }
            if ($change_fontcolor{$c} eq 'smsiteend') {
                $seq[$c]=$seq[$c].'</font>';
                $fontcolor=0;
            }
            ####### AU-rich region #######
            if ($change_fontcolor{$c} eq 'aurichregionstart') {
                $seq[$c]='<font color=blue>'.$seq[$c];
                $fontcolor=2;
            }
            if ($change_fontcolor{$c} eq 'aurichregionend') {
                $seq[$c]=$seq[$c].'</font>';
                $fontcolor=0;
            }
            ####### Stem-GG-Pairs #######
            if ($change_fontcolor{$c} eq 'stemggpairsstemstart') {
                $seq[$c]='<font color=fuchsia>'.$seq[$c];
                $fontcolor=4;  
            }
            if ($change_fontcolor{$c} eq 'stemggpairsstemend') {
                $seq[$c]=$seq[$c].'</font>';
                $fontcolor=0;
            }
            if ($change_fontcolor{$c} eq 'stemggpairsbulge') {
                $seq[$c]='<font color=green>'.$seq[$c].'</font>';
            }
        }
        
     
        #############################
        ###Formatting the Output ########
        ### into suitable lines ###########
        #############################
        
        
        if ($c % 50 == 0) {   
            if ($boldonout==1) {
                $seq[$c]=$seq[$c].'</b>';
                $seq[$c+1]='<b>'.$seq[$c+1];
            }
	    if ($italiconout==1) {
		$seq[$c]=$seq[$c].'</i>';
		$seq[$c+1]='<i>'.$seq[$c+1]; 
	    }
            if ($underlined==1) {
                $seq[$c]=$seq[$c].'</u>';
                $seq[$c+1]='<u>'.$seq[$c+1];
            }
            if ($fontcolor != 0) {
                $seq[$c]=$seq[$c].'</font>';
                $seq[$c+1]='<font color = red>'.$seq[$c+1] if ($fontcolor==1);
                $seq[$c+1]='<font color = blue>'.$seq[$c+1] if ($fontcolor==2);
                $seq[$c+1]='<font color = fuchsia>'.$seq[$c+1] if ($fontcolor==4);
		$seq[$c+1]='<font color = lime>'.$seq[$c+1] if ($fontcolor==5);
		$seq[$c+1]='<font color = purple>'.$seq[$c+1] if ($fontcolor==8);
            }
        }
    }

    #for debugging only
    #$seq[30]='<span style="underline-color: yellow">'.$seq[30];
    #$seq[45]=$seq[45].'</span>';
    ################


    for ($c=1;$c<=@seq-1;$c++){
        print "$seq[$c]";
        if ($c % 50 == 0) {
            print "     ".$c."<br>";
        }
        if ($c %10 == 0) {
            print "";
        }
    }
    #print;
    @seq=();
    print "<br>Output not possible:<br>@nooutputpossible<br>" if (@nooutputpossible>0);
    #print "<br><br>Legend:<br><b>BOLD</b> marks EXONS<br><u>UNDERLINED</u> marks IRE or TRANSSPLICNG hits<br>";
    #print "<i>ITALIC</i> marks putative UTRs<br>";
    #print "<font color=red>RED </font> marks SMSITES or snRNP-binding motifs<br><font color=blue>BLUE</font> marks AU-rich regions<br>";
    #print "<font color=fuchsia>FUCHSIA</font> with <font color=green>GREEN</font> marks Stemm-GG-Pairs<br>";
    #print "<font color=lime>LIME</font> marks PolyA-signal<br><font color=purple>PURPLE</font> marks a Promotor<br>";
    #print "<br>Pr.A1.bin.site = Protein A1 binding site<br>";
    
    #print "<br><br>%change_fontcolor<br><br>";
}


>>>>>>> origin/master
