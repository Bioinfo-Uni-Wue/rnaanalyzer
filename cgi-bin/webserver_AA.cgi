#!/usr/bin/perl

use lib ".";
use lib "./RNASERVER/";
use CGI;
use RNASERVER::TRANS2;
use RNASERVER::IRE;
use Bio::Tools::Genscan;
use Cwd;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use RNASERVER::JobUtil qw(get_next_job_id);
use JSON;
use File::Slurp;
use File::Path qw(make_path);
use Cwd qw(abs_path);
use CGI qw(:standard escapeHTML);
use Bio::Seq;


$debug=0;

#Dir-Localisations
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
$MAXFOLDINGLEN=5000;
$MAXFOLDINGLENUTR=5000;
$MAXFORNALENGTH=5000;
$maxcoloredseqlen=10000;

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
my $dnarna              = $params->{dnarna};
my $SEQNAMECHECKED      = $params->{sequence_name};
my $SEQUENCECHECKED     = $params->{sequence_clean};
my $SEQUENCELENGTH      = $params->{sequence_length};

my $html_file = "$TEMPDIR/result.html";

open my $out, ">", $html_file or die "Can't write result page: $!";
select $out;  # Redirect STDOUT to file

# We print into one html which can be accessed again or bookmarked 
# print "<!DOCTYPE html>\n";
# print "<html><head><meta charset='UTF-8'><title>Job $job Results</title></head><body>";
print <<'HTML';
<html><head>
  <title>Batch Results</title>
  <link rel="stylesheet" href="/css/newresults.css">
</head><body>
<header>
    <a href="http://localhost">    <!--- change after putting on server -->
      <img src="http://localhost/images/logo.png" alt="RNA Analyzer Logo" class="logo" />
    </a>
    <div class="header-text">
      <h1>RNA Analyzer 2025</h1>
      <p>Webserver for RNA Sequence Overview</p>
    </div>
    <div class="header-links">
      <a href="http://localhost/about.html" target="_blank">About</a> |
      <a href="http://localhost/contact.html" target="_blank">Contact</a> |
      <a href="https://www.biozentrum.uni-wuerzburg.de/bioinfo" target="_blank">Dandekar Lab</a>
    </div>
  </header>

<main>
<div class="container">
HTML
# ... print your actual results here ...


    
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

    
    print $cgi->h2("Here are the results for JOB ID: $job with sequence name: ". CGI::escapeHTML($SEQNAMECHECKED));

    # running analysis
    &analysis;
}

##calls all the new subrotines
sub analysis {
    chdir $TEMPDIR;

    print "<pre>";

    print "<div class='box'>";
    print "<div class='box-header'>Structural information</div>";
    print "<div class='box-content'>";
    if (length $SEQUENCECHECKED <= $maxcoloredseqlen) {
        &createfolding;
        &checkstems;
        &stemggpairs;
    } else {
        print "<br><b>Length:</b>\t$SEQUENCELENGTH";
        print "     *some information is only available up to $MAXFOLDINGLEN nt\n" if ($SEQUENCELENGTH > $MAXFOLDINGLEN);
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

    # print "<div class='box'>";
    # &csfce;
    # print "</div>";

    print "<div class='box'>";
    &createfoldingpicture;
    print "</div>";

    print "<div class='table-container'>";
    &location_table;
    print "</div>";

    # &drawcoloredsequence;

    print "</pre>";
}

sub TRANS{
	#this checkes the ciona-consensus !		
			@transcionareturnvalues=RNASERVER::TRANS2::ciona($SEQUENCECHECKED); ## Ist das der CIONA ??? Auf jeden Fall aber SCHISOTSOMA
			print "<div class='box-header'>Trans-Splicing:</div>";
			#print "<b> <big>Putative trans-splicing Schistosoma-consensus</b></big> search:";
            print "<div class='box-content'>";
			if (@transcionareturnvalues==1) {
				#print "No hit detected<br>";
				print "<h3><b>Schistosoma:</b>\tNone</h3><br>";
			}
			else {
				$hits=pop @transcionareturnvalues;
				for ($count=0;$count<$hits;$count++){
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
				}
			}
			#this checkes the c. elegans consensus !!
			@transcelegansvalues=RNASERVER::TRANS2::celegans($SEQUENCECHECKED);	
			if (@transcelegansvalues==1){
				#print "No hit detected<br>";
				print "<h3><b>C. elegans:</b>\tNone</h3>";
			}
			else {
				$hits=pop @transcelegansvalues;
				for ($count=0;$count<$hits;$count++){
					$transcelegansvalues[$count*10+5]=uc($transcelegansvalues[$count*10+5]);
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
    print "<div class='box-header'>Iron-resp Elements</div>";
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
            
            print "</table>\n";
            print "<br>\n"; 
        }
        
        $irelineprintout = 1;
        
    } else {
        print "<h3><b>Iron-resp Ele.:</b> None</h3>";
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
    print "<div class='box-content'>";

    ELEMENT1: for ($count1=0;$count1<@seq-14;$count1++) {
        my $mismatch=0;
        for ($count2=0;$count2<=14;$count2++){
            $mismatch ++ if ($seq[$count1+$count2] ne $element1[$count2]);
            next ELEMENT1 if ($mismatch >=2);
        }
        if ($putativeCVfound == 0) {
            print "<h4>CstF Motif Hits</h4>\n";
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

    print "</div>";
}

sub stemggpairs {
    # print "<div class='box-content'>";
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
                pos($str)=pos($str)-2;
            }        
        }
        pos($str)=$anfang-2;
    }
    if ($stemggleadlineprinted == 1) {
                        print "</table>\n";
                    }

    print "<h3>No Stem GG Pair Found.</h3><br>" if ($stemggpairfound==0); 
    @sequ=();
    $str='';

    # print "</div>";
}


sub rbp {
    print "<div class='box-header'>RNA Binding Protein Scan:</div>";
    print "<div class='box-content'>";

    my $fimo_outfile = "$TEMPDIR/fimo.txt";
    my $input_seq = "$TEMPDIR/$job.seq";  # path to your .meme file

    # Run FIMO
    my $fimo_cmd = "$FIMO/fimo --text --thresh 1e-5 $RBPDB $input_seq > $fimo_outfile";
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

            push @rbp_locs, $start, $stop;

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

    # Display results
    my $total = scalar @fimo_results;

    if ($total > 0) {
        print "<table class='table-result'>\n";
        print "<tr><th>Protein</th><th>Motif</th><th>Start</th><th>End</th><th>Score</th><th>p-value</th><th>Matched Sequence</th></tr>\n";

        foreach my $hit (@fimo_results) {
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
    print "<div class='box-header'>Au-rich regions:</div>";

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
        print "<h3>None found.</h3>     *(AU-rich region of at least 30 nt)<br>";
    }
    print "</div>";
}


sub tRNA {
    # Looking for tRNAs using tRNAscan-SE

    print "<div class='box-header'>tRNA Scan:</div>";
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
        print "<h3>None Found.</h3>";
    }
    
    print "</div>";
    # print "DEBUG: @trna_loc" if @trna_loc;
}


sub smsite {
    print "<div class='box-header'>Catalytic RNA:</div>";
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
    
    print "<h3>No snRNP-motifs found.</h3>" if ($leadlineprinted==0);
    
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
        print "\t\tHighly structured RNA can be likely rRNA, tRNA, or any other regulatory RNA\n";
    } elsif ($hairpin_count >= 1 && $stem_count <= 3) {
        print "\t\tSimple structured RNA, possible miRNA, siRNA, or regulatory element\n";
    } elsif ($stem_count >= 1) {
        print "\t\tSome secondary structure detected, may have biological significance\n";
    } else {
        print "\t\tMinimal secondary structure\n";
    }
    
    print "\t\t****It might be interesting to have a closer look at the structures.\n";
    print "\t\tYou might find it useful to look in the book\n";
    print "\t\t'RNA Motifs and Regulatory Elements'\n";
    print "\t\tThomas Dandekar (Ed.)\n";
    print "\t\tISBN 3-540-41701\n";

    if ($hairpin_count > 0 && $stem_count > 0 && $hairpin_count / $stem_count > 0.7) {
        print "\t\tHigh hairpin content — characteristic of miRNA precursors\n";
    }

    print "<br>";
    
    return ($stem_count, $hairpin_count);
}


sub checkstemsonly {
    # revamped checkstemsonly based on the logic of the new checkstems
    my $inseq = $_[0];
    my $inwhattodo = $_[1];   # 0 means is structure, evaluate, 1--> sequence, please fold first!
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
        next if $pairing[$i] < $i;  # avoid counting the same pair twice
        
        my $left = $i;
        my $right = $pairing[$i];
        
        # BIOLOGICAL VALIDATION: Check minimum loop length
        my $loop_length = $right - $left - 1;
        next if $loop_length < 3;  # Skip sterically impossible loops
        
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
	$predprotforAnDom=0;
	print "<div class='box-header'>Predicted Protein:</div>";
    print "<div class='box-content'>";
    if (defined $cpc){
        print "<i>Protein is predicted from CPC2 output.<i><br>";
    }
	if (@predprot>1){
		print "<br>";	
		for ($count1=1;$count1<=@predprot;$count1++) {
			print $predprot[$count1-1];
			# print "<br>" if ($count1%120==0);
		}
	print "<br>";
	}
	else {
		print " none";
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

	print "<div class='box-header'>RNA motif search:</div>";
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
        $forna_html .= "<div id='rna_ss_0' style='width: 250px; height: 250px; display: none; 10px; overflow: hidden;'></div>\n";
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
		print "No RNA motif recognized\n";
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

	print "<div class='box-header'>Coding potential:</div>";
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
        print "No result in CPC2 output\n";
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
        print "<i>Transcript predicted to be noncoding. Running structural region scan instead.</i><br>";
        print "<div class='box'>";
        &scan_structured_regions;
        "</div>";
    }
    print "</div>";
    return \%transcripts;
    
}

sub microRNA {

		my $mirbase_output = "$TEMPDIR/$job.mirtbl";
		my $mirbase_out    = "$TEMPDIR/$job.mir";

		my $mirna_search = "$HMMER/nhmmer --rna --watson -Z 3.73 -E 1 --tblout $mirbase_output -o $mirbase_out $TEMPDIR/$job.seq $MIRBASE/hairpin.fa";

		system($mirna_search);

        print "<div class='box-header'>miRNA scan:</div>";
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
			print "<b>Total microRNA hits found:</b> $total\n";
		} else {
			print "No regions matching a mircroRNA was found.\n";
		}

		print "</div>";


}

# microRNA scan using miRanda
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


    print "<div class='box-header'>miRNA target prediction</div>\n";
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


    print "<i>Extracted $found_valid valid 3' UTR regions for target prediction.</i><br>";

    # no valid UTRs 
    if (!$found_valid) {
        print "<i>No valid 3' UTRs found. Scanning full sequence instead.</i><br>";
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

    print "<div class='box-header'>Gene Prediction:</div>";
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
        print "No gene predictions were found.<br>";
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

    print "<div class='box-header'>UTR(s) Prediction:</div>";
    print "<div class='box-content'>";

    print "<i>UTR prediction source: $source</i><br>";

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

    print "<b>UTR:</b>           start  -   end   -  stems - energy<br>";
    for (my $i = 0; $i < @utrprintout; $i += 3) {
        my ($type, $start, $end) = @utrprintout[$i, $i+1, $i+2];

        next if $end < $start || $start < 1 || $end > $seq_length;

        my $utr_seq = substr($seq, $start - 1, $end - $start + 1);
        next unless length($utr_seq) > 0;

        my @returnout = checkstemsonly($utr_seq, 1);

        printf(" %d'            %-6d - %6d", $type, $start, $end);
        print "       $returnout[0]-$returnout[1]" if $returnout[0] != $returnout[1];
        print "       $returnout[0]" if $returnout[0] == $returnout[1];
        print "     $returnout[2]<br>" if $returnout[2] != 1;

        if ($type == 5) {
            print "         SD motif<br>" if $utr_seq =~ /AGGAGG/i;
            print "         Kozak motif<br>" if $utr_seq =~ /gcc[AG]ccATGG/i;
        }

        if ($type == 3) {
            foreach my $motif (qw(AATAAA ATTAAA TATAAA AAGAAA AGTAAA AATATA)) {
                if ($utr_seq =~ /$motif/i) {
                    my $pos = $-[0] + $start;
                    print "         PolyA signal $motif at $pos<br>";
                    my $signal_end = $pos + length($motif) - 1;
                    push @polyasignal, $pos, $signal_end;
                    last;
                }
            }
            if ($utr_seq =~ /A{10,}/i) {
                my $tail_pos = $-[0] + $start;
                print "         PolyA tail near $tail_pos<br>";
                my $tail_end = $+[0] - 1 + $start;
                push @polyatail, $tail_pos, $tail_end;
            }
        }
    }

    print "<i>No valid UTRs inferred.</i><br>" unless @utr;

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
            print "PolyA signal ($motif) detected at $signal_pos<br>";
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
            print "PolyA tail detected near $tail_pos<br>";
            last;
        }
    }

    if ($signal_pos > 0 || $tail_pos > 0) {
        my $utr_start = $cds_end + 1;
        my $utr_end   = $tail_pos > 0 ? $tail_pos : ($signal_pos + 20);

        push @new3primeutr, $utr_start, $utr_end;
        push @utrprintout, 3, $utr_start, $utr_end;
        push @utr, $utr_start, $utr_end;

        print "Inferred 3' UTR based on polyA: $utr_start - $utr_end<br>";
    } else {
        print "No strong polyA signal/tail detected in 3' region.<br>";
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
        print "No Shine-Dalgarno motifs detected in 5′ UTR regions.<br>";
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
    print "<br><br>";
    print "<div class='box-header'>Structured region scan</div>";
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
        print "<br>** Highly structured regions found. Consider tRNA, rRNA, or ncRNA elements.<br>";
    } else {
        print "<i>No regions with significant RNA structure detected.</i><br>";
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

	# if ($SEQUENCELENGTH <= $MAXFOLDINGLEN && $SEQUENCELENGTH > $MAXFORNALENGTH) {
    #     my $svg_file = "$TEMPDIR/${SEQNAMECHECKED}_ss.svg";
    #     my $ps_url  = "/tmp/jobs/job_$job/${SEQNAMECHECKED}_ss.ps";
    #     my $svg_url = "/tmp/jobs/job_$job/${SEQNAMECHECKED}_ss.svg";

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
    #     print "<b>Download As: </b>\n";
    #     print "<a href='$svg_url' target='_blank'><button>SVG File</button></a>";
    #     print "<a href='$ps_url' target='_blank'><button>PS File</button></a>\n";

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

        # Make sure RNAplot output is ready
        system("$VIENNARNAFOLDDIR/RNAplot --infile=$TEMPDIR/$job.foldout -f svg --filename-full");

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

        # Print HTML content
        print "<div class='box-header'>RNA Structure Visualization:</div>";
        print "<div class='box-content'>";

        if ($SEQUENCELENGTH >= 2500){
            print "Structure visualization can take a few seconds to load due to sequence length.";
        }
        
        print "<div id='rna_ss'>RNA Structure</div>";
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

        for (my $j = 0; $j < @transsplicing; $j += 2) {
            my ($from, $to) = @transsplicing[$j, $j + 1];

            ($from, $to) = ($to, $from) if $from > $to;

            for (my $i = $from; $i <= $to; $i++) {
                $color_text .= "$i:cyan ";
            }
        }

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
        print "  var container;\n";  # Make container accessible globally\n";

        print "  window.onload = function () {\n";
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
        print "  };\n";
        print "</script>";

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
        print "<span style='display: inline-block; width: 12px; height: 12px; background: cyan; margin-left: 10px;'></span> TRANS-splicing\t";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: lime; margin-left: 10px;'></span> PolyA Signal/PolyA Tail\t";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: pink; margin-left: 10px;'></span> Protein Binding Site(s)";

        print "</div>";
        print "</div>";

    } else {
        print "<br><b>Maximum folding limit reached</b><br>";
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

    print "<h2>Locations of the Detected Structures:</h2>\n";
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

write_file("$TEMPDIR/result.txt", "done\n");

print "</div></main></body></html>";
close $out;