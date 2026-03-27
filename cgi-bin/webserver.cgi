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
use List::Util qw(min);
use RNASERVER::SanitizedOut qw(sanitized_output_link);
use RNASERVER::FoldDiagram qw(
    render_raw_fold_svg
    render_annotated_fold_svg
);

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
$YAMTK=abs_path('../bin/yamtk'); #yamtk scan folder
$RBPDB=abs_path('../databases/rbpdb/motifs.meme');
$RIBOSWDB=abs_path('../databases/riboswitches/riboswitch.cm'); #riboswitch location 
$BINDIR=abs_path('../cgi-bin'); #bin directory
$JAVA=abs_path('../bin/jdk-11/bin/java'); #java location
$VARNA_JAR= abs_path('../bin/jdk-11/varna/VARNAv3-93.jar'); #varna location
$IPKNOT=abs_path('../bin/ipknot/bin'); #IPKnot location
$MAXFOLDINGLEN=5000;
$MAXFOLDINGLENUTR=5000;
$MAXFORNALENGTH=5000;

# aprsing cgi
my $cgi = CGI->new;


$job_id = $ARGV[0] // $cgi->param("job_id");  # passed from batch script
$job = $job_id;


$TEMPDIR = abs_path("../tmp/jobs/job_$job");
$download_dir = ("/tmp/jobs/job_$job");

# a directory to store sanitzed out for users
$raw_out = ("$TEMPDIR/raw_out");
mkdir $raw_out unless -d $raw_out;

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
my $do_centroid         = $params->{centroid};
my $do_annotated        = $params->{annotated};
my $do_fornac            = $params->{fornac};
my $do_pseudoknot         = $params->{pseudoknot};
my $do_subopt             = $params->{subopt};
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
print "    <a href='localhost/'>";    # change after putting on server
print "      <img src='localhost//images/logo.png' alt='RNA Analyzer Logo' class='logo' />";
print "    </a>";
print "    <div class='header-text'>";
print "      <h1>RNA Analyzer<sup>3</sup></h1>";
print "      <p>Webserver for RNA Sequence Overview</p>";
print "    </div>";
print "    <div class='header-links'>";
print "      <a href='localhost/about.html' target='_blank'>About</a> |";
print "      <a href='localhost/contact.html' target='_blank'>Contact</a> |";
print "      <a href='' target='_blank'>Dandekar Lab</a>";
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
            print "<i>\ti. Annotated Structure visualization can take a few seconds to load due to sequence length. Please wait till the page loads.</i>\n";
            print "</div>";
        }
    
    print "<br>";
    
        # running analysis
    &analysis;
}

##calls all the new subrotines
sub analysis {
    chdir $TEMPDIR;

    # creates folding and does structural analysis like stems etc
    print "<div class='box'>";
    print "<div class='box-header' onclick='toggleBoxContent(this)'>Structural information</div>";
    print "<div class='box-content-structure'>";
    if (length $SEQUENCECHECKED <= $MAXFOLDINGLEN) {
        &createfolding;
        if ($do_subopt) {
            &subopt;
        }
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

    # base composition analysis 
    &composition;

    # trans-splicing analysis
    if ($do_TRANS) {
    print "<div class='box'>";
    &TRANS;
    print "</div>";
    }

    # iron resposive element check
    if ($do_IRE) {
    print "<div class='box'>";	
	&IRE;
    print "</div>";
    }

    # check for riboswitches
    if ($do_ribo) {
    print "<div class='box'>";	
	&riboswitch;
    print "</div>";
    }

    # check for Au rich regions
    print "<div class='box'>";
    &ARE;
    print "</div>";

    # check for catalytic RNA sites, including snRNP sites and sm Sites
    print "<div class='box'>";
    &smsite;
    print "</div>";

    # RNA motif scan agaisnt Rfam database  
    if ($do_rnamotif) {
    print "<div class='box'>";
    &RNAMOTIF;
    print "</div>";
    }

    # tRNA scan using tRNAscan-SE
    if ($do_trna) {
    print "<div class='box'>";
    &tRNA;
    print "</div>";
    }

    # RNA binding Protein motif scan
    if ($do_rbp) {
    print "<div class='box'>"; 
    &rbp;
    print "</div>";
    }

    # microRNA scan 
    if ($do_mirna) {
    print "<div class='box'>";
    &microRNA;
    print "</div>";
    }

    # coding potential analysis and UTR prediction
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

    # outputs protein seq when coding transcript
    print "<div class='box'>";
    &predprotein;
    print "</div>";

    # miRNA target scan using miranda
    if ($do_mirnatarget) {
    print "<div class='box'>";
    miRNAtarget(\%transcripts);
    print "</div>";
    } 
    # Optional: enable if scan_rbs supports multi-transcript
    # scan_rbs(\%transcripts);

    # Cstf motif scan
    print "<div class='box'>";
    &csfce;
    print "</div>";

    # creates annotated folding picture using Varna
    if (!($do_annotated && $do_fornac)){
    print "<div class='box'>";
    &createfoldingpicture;
    print "</div>";
    }

    if ($do_annotated && $do_fornac) {
    print "<div class='box'>";
    &createfoldingpictureFornac;
    print "</div>";
    }

    # feature table
    print "<div class='table-container'>";
    &location_table;
    print "</div>";

    &make_txt;
    # &drawcoloredsequence;

    # print "</pre>";
}

sub TRANS{
	#this checkes the ciona-consensus !		
			@transcionareturnvalues=RNASERVER::TRANS2::ciona($SEQUENCECHECKED); ## Ist das der CIONA ??? Auf jeden Fall aber SCHISOTSOMA
			print "<div class='box-header' onclick='toggleBoxContent(this)'>Trans-Splicing Analysis:<a href='localhost/TRANSinfo.html' target='_blank' onclick='event.stopPropagation();' style='display:inline-block; margin-left:6px; width:16px; height:16px; line-height:16px; text-align:center; border-radius:50%; background:#ccc; color:#000; font-size:14px; text-decoration:none;' title='Help'>?</a></div>";
			#print "<b> <big>Putative trans-splicing Schistosoma-consensus</b></big> search:";
            print "<div class='box-content'>";
			if (@transcionareturnvalues==1) {
				print "<div class='info-warning'>";
				print "<b>Schistosoma:</b>\tNo Trans-splicing element(s) found.<br>";
                print "</div>";
			}
			else {
				$hits=pop @transcionareturnvalues;
				for ($count=0;$count<$hits;$count++){

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
                print "<div class='info-warning'>";
				print "<b>C. elegans:</b>\tNo Trans-splicing element(s) found.";
                print "</div>";
			}
			else {
				$hits=pop @transcelegansvalues;
				for ($count=0;$count<$hits;$count++){
					$transcelegansvalues[$count*10+5]=uc($transcelegansvalues[$count*10+5]);

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
    my @ire_return_values  = RNASERVER::IRE::suboptimalfindire($SEQUENCECHECKED);
    my @epas1_return_values = RNASERVER::IRE::find_epas1_ire($SEQUENCECHECKED);

    my @all_ire_hits = (@ire_return_values, @epas1_return_values);

    $irelineprintout = 0;

    print "<link rel='stylesheet' href='/css/fornac.css'>";
    print "<script src='https://d3js.org/d3.v3.min.js'></script>";
    print "<script src='/js/fornac.js'></script>";

    # print "<div class='box-header' onclick='toggleBoxContent(this)'>Iron-resp Element(s):</div>";
    print "<div class='box-header' onclick='toggleBoxContent(this)'>Iron-resp Element(s): <a href='localhost/IREinfo.html' target='_blank' onclick='event.stopPropagation();' style='display:inline-block; margin-left:6px; width:16px; height:16px; line-height:16px; text-align:center; border-radius:50%; background:#ccc; color:#000; font-size:14px; text-decoration:none;' title='Help'>?</a></div>";
    print "<div class='box-content'>";

    if (@all_ire_hits > 1) {

        my %positions;

        for (my $count = 0; $count <= $#all_ire_hits; $count += 8) {
            my $pos = $all_ire_hits[$count];

            push @{ $positions{$pos} }, {
                position  => $all_ire_hits[$count],
                quality   => $all_ire_hits[$count + 1],
                sequence  => $all_ire_hits[$count + 2],
                structure => $all_ire_hits[$count + 3],
                energy    => $all_ire_hits[$count + 4],
                upper     => $all_ire_hits[$count + 5],
                lower     => $all_ire_hits[$count + 6],
                label     => $all_ire_hits[$count + 7],
            };
        }

        my $total_positions  = scalar keys %positions;
        my $total_structures = 0;
        $total_structures += @{ $positions{$_} } for keys %positions;

        print "<div class='info-info'>";
        print "<strong>Found:</strong> $total_structures structure(s) at $total_positions position(s)";
        print "</div>";

        print "<table class='table-result'>\n";
        print "<tr><th>Position</th><th>Sequence</th><th>Label</th><th>Structure</th><th>Energy</th><th>Quality</th><th class='no-print'>View Structure</th></tr>\n";

        my $structure_counter = 0;

        foreach my $pos (sort { $a <=> $b } keys %positions) {
            my @structures = @{ $positions{$pos} };
            my $row_count  = @structures;

            push @ire, ($pos - 16, $pos + 22);

            # sort by quality first, then energy
            @structures = sort {
                   $a->{quality} <=> $b->{quality}
                || $a->{energy}  <=> $b->{energy}
            } @structures;

            for my $i (0 .. $#structures) {
                my $struct = $structures[$i];
                my $quality_text = $struct->{quality} == 1 ? "Good" : "Bad";

                $structure_counter++;
                my $div_id    = "rna_ss_$structure_counter";
                my $sequence  = $struct->{sequence};
                my $structure = $struct->{structure};

                print "<tr>";

                if ($i == 0) {
                    print "<td rowspan='$row_count'>$pos</td>";
                    print "<td rowspan='$row_count'>$struct->{sequence}</td>";
                    print "<td rowspan='$row_count'>$struct->{label}</td>";
                }

                print "<td>$struct->{structure}</td>";
                printf "<td>%.2f</td>", $struct->{energy};
                print "<td>$quality_text</td>";
                
                print "<td class='no-print'>";
                print "<button onclick=\"toggleStructure$structure_counter()\" style='padding: 5px 10px;'>View</button>";
                print "<div id='$div_id' style='width: 250px; height: 250px; display: none; margin: 10px; overflow: hidden;'></div>\n";
                print "</td>";

                print "</tr>\n";

                print "<script>\n";
                print "function toggleStructure$structure_counter() {\n";
                print "  var el = document.getElementById('$div_id');\n";
                print "  if (el.style.display === 'none') {\n";
                print "    el.style.display = 'block';\n";
                print "    var container = new fornac.FornaContainer('#$div_id', {\n";
                print "      animation: false,\n";
                print "      applyForce: false,\n";
                print "      labelInterval: 0,\n";
                print "      allowPanningAndZooming: true,\n";
                print "      structurePadding: 0,\n";
                print "      drawBackground: false,\n";
                print "      layout: 'naview'\n";
                print "    });\n";
                print "    var options = {\n";
                print "      structure: '$structure',\n";
                print "      sequence: '$sequence'\n";
                print "    };\n";
                print "    container.addRNA(options.structure, options);\n";
                print "  } else {\n";
                print "    el.innerHTML = '';\n";
                print "    el.style.display = 'none';\n";
                print "  }\n";
                print "}\n";
                print "</script>\n";
            }
        }

        print "</table>\n";
        $irelineprintout = 1;

    } else {
        print "<div class='info-warning'>";
        print "<b>No Iron Response Element(s) recognized.</b>";
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
    print "<div class='info-info'>Those elements are an indication for a processing protein binding motif.</div>" if ($putativeCVfound==1);
    
    if ($putativeCVfound == 0) {
        print "<div class='info-warning'>";
        print "No CstF Motif(s) found.";
        print "</div>";
    }


    print "</div>";
}

sub stemggpairs {
    print "<div class='box-header' onclick='toggleBoxContent(this)'>Stem GG pairs:</div>";
    print "<div class=box-content>";
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

    print "<div class='info-warning'>No Stem GG Pair(s) Found.</div>" if ($stemggpairfound==0); 
    @sequ=();
    $str='';

    print "</div>";
}

sub rbp {
    print "<div class='box-header' onclick='toggleBoxContent(this)'>RNA Binding Protein Motif(s) Scan:</div>";
    print "<div class='box-content'>";

    my $yamtk_outfile = "$TEMPDIR/rbp.txt";
    my $input_seq = "$TEMPDIR/$job.seq";  

    # Run YAMTK
    my $yamtk_cmd = "$YAMTK/yamtk scan -t 0.0001 -m $RBPDB -s $input_seq -o $yamtk_outfile";
    system($yamtk_cmd);

    my @yamtk_results;
    @rbp_locs = ();  #array to store start/end positions

    # Parse yamtk output
    if (-e $yamtk_outfile) {
        open my $fh, '<', $yamtk_outfile or die "Cannot open YAMTK output: $!";
        while (my $line = <$fh>) {
            next if $line =~ /^#/;
            # next if $. == 1;
            chomp $line;
            my ($seq_name, $start, $stop, $strand, $motif, $pvalue, $score, $score_pct, $matched_seq) = split /\t/, $line;

            my ($protein_name, $motif_name);

            # extract protein and motif name
            if ($motif && index($motif, '_') != -1) {
                ($protein_name, $motif_name) = split(/_/, $motif, 2);
            } else {
                $protein_name = $motif_id;
                $motif_name   = '';
            }

            my $protein_link = "<a href=\"http://rbpdb.ccbr.utoronto.ca/proteins.php?search_term=$protein_name\" target=\"_blank\">$protein_name</a>";

            push @yamtk_results, {
                motif        => $motif,
                protein      => $protein_link,
                motif        => $motif_name,
                seq_name     => $seq_name,
                start        => $start,
                end          => $stop,
                strand       => $strand,
                score        => $score,
                pvalue       => $pvalue,
                matched_seq  => $matched_seq,
            };
        }
        close $fh;
    }

    # remove reduntant hits for each protein+domain combination
    my %best_hits;
    
    foreach my $hit (@yamtk_results) {
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
        print "<tr><th>Protein</th><th>Motif</th><th>Start</th><th>End</th><th>Strand</th><th>Score</th><th>p-value</th><th>Matched Sequence</th></tr>\n";

        foreach my $hit (@filtered_results) {
            print "<tr>";
            print "<td>$hit->{protein}</td>";
            print "<td>$hit->{motif}</td>";
            print "<td>$hit->{start}</td>";
            print "<td>$hit->{end}</td>";
            print "<td>$hit->{strand}</td>";
            print "<td>$hit->{score}</td>";
            print "<td>$hit->{pvalue}</td>";
            print "<td><tt>$hit->{matched_seq}</tt></td>";
            print "</tr>";
        }

        print "</table>";
        # print "<br>";
        # print "<div class='download-bttn'>";
        # print "<a class='no-print' href='$download_dir/rbp.txt'><button>Check Raw Output</button></a>";
        # print "</div>";

        # testing output link for redacted file
        print "<br>";
        print "<div class='sanitized-out-wrap'>";
        print sanitized_output_link(
            raw_file => $yamtk_outfile,
            out_dir  => $raw_out,
            out_path => "/tmp/jobs/job_$job_id/raw_out",
            out_name => "rbp_yamtk_raw.out",
        );
        print "</div>";
    } else {
		print "<div class='info-warning'>";
        print "<b>No Protein Binding Motif(s) found above threshold.</b>\n";
		print "</div>";
	}

    print "</div>";
}

# new ARE routine ti include all new classes of AREs (1,2,3 and core)
# based on https://pubmed.ncbi.nlm.nih.gov/22242014/; https://pmc.ncbi.nlm.nih.gov/articles/PMC2682044; https://www.pnas.org/doi/10.1073/pnas.1808696116;
# https://pmc.ncbi.nlm.nih.gov/articles/PMC8728028
sub ARE {
    $arepresent = 0;
    $are_pos = 1;
    my $are_table_opened = 0;
    my @all_ares = ();
    
    print "<div class='box-header' onclick='toggleBoxContent(this)'>AU-rich regions:</div>";
    print "<div class='box-content'>";
    
    # check papers for patterns (cmyc fore xample)
    my %are_patterns = (
        'Class I'   => qr/([au]{0,3}u{2}au{3,5}a[au]{2}[au]{0,3})/i,              # Discontinuous UUAUUUAWW
        'Class II'  => qr/((?:u{2}au{3}a[au]{0,2}){2,})/i,                        # Overlapping AUUUA motifs
        'Class III' => qr/(u{15,})/i,                                              # U-rich regions without AUUUA
        'Core'      => qr/([au]{0,3}uau{3}au[au]{0,3})/i,                         # 13-bp consensus: WWWU(AUUUA)UWWW
        'Classic'   => qr/(auuua(auuua)+)/i                                        # Original pentamer repeats (strict AUUUA)
    );
    
    # scan ARE 
    foreach my $are_type (keys %are_patterns) {
        my $pattern = $are_patterns{$are_type};
        
        while ($SEQUENCECHECKED =~ /$pattern/g) {
            my $match = $1;
            $are_len = length($match);
            $are_pos = pos($SEQUENCECHECKED) - $are_len;
            
            # Calculate AU content # validation
            my $au_count = ($match =~ tr/aAuU//);
            my $total_len = length($match);
            my $au_percent = ($au_count / $total_len) * 100;
            
            # Filter based on ARE type
            my $passes_filter = 0;
            if ($are_type eq 'Class III') {
                # Class III: U-rich only, no AUUUA required, stricter U content
                my $u_count = ($match =~ tr/uU//);
                my $u_percent = ($u_count / $total_len) * 100;
                $passes_filter = 1 if ($are_len >= 15 && $u_percent >= 70 && $match !~ /a[ut]{3}a/i);
            } else {
                # Other classes: minimum length of 9 nt and at least 60% AU content
                $passes_filter = 1 if ($are_len >= 9 && $au_percent >= 60);
            }
            
            if ($passes_filter) {
                # overlapping regions to avoid duplicates
                my $is_duplicate = 0;
                foreach my $existing (@all_ares) {
                    my ($ex_start, $ex_end) = @{$existing}[1,2];
                    if (!(($are_pos + $are_len) < $ex_start || $are_pos > $ex_end)) {
                        # overlap found< - always keep the longer one
                        if ($are_len > ($ex_end - $ex_start)) {
                            @{$existing} = ($are_type, $are_pos, $are_pos + $are_len - 1, $match, $au_percent);
                        }
                        $is_duplicate = 1;
                        last;
                    }
                }
                
                unless ($is_duplicate) {
                    push @all_ares, [$are_type, $are_pos, $are_pos + $are_len - 1, $match, $au_percent];
                }
            }
        }
    }
    
    # sort AREs by position
    @all_ares = sort { $a->[1] <=> $b->[1] } @all_ares;
    
    # Display results
    if (@all_ares > 0) {
        print "<table class='table-result'>\n";
        print "<tr><th><b>ARE</b></th><th>Start</th><th>End</th><th>Sequence</th><th>Mismatches</th><th>Score</th><th>Quality</th></tr>\n";
        $are_table_opened = 1;
        $arepresent = 1;
        
        foreach my $are (@all_ares) {
            my ($are_type, $are_start, $are_end, $are_sequence, $are_au_percent) = @$are;
            $are_len = $are_end - $are_start + 1;
            
            # Count mismatches (G and C in the sequence)
            $mismatchinare = 0;
            @aretemp = split('', $are_sequence);
            for ($arecount = 0; $arecount < @aretemp; $arecount++) {
                $mismatchinare++ if ($aretemp[$arecount] eq 'g' || $aretemp[$arecount] eq 'c');
            }
            
            # Calculate ARE score (simplified version based on AREScore algorithm)
            my $are_score = 0;
            # Count AUUUA pentamers (core motif)
            my $pentamer_count = () = $are_sequence =~ /auuua/gi;
            $are_score += $pentamer_count * 3;
            # Add points for AU content
            $are_score += ($are_au_percent / 100) * 2;
            # Bonus for U-rich flanking
            $are_score += 1 if $are_sequence =~ /^u+/i;
            $are_score += 1 if $are_sequence =~ /u+$/i;
            # Penalty for G/C content
            my $gc_count = ($are_sequence =~ tr/gc//);
            $are_score -= $gc_count * 0.5;
            $are_score = 0 if $are_score < 0;
            
            # quality based on score
            my $quality;
            if ($are_score >= 13) {
                $quality = "Excellent";
            } elsif ($are_score >= 7) {
                $quality = "Good";
            } elsif ($are_score >= 4) {
                $quality = "Moderate";
            } else {
                $quality = "Low";
            }

            printf "<tr><td>Hit</td><td>%d</td><td>%d</td><td>%s</td><td>%d</td><td>%.2f</td><td>%s</td></tr>\n",
                   $are_start, $are_end, $are_sequence, $mismatchinare, $are_score, $quality;
            
            # Store for downstream analysis (using original variable names)
            push @aurichregion, ($are_start + 1, $are_end + 1);
        }
        
        print "</table>\n";
    }
    
    if ($arepresent == 0) {
        print "<div class='info-warning'>";
        print "No AU-rich region found.";
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
        
        print "<br>";
        print "<div class='sanitized-out-wrap'>";
        print sanitized_output_link(
            raw_file => $trnascan_output,
            out_dir  => $raw_out,
            out_path => "/tmp/jobs/job_$job_id/raw_out",
            out_name => "tRNAscan_raw.out",
        );
        print "</div>";

    } else {
        print "<div class='info-warning'>";
        print "No region(s) matching a tRNA was found.";
        print "</div>";
    }
    
    print "</div>";
    # print "DEBUG: @trna_loc" if @trna_loc;
}

sub smsite {
    print "<div class='box-header' onclick='toggleBoxContent(this)'>Sm site / snRNP motif scan:</div>";
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
    
    print "<div class='info-warning'><b>No snRNP-motif(s) found.</b></div>" if ($leadlineprinted==0);
    
    ##### OUTPUT if Seq is RNA has no cds but smsite --> structured perhaps catalytic RNA !

    # if ($grepanswer=~/NO EXONS/ && $ORIGINchecked==1 && @smsite>0){
    #     print "<br>Could not detect a coding sequence on this RNA, but there are 1 or more sn-RNP motifs (sm-sites),<br>it might be possible that this is a catalytic RNA.<br>";
    # }
    print "</div>";
}

sub subopt{
    my $infile = "$TEMPDIR/$job.seq";
    my $subopt_file = "$TEMPDIR/$job.raw.subopt";
    my $sorted_subopt_file = "$TEMPDIR/$job.sorted10.subopt";

    my $cmd = "cat $infile | $VIENNARNAFOLDDIR/RNAsubopt --noLP > $subopt_file 2>/dev/null";

    system($cmd) == 0 or die "RNAsubopt failed, exit code: $?";

    my $top_n = 10;

    open(my $fh, '<', $subopt_file) or die "Cannot open $subopt_file: $!\n";

    my @rows;

    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ /^\s*$/;

        # store header line
        if ($line =~ /^>/) {
            $header = $line;
            next;
        }

        # store sequence line, e.g.
        # GAGC... -115.30 1.00
        if ($line =~ /^[ACGUTNacgutn]+\s+-?\d+(?:\.\d+)?\s+\d+(?:\.\d+)?\s*$/) {
            $seqline = $line;
            next;
        }

        # parse structure + energy
        if ($line =~ /^([().]+)\s+(-?\d+(?:\.\d+)?)\s*$/) {
            push @rows, [$1, $2 + 0];
        }
    }

    close $fh;

    @rows = sort { $a->[1] <=> $b->[1] } @rows;

    my $last = $#rows < $top_n - 1 ? $#rows : $top_n - 1;
    my @top_results = $last >= 0 ? @rows[0 .. $last] : ();

    open(my $out, '>', $sorted_subopt_file) or die "Cannot write $sorted_subopt_file: $!\n";

    print $out "$header\n"  if $header ne '';
    print $out "$seqline\n" if $seqline ne '';

    for my $r (@top_results) {
        print $out $r->[0], "\t", $r->[1], "\n";
    }

    close $out;
}

sub createfolding {
    # comment out by liang
    # new ViennaRNA does not have the old FOLD program but rather incorporated in the same program RNAfold AA
    # this is the new routine to run RNAfold and RNAplot to get the structure and the svg for display
    #### revision AA: changed rnaplot with varna to keep cosnistensy with annotated structure #########

    my $infile = "$TEMPDIR/$job.seq";
    my $foldout_file = "$TEMPDIR/$job.foldout";
    my $err_file     = "$TEMPDIR/$job.ipknot.err";

    my $seq = "";

    our $structure = "";
    our $energy    = "";
    our @structure = ();

    if ($do_pseudoknot) {
        my $cmd = "$IPKNOT/ipknot $infile > $foldout_file";
        system($cmd) == 0 or die "IPknot failed";

        open(my $fh, '<', $foldout_file) or die "Can't open output: $!";
        my $header      = <$fh>;
        $seq            = <$fh>;
        my $struct_line = <$fh>;
        close($fh) or die "Can't close output: $!";

        defined $header      or die "Missing header line in $foldout_file\n";
        defined $seq         or die "Missing sequence line in $foldout_file\n";
        defined $struct_line or die "Missing structure line in $foldout_file\n";

        chomp($header, $seq, $struct_line);

        $energy = "No energy info available when using pseudoknot prediction";
        $structure  = $struct_line;
        @structure  = split('', $structure);
    }
    else {
        my $rnafold_opts = "";
        $rnafold_opts .= "-p" if $do_centroid;

        die "Infile not found: $infile" unless -e $infile;
        die "RNAfold binary not found: $VIENNARNAFOLDDIR/RNAfold"
            unless -x "$VIENNARNAFOLDDIR/RNAfold";

        my $cmd = "cat $infile | $VIENNARNAFOLDDIR/RNAfold $rnafold_opts -d2 --noLP > $foldout_file 2>&1";
        system($cmd) == 0 or die "RNAfold failed";

        open(my $fh, '<', $foldout_file) or die "Can't open output: $!";
        my @lines = <$fh>;
        close($fh) or die "Can't close output: $!";

        chomp @lines;

        my $centroid_struct = "";
        my $centroid_energy = undef;
        my $mfe_found       = 0;

        for my $line (@lines) {
            next if $line =~ /^>/;

            if (!$seq && $line =~ /^[ACGUTNacgutn]+$/) {
                $seq = $line;
                next;
            }

            if (!$mfe_found && $line =~ /^([().]+)\s+\(\s*([^)]+?)\s*\)$/) {
                $structure = $1;
                $energy    = $2;
                @structure = split('', $structure);
                $mfe_found = 1;
                next;
            }

            if ($do_centroid && !$centroid_struct && $line =~ /^([().]+)\s+\{\s*([^}\s]+)\s+d=([^\}\s]+)\s*\}$/) {
                $centroid_struct = $1;
                $centroid_energy = $2;
                next;
            }
        }

        die "Sequence not found in RNAfold output" unless $seq;
        die "MFE structure not found in RNAfold output" unless $mfe_found;

        if ($do_centroid && $centroid_struct) {
            $structure = $centroid_struct;
            $energy    = $centroid_energy if defined $centroid_energy;
            @structure = split('', $structure);
        }
    }


    my $svg_file = "$TEMPDIR/${SEQNAMECHECKED}_ss.svg";

    my $raw_svg = render_raw_fold_svg(
        java      => $JAVA,
        jar       => $VARNA_JAR,
        sequence  => $seq,
        structure => $structure,
        svg_out   => $svg_file,
        layout    => $LAYOUT,
        padding   => 60,
        width     => 2200,
        height    => 1400,
    );

    # Read SVG file content
    open(my $svgfh, '<', $svg_file) or die "Cannot open SVG file: $!";
    my $svg_content = do { local $/; <$svgfh> };
    close($svgfh);

    # Add ID to <svg> tag if not present
    $svg_content =~ s/<svg\b/<svg id="static_rna_ss" width="100%" height="600" style="background:#FAFBFC;"/;

    # Output
    print "<h3>RNA Structure (Annotated Structure Below):</h3>\n";
    print "$svg_content\n";
    print "<p style='font-size: 0.9em; color: gray; text-align: center;'> Drag to pan, scroll to zoom</p>\n";

    # Add svg-pan-zoom script
    print "<script src='https://cdn.jsdelivr.net/npm/svg-pan-zoom\@3.6.2/dist/svg-pan-zoom.min.js'></script>\n";
    print "<script>\n";
    print "  svgPanZoom('#static_rna_ss', {\n";
    print "    zoomEnabled: true,\n";
    print "    controlIconsEnabled: true,\n";
    print "    fit: true,\n";
    print "    center: true\n";
    print "  });\n";
    print "</script>";

    
    print "<div class='sanitized-out-wrap'>";
    print sanitized_output_link(
        raw_file => $foldout_file,
        out_dir  => $raw_out,
        out_path => "/tmp/jobs/job_$job_id/raw_out",
        out_name => "RNAfold_raw.out",
    );
    print "</div>";

		
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
    if ($do_pseudoknot) {
        print "<tr><th>Energy</th><td>$energy</td></tr>\n";
    } else {
        print "<tr><th>Energy</th><td>$energy kcal/mol</td></tr>\n";
        print "<tr><th>Stems</th><td>$stem_count stem structure(s)</td></tr>\n";
        print "<tr><th>Hairpins</th><td>$hairpin_count hairpin(s)</td></tr>\n";
    }

    # the line to download suboptimal output if selected by user
    if ($do_subopt) {
        print "<tr><th>Alternate Structure(s)</th><td><tt><a class='sanitized-out-btn' href='/tmp/jobs/job_$job_id/$job.sorted10.subopt' target='_blank' rel='noopener noreferrer'>View Predicted Alternate Structures</a></tt></td></tr>\n";
    }
    
    
    print "</table>";
    # Prediction logic
    my $structure_length = scalar @structure;
    my $avg_spacing = $structure_length / ($stem_count || 1);
    
    if ($stem_count >= 15 && $avg_spacing < 100) {
        $comment = "Highly structured RNA.\n";
    } elsif ($stem_count >= 1) {
        $comment = "Some secondary structure detected, may have biological significance.\n";
    } else {
        $comment = "Minimal secondary structure.\n";
    }
    print "<div class='info-info'>$comment</div>";
    
    if ($hairpin_count > 0 && $stem_count > 0 && $hairpin_count / $stem_count > 0.7) {
        print "\t\tHigh hairpin content characteristic of miRNA precursors\n";
    }
    
    
    print "<br>";
    
    return ($stem_count, $hairpin_count, $energy);
}

sub composition {
    
    print "<div class='box'>";
    print "<div class='box-header' onclick='toggleBoxContent(this)'>Composition:</div>";
    print "<div class=box-content>";
    # Count bases
    my $a_count = ($SEQUENCECHECKED =~ tr/a/a/);
    my $u_count = ($SEQUENCECHECKED =~ tr/u/u/);
    my $g_count = ($SEQUENCECHECKED =~ tr/g/g/);
    my $c_count = ($SEQUENCECHECKED =~ tr/c/c/);
    
    # Calculate percentages
    my $a_percent = $SEQUENCELENGTH > 0 ? sprintf("%.2f", ($a_count / $SEQUENCELENGTH) * 100) : 0;
    my $u_percent = $SEQUENCELENGTH > 0 ? sprintf("%.2f", ($u_count / $SEQUENCELENGTH) * 100) : 0;
    my $g_percent = $SEQUENCELENGTH > 0 ? sprintf("%.2f", ($g_count / $SEQUENCELENGTH) * 100) : 0;
    my $c_percent = $SEQUENCELENGTH > 0 ? sprintf("%.2f", ($c_count / $SEQUENCELENGTH) * 100) : 0;
    
    # Calculate GC and AU content
    my $gc_content = $SEQUENCELENGTH > 0 ? sprintf("%.2f", (($g_count + $c_count) / $SEQUENCELENGTH) * 100) : 0;
    my $au_content = $SEQUENCELENGTH > 0 ? sprintf("%.2f", (($a_count + $u_count) / $SEQUENCELENGTH) * 100) : 0;

    
    
    print "<table class='table-result'>\n";
    print "<tr><th>Base</th><th>Adenine (A)</th><th>Uracil (U)</th><th>Guanine (G)</th><th>Cytosine (C)</th></tr>";
    print "<tr><th>Count</th><td>$a_count</td><td>$u_count</td><td>$g_count</td><td>$c_count</td></tr>";
    print "<tr><th>Percentage (%)</th><td>$a_percent%</td><td>$u_percent%</td><td>$g_percent%</td><td>$c_percent%</td></tr>";
    print "</table>";

    print "<br>";

    print "<table class='table-result'>\n";
    print "<tr><th>Content Type</th><th>GC Content</th><th>AU Content</th></tr>";
    print "<tr><th>Percentage (%)</th><td style='color:green;'>$gc_content%</td><td style='color:red;'>$au_content%</td></tr>";
    print "</table>";

    print "</div>";
    print "</div>";
}

sub checkstemsonly {
    # revamped checkstemsonly based on the logic of the new checkstems
    my $inseq = $_[0];
    my $inwhattodo = $_[1];   # check if structure given, 0 = strcuture, 1 = sequence which will be folded
    my $struct;
    my @structure;
    my $energy;
    my $structure_string = '';
    
    # Handle RNA folding or structure input
    if ($inwhattodo == 1 && length $inseq <= $MAXFOLDINGLENUTR) {
        my $rnafold_opts = $do_centroid ? "-p" : "";
        my @struct = `echo $inseq | $VIENNARNAFOLDDIR/RNAfold $rnafold_opts -d2 --noLP 2>&1`;
        if ($? != 0) {
            return ("RNAfold failed", 0);
        }

        chomp @struct;

        my $centroid_struct = "";
        my $centroid_energy = undef;
        my $mfe_found       = 0;

        for my $line (@struct) {
            next if $line =~ /^>/;

            if (!$mfe_found && $line =~ /^([().]+)\s+\(\s*([^)]+?)\s*\)$/) {
                @structure = split('', $1);
                $structure_string = $1;
                $energy    = $2;
                $mfe_found = 1;
                next;
            }

            if ($do_centroid && !$centroid_struct && $line =~ /^([().]+)\s+\{\s*([^}\s]+)\s+d=([^\}\s]+)\s*\}$/) {
                $centroid_struct = $1;
                $centroid_energy = $2;
                next;
            }
        }

        if ($do_centroid && $centroid_struct) {
            @structure = split('', $centroid_struct);
            $structure_string = $centroid_struct;
            $energy    = $centroid_energy if defined $centroid_energy;
        }

        return ("RNAfold parse failed", 0) unless @structure;
    }
    if ($inwhattodo == 1 && length $inseq > $MAXFOLDINGLENUTR) {
        return ("Too long for detection", 1); 
    }
    if ($inwhattodo == 0) {
        @structure = split ('', $inseq);
        $structure_string = $inseq;
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
    # added option to return structure as well
    # print "structure string: $structure_string\n";
    return ($fldstemsauf, $fldstemszu, $energy, $structure_string);
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
    print "<script src='https://d3js.org/d3.v3.min.js'></script>";
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
        my $start = min($from, $to);
        my $length = abs($to - $from) + 1;
        my $motif_seq = substr($SEQUENCECHECKED, $start - 1, $length);
	
        if ($from > $to) {
            
            $motif_seq = reverse($motif_seq);
            $motif_seq =~ tr/acgu/ugca/;
        }
        # print ("Motif Seq: $motif_seq\n");

        ### Fold subsequence using RNAfold
        # my $cmd = qq{echo "$motif_seq" | $VIENNARNAFOLDDIR/RNAfold --noPS 2>/dev/null};
        # my $foldout = `$cmd`;  
        # my @fold_lines = split("\n", $foldout);

        # my $dot_bracket = '';
        # $dot_bracket = $1 if $fold_lines[1] && $fold_lines[1] =~ /([().]+)\s+\([^)]+\)/;
        
        # replaced running rnafold again to use chemstems only
        my @check_stem = &checkstemsonly($motif_seq, 1);

        my $dot_bracket = $check_stem[3];

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
        print "<tr><th>Match</th><th>Family</th><th>From</th><th>To</th><th>Score</th><th>E-value</th><th>Description</th><th class='no-print'>View Structure</th></tr>\n";

        foreach my $hit (@results) {
            print "<tr>";
            print "<td>$hit->{match}</td>";
            print "<td>$hit->{family}</td>";
            print "<td>$hit->{from}</td>";
            print "<td>$hit->{to}</td>";
            print "<td>$hit->{score}</td>";
            print "<td>$hit->{evalue}</td>";
            print "<td>$hit->{description}</td>";
            print "<td class='no-print'>$hit->{structure}</td>";
            print "</tr>";
        }
        
        print "</table>";

        print "<br>";
        print "<div class='sanitized-out-wrap'>";
        print sanitized_output_link(
            raw_file => $output_file,
            out_dir  => $raw_out,
            out_path => "/tmp/jobs/job_$job_id/raw_out",
            out_name => "cmscan_rfam_raw.out",
        );
        print "</div>";        

	} else {
        print "<div class='info-warning'>";
		print "No RNA motif(s) recognized\n";
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

    my ($orf_start, $peptide_length, $label, $strand);

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
        print "</table>";

        print "<br>";
        print "<div class='sanitized-out-wrap'>";
        print sanitized_output_link(
            raw_file => "$cpc_output.txt",
            out_dir  => $raw_out,
            out_path => "/tmp/jobs/job_$job_id/raw_out",
            out_name => "cpc2_raw.out",
        );
        print "</div>";

    } else {
        print "<div class='info-warning'>";
        print "No result in CPC2 output\n";
        print "</div>";
    }

    # Inline Kozak (checks only if we have a valid ORF start)
    my $ok_kozak = 1;  # assume OK unless we detect a strong Kozak
    if (defined $orf_start && $orf_start > 6 && ($orf_start + 3) <= length($SEQUENCECHECKED)) {
        my $win = substr($SEQUENCECHECKED, $orf_start - 6, 10);  # extract 10 nt widnow for koyak check -6..+3 around AUG
        $win = uc $win;

        # flip to coding orientation if negative strand
        if (defined $strand && $strand eq '-') {
            $win =~ tr/ACGU/TGCA/; $win = reverse $win;
        }

        # indices in coding orientation: -3 = [3], +4 = [9]
        my $minus3 = substr($win, 3, 1);
        my $plus4  = substr($win, 9, 1);

        my $strong_kozak = (($minus3 eq 'A' || $minus3 eq 'G') && $plus4 eq 'G') ? 1 : 0;
        $ok_kozak = !$strong_kozak;  # block lncRNA if strong Kozak is present
    }

    # check for lncRNA
    if (defined $label && $label eq 'noncoding' && $SEQUENCELENGTH > 200 && $peptide_length < 100 && $ok_kozak == 1) {
        print "<div class='info-info' style='overflow: hidden; text-overflow: ellipsis;'>";
        print "<p style='white-space: normal; word-break: break-word;'>Based on its non-coding potential, sequence length, small open reading frame (ORF) length, and the absence or weakness of a Kozak sequence, the given sequence may be classified as a long non-coding RNA (lncRNA). Please check the folding below.</p>";
        # upper page shows the folding # to be added
        print "</div>";
    }

    # Extract cds sequnece 
    my $cds_nt = substr($SEQUENCECHECKED, $orf_start - 1, $peptide_length * 3);

    if (defined $label && $label eq 'coding' && $orf_start > 0 && $peptide_length > 0) {

    my $orf_end = $orf_start + ($peptide_length * 3) - 1;
    

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

            print "<br>";
            print "<div class='sanitized-out-wrap'>";
            print sanitized_output_link(
                raw_file => $mirbase_out,
                out_dir  => $raw_out,
                out_path => "/tmp/jobs/job_$job_id/raw_out",
                out_name => "nhmmer_mirscan_raw.out",
            );
            print "</div>";

            
		} else {
            print "<div class='info-warning'>";
			print "No region(s) matching a mircroRNA found.\n";
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

        print "<br>";
        print "<div class='sanitized-out-wrap'>";
        print sanitized_output_link(
            raw_file => $ribosw_out,
            out_dir  => $raw_out,
            out_path => "/tmp/jobs/job_$job_id/raw_out",
            out_name => "cmscan_riboswitch_raw.out",
        );
        print "</div>"; 

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
    my $mirna_db             = "$MIRBASE/mature.fa";
    my $miranda_raw_out      = "$TEMPDIR/$job.miranda.out";
    my $miranda_out          = "$TEMPDIR/$job.miranda.tsv";
    my $utr_fasta            = "$TEMPDIR/$job.utr3.fa";

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

    my $cmd = "python3 $MIRANDA/miranda_wrapper.py --parsed_out $miranda_out --miranda_bin $MIRANDA/miranda --tmpdir $TEMPDIR --job_id $job  $mirna_db $utr_fasta $miranda_raw_out";
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
    
    print "<br>";
    print "<div class='sanitized-out-wrap'>";
    print sanitized_output_link(
        raw_file => $miranda_raw_out,
        out_dir  => $raw_out,
        out_path => "/tmp/jobs/job_$job_id/raw_out",
        out_name => "miranda_mirtarget_raw.out",
    );
    print "</div>";

    print "</div>";
    return @regions;
}

########################
# augustus replacing the old genscan; need to check the whole sub for errors
# Flag inferred UTRs as below confidence when:
# ORF is very short
# ORF start is <15 nt from sequence start (no room for 5' UTR)
# Sequence ends shortly after ORF (truncated transcript)
# Provide scoring or confidence per UTR:
# E.g., "Predicted 5' UTR: 37 bp (contains weak Shine-Dalgarno motif, deltaG = -3.2 kcal/mol)"
# Refactored AUGUSTUS + UTR + PolyA logic to handle multiple transcripts

sub AUGUSTUS {

    print "<div class='box-header' onclick='toggleBoxContent(this)'>Gene Prediction Analysis:</div>";
    print "<div class='box-content'>";
    my ($species) = @_;
    $species ||= "human";
    my $utr_flag = ($species =~ /^(human|fly|zebrafish)$/i) ? "--UTR=on" : "";

    my $output_gff = "$TEMPDIR/$job.augustus";
    my $input_dna  = "$TEMPDIR/$job.dna";

    my $augustus_cmd = "$AUGUSTUS --softmasking=0 --protein=on --/augustus/verbosity=0 $utr_flag --species=$species $input_dna > $output_gff 2>&1";
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

        # lncRNA check for noncoding transcript
        if ($SEQUENCELENGTH > 200) {
            print "<div class='info-info' style='word-wrap: break-word; overflow-wrap: break-word; white-space: normal;'>";
            print "<p style='word-wrap: break-word;'>Based on its non-coding potential, and overall sequence length, the given sequence may be classified as a long non-coding RNA (lncRNA). Please check the folding of the sequence below.</p>";
            print "</div>";
        }
        
    } else {
        # Process each transcript
        foreach my $tid (sort keys %transcripts) {
            my $model = $transcripts{$tid};
            
            my ($gene_id) = $tid =~ /^(g\d+)\./;
            print "<h3>Transcript $i (ID: $tid";
            print ", gene: $gene_id" if $gene_id;
            print ")</h3><br>";
            
            # Kozak check for this transcript
            my $orf_start = $model->{cds_start};
            my $strand = $model->{strand};
            my $ok_kozak = 1;  # assume OK unless we detect a strong Kozak
            
            if (defined $orf_start && $orf_start > 6 && ($orf_start + 3) <= length($SEQUENCECHECKED)) {
                my $win = substr($SEQUENCECHECKED, $orf_start - 6, 10);  # extract 10 nt widnows to check kozak -6..+3 around AUG (start codon)
                $win = uc $win; 
                
                # check negative strand
                if (defined $strand && $strand eq '-') {
                    $win =~ tr/ACGU/UGCA/; 
                    $win = reverse $win; 
                }
                
                # Indices in coding orientation: -3 = [3], +4 = [9]
                my $minus3 = substr($win, 3, 1);
                my $plus4  = substr($win, 9, 1);
                
                my $strong_kozak = (($minus3 eq 'A' || $minus3 eq 'G') && $plus4 eq 'G') ? 1 : 0;
                $ok_kozak = !$strong_kozak;  # flag is 1 if weak Kozak
                
                
                # print "Kozak check = $ok_kozak / $strong_kozak" # debug
            }
            
            $i++;
            
            # Check if transcript has CDS (coding) or not (noncoding)
            my $has_cds = exists $model->{cds} && @{ $model->{cds} } > 0;
            
            # Calculate peptide length
            my $peptide_length = 0;
            if (exists $model->{protein}) {
                my $prot = $model->{protein};
                $prot =~ s/\s//g;
                $peptide_length = length($prot);
            }

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
            
            print "<br>";
            print "<div class='sanitized-out-wrap'>";
            print sanitized_output_link(
                raw_file => $output_gff,
                out_dir  => $raw_out,
                out_path => "/tmp/jobs/job_$job_id/raw_out",
                out_name => "augustus_raw.out",
            );
            print "</div>";   
            
            # Check if noncoding (no CDS predicted) AND no protein sequence
            if (!$has_cds && !exists $model->{protein} && $SEQUENCELENGTH > 200) {
                print "<div class='info-info' style='word-wrap: break-word; overflow-wrap: break-word; white-space: normal;'>";
                print "<p style='word-wrap: break-word;'><i>Based on its non-coding potential, and sequence length, the given sequence may be classified as a long non-coding RNA (lncRNA). Please check the folding of the sequence below.</i></p>";
                print "</div>";
                
            } elsif ($has_cds && $peptide_length < 100 && $ok_kozak == 1 && $SEQUENCELENGTH > 200) {
                # lncRNA check for coding transcript with short ORF and weak Kozak
                print "<div class='info-info' style='word-wrap: break-word; overflow-wrap: break-word; white-space: normal;'>";
                print "<p style='word-wrap: break-word;'>Based on its sequence length, small open reading frame (ORF) length, and the absence or weakness of the Kozak sequence, the given sequence may be classified as a long non-coding RNA (lncRNA) Please check the folding of the sequence below.</p>";
                print "</div>";
            }
            
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

# planned upgrade to UTR subroutine
# new subroutine for UTR prediction with fallback to longest ORF inference when AUGUSTUS predictions are not available.
# Also adds motif scanning and confidence flagging based on ORF length and proximity to sequence ends.


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

sub createfoldingpictureFornac {
    my $seq_file       = "$TEMPDIR/$job.seq";
    my $foldout_file   = "$TEMPDIR/$job.foldout";

    # Print HTML content
    print "<div class='box-header' onclick='toggleBoxContent(this)'>RNA Structure Analysis:</div>";
    print "<div class='box-content-structure'>";
    
    if ($SEQUENCELENGTH <= $MAXFOLDINGLEN) {
        
        my $svg_file = "$TEMPDIR/${SEQNAMECHECKED}_ss.svg";
        my $ps_url  = "/tmp/jobs/job_$job/${SEQNAMECHECKED}_ss.ps";
        my $svg_url = "/tmp/jobs/job_$job/${SEQNAMECHECKED}_ss.svg";

        open(my $fh, '<', $foldout_file) or die "Can't open output: $!";
        my @lines = <$fh>;
        close($fh);

        chomp @lines;

        my $seq             = "";
        my $centroid_struct = "";
        my $mfe_found       = 0;

        for my $line (@lines) {
            next if $line =~ /^>/;

            if (!$seq && $line =~ /^[ACGUTNacgutn]+$/) {
                $seq = $line;
                next;
            }

            # first normal structure line = MFE
            if (!$mfe_found && $line =~ /^([().]+)\s+\(([-\d.]+)\)$/) {
                $structure = $1;
                $energy    = $2;
                @structure = split('', $structure);
                $mfe_found = 1;
                next;
            }

            # centroid line from RNAfold -p
            if ($do_centroid && !$centroid_struct && $line =~ /^([().]+)\s+\{.+\}$/) {
                $centroid_struct = $1;
                next;
            }
        }

        die "Sequence not found in RNAfold output" unless $seq;
        die "MFE structure not found in RNAfold output" unless $mfe_found;

        # override displayed structure if centroid requested and found
        if ($do_centroid && $centroid_struct) {
            $structure  = $centroid_struct;
            @structure  = split('', $structure);
        }
          
        print "<div id='rna_ss'></div>";
        print "<p style='font-size: 0.9em; color: gray; text-align: center;'> Drag to pan, scroll to zoom</p>";
        print "<b>Download Folding As: </b>\n";
        print "<a class='no-print' href='$svg_url' target='_blank'><button>SVG File</button></a>";
        print "<a class='no-print' href='$ps_url' target='_blank'><button>PS File</button></a>\n";

        # Include required scripts
        print "<link rel='stylesheet' href='/css/fornac.css'>\n";
        print "<script src='https://d3js.org/d3.v3.min.js'></script>\n";
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

		#ire
		for (my $j = 0; $j < @ire; $j += 2) {
            my ($from, $to) = @ire[$j, $j + 1];

            ($from, $to) = ($to, $from) if $from > $to;

            for (my $i = $from; $i <= $to; $i++) {
                $color_text .= "$i:lightcoral ";
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
        print "</script>";

        # creating legend for visual interpretation 
        print "<div class='legend'>";
        print "<div class='legend-title'>Legend:</div>";
        print "<div class='legend-items'>";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: lightblue; margin-left: 10px;'></span> UTRs\t";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: green; margin-left: 10px;'></span> CDS\t";
		print "<span style='display: inline-block; width: 12px; height: 12px; background: lime; margin-left: 10px;'></span> PolyA Signal/PolyA Tail\t";
		print "<span style='display: inline-block; width: 12px; height: 12px; background: red; margin-left: 0px;'></span> RNA Motifs\t";
		print "<span style='display: inline-block; width: 12px; height: 12px; background: yellow; margin-left: 10px;'></span> MiRNA\t";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: orange; margin-left: 10px;'></span> TRNA\t";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: purple; margin-left: 10px;'></span> SM-site/snRNP-motif\t";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: royalblue; margin-left: 10px;'></span> Riboswitches(s)";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: cyan; margin-left: 10px;'></span> TRANS-splicing\t";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: lightcoral; margin-left: 10px;'></span> IRE(s)\t";
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

sub createfoldingpicture {
    #### revision AA: updated folding structure , uses VARNA for visualization and annotation, with custom feature highlighting. 
    #### Also added download options for both raw and annotated structures which was not possible with FORNA
    #### This also helps in the structure simialrity of the both the figures adn visualization looks is convinient to interpret

    my $seq_file       = "$TEMPDIR/$job.seq";
    my $foldout_file   = "$TEMPDIR/$job.foldout";
    my $annotated_svg    = "$TEMPDIR/${SEQNAMECHECKED}_ss_ann.svg";

    # Print HTML content
    print "<div class='box-header' onclick='toggleBoxContent(this)'>RNA Structure Analysis:</div>";
    print "<div class='box-content-structure'>";


    if ($SEQUENCELENGTH <= $MAXFOLDINGLEN) {
        
        my $svg_file = "$TEMPDIR/${SEQNAMECHECKED}_ss.svg";
        my $ps_url  = "/tmp/jobs/job_$job/${SEQNAMECHECKED}_ss.ps";
        my $annotated_svg_url = "/tmp/jobs/job_$job/${SEQNAMECHECKED}_ss_ann.svg";
        my $svg_url = "/tmp/jobs/job_$job/${SEQNAMECHECKED}_ss.svg";

        open(my $fh, '<', $foldout_file) or die "Can't open output: $!";
        my @lines = <$fh>;
        close($fh);

        chomp @lines;

        my $seq             = "";
        my $centroid_struct = "";
        my $mfe_found       = 0;

        if ($do_pseudoknot) {

            my $header      = $lines[0];
            my $seq_line    = $lines[1];
            my $struct_line = $lines[2];

            $seq = $seq_line;

            # if you still use -E and it works sometimes:
            if ($header =~ /\(e=([-\d.]+)\)/) {
                $energy = $1;
            } else {
                $energy = "No energy info available when using pseudoknot prediction";
            }

            $structure = $struct_line;
            @structure = split('', $structure);
            $mfe_found = 1;
        }
    
        else {
            for my $line (@lines) {
                next if $line =~ /^>/;

                if (!$seq && $line =~ /^[ACGUTNacgutn]+$/) {
                    $seq = $line;
                    next;
                }

                # first normal structure line = MFE
                if (!$mfe_found && $line =~ /^([().]+)\s+\(\s*([^)]+?)\s*\)$/) {
                    
                    $structure = $1;
                    $energy    = $2;
                    @structure = split('', $structure);
                    $mfe_found = 1;
                    next;
                }

                # centroid line from RNAfold -p
                if ($do_centroid && !$centroid_struct && $line =~ /^([().]+)\s+\{\s*([^}\s]+)\s+d=([^\}\s]+)\s*\}$/) {
                    $centroid_struct = $1;
                    next;
                }
            }

            die "Sequence not found in RNAfold output" unless $seq;
            die "MFE structure not found in RNAfold output" unless $mfe_found;

            # override displayed structure if centroid requested and found
            if ($do_centroid && $centroid_struct) {
                $structure  = $centroid_struct;
                @structure  = split('', $structure);
            }
        }

        ## keeping debug lines for later use
        # print "DEBUG: MFE structure = $structure, energy = $energy\n";
        # print "DEBUG: Centroid structure = $centroid_struct\n" if $do_centroid;
        # print "DEBUG: Sequence length = " . length($seq) . "\n";

        my %feature_styles = (
            utr          => { color => '#ADD8E6', mode => 'region' },
            exons        => { color => '#268b26', mode => 'region' },
            polyA_signal => { color => '#00FF00', mode => 'region' },
            polyA_tail   => { color => '#32CD32', mode => 'region' },
            motif        => { color => '#FF0000', mode => 'region' },
            mirna        => { color => '#FFD700', mode => 'region' },
            trna         => { color => '#FFA500', mode => 'region' },
            sm           => { color => '#FF00FF', mode => 'region' },
            riboswitch   => { color => '#4169E1', mode => 'region' },
            transsplice  => { color => '#00FFFF', mode => 'region' },
            ire          => { color => '#F08080', mode => 'region' },
            rbp          => { color => '#FFC0CB', mode => 'region' },
        );

        my %features = (
            utr          => \@utr,
            exons        => \@exons,
            polyA_signal => \@polyasignal,
            polyA_tail   => \@polyatail,
            motif        => \@rna_motif,
            mirna        => \@mirna_loc,
            trna         => \@trna_loc,
            sm           => \@sm,
            riboswitch   => \@ribosw,
            transsplice  => \@transsplicing,
            ire          => \@ire,
            rbp          => \@rbp_locs,
        );

        my $annotated_svg = render_annotated_fold_svg(
            java      => $JAVA,
            jar       => $VARNA_JAR,
            #foldout   => $foldout_file,
            sequence  => $seq,
            structure => $structure,
            svg_out   => $annotated_svg,
            layout    => $LAYOUT,
            features  => \%features,
            styles    => \%feature_styles,
            padding   => 60,
            width     => 2200,
            height    => 1400,
        );

        #print "DEBUG: Annotated SVG generated at $annotated_svg\n";

        # same logic as above to read and display the annotated SVG (sub createfolding)

        # Read SVG file content
        open(my $anno_svgfh, '<', $annotated_svg) or die "Cannot open annotated SVG file: $!";
        my $anno_svg_content = do { local $/; <$anno_svgfh> };
        close($anno_svgfh);

        # Add ID to <svg> tag if not present
        $anno_svg_content =~ s/<svg\b/<svg id="rna_ss" width="100%" height="800" style="background:#FAFBFC;"/;

        # Output
        print "<h3>RNA Structure (Annotated Structure Below):</h3>\n";
        print "$anno_svg_content\n";
        print "<p style='font-size: 0.9em; color: gray; text-align: center;'> Drag to pan, scroll to zoom</p>\n";

        # Add svg-pan-zoom script
        print "<script src='https://cdn.jsdelivr.net/npm/svg-pan-zoom\@3.6.2/dist/svg-pan-zoom.min.js'></script>\n";
	    print "<script>\n";
        print "  svgPanZoom('#rna_ss', {\n";
        print "    zoomEnabled: true,\n";
        print "    controlIconsEnabled: true,\n";
        print "    fit: true,\n";
        print "    center: true\n";
        print "  });\n";
        print "</script>\n";


        print "<b>Download Annotated Structure: </b>\n";
        print "<a class='no-print' href='$annotated_svg_url' target='_blank'><button>Annotated RNA Structure</button></a>";
        print "<a class='no-print' href='$svg_url' target='_blank'><button>Raw RNA Structure</button></a>";


        # creating legend for visual interpretation 
        print "<div class='legend'>";
        print "<div class='legend-title'>Legend:</div>";
        print "<div class='legend-items'>";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: lightblue; margin-left: 10px;'></span> UTRs\t";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: green; margin-left: 10px;'></span> CDS\t";
		print "<span style='display: inline-block; width: 12px; height: 12px; background: lime; margin-left: 10px;'></span> PolyA Signal/PolyA Tail\t";
		print "<span style='display: inline-block; width: 12px; height: 12px; background: red; margin-left: 0px;'></span> RNA Motifs\t";
		print "<span style='display: inline-block; width: 12px; height: 12px; background: yellow; margin-left: 10px;'></span> MiRNA\t";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: orange; margin-left: 10px;'></span> TRNA\t";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: purple; margin-left: 10px;'></span> SM-site/snRNP-motif\t";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: royalblue; margin-left: 10px;'></span> Riboswitches(s)";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: cyan; margin-left: 10px;'></span> TRANS-splicing\t";
        print "<span style='display: inline-block; width: 12px; height: 12px; background: lightcoral; margin-left: 10px;'></span> IRE(s)\t";
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

    print "<tr><th>Structure</th><th>Location(s)</th></tr>\n";

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
	
	if (@ire) {
        my $ire_str = format_flat_ranges(@ire);
        print "<tr><th>IRE(s)</th><td>$ire_str</td></tr>\n";
    }
	
    if (@polyasignal) {
        my $polyasignal_str = format_flat_ranges(@polyasignal);
        print "<tr><th>PolyA motif(s)</th><td>$polyasignal_str</td></tr>\n";
    }

    if (@polyatail) {
        my $polyatail_str = format_flat_ranges(@polyatail);
        print "<tr><th>PolyA tail</th><td>$polyatail_str</td></tr>\n";
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

# Download options
print "<br>";
print "<div class='download-bttn'>";
print "<a class='no-print' href='localhost/cgi-bin/download_pdf.cgi?job=$job'><button>Download Results as PDF</button></a>";
print "<a class='no-print' href='$download_dir/result_file.txt'><button>Download Results as TEXT</button></a>";
print "</div>";

write_file("$TEMPDIR/result.txt", "done\n");


print "</div></main></body></html>";

sub make_txt {

    my $html_file = "$TEMPDIR/result.html";
    my $txt_file = "$TEMPDIR/result_file.txt";

    my $cmd = "python3.11 $BINDIR/download_txt.py -i $html_file -o $txt_file";
    my $output = `$cmd`;
    # print "Python Error: $output\n" if $?;

}

close $out;
