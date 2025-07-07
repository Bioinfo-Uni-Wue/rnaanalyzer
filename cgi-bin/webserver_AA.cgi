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
$MAXFOLDINGLEN=5000;
$MAXFOLDINGLENUTR=5000;
$MAXFORNALENGTH=5000;
$maxcoloredseqlen=10000;

# Create a new CGI object
my $cgi = CGI->new;


$job_id = $ARGV[0] // $cgi->param("job_id");  # passed from batch script
$job = $job_id;
# Store input in job folder

$TEMPDIR = abs_path("../tmp/jobs/job_$job");


unless (-d $TEMPDIR) {
    make_path($TEMPDIR) or die "Can't create job dir: $TEMPDIR";
}
print "<p>Job submitted: $job</p>";


open my $log, ">>", "../tmp/backend_log.txt";
print $log scalar localtime() . " Received job: $job\n";
close $log;

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
my $dnarna              = $params->{dnarna};
my $SEQNAMECHECKED      = $params->{sequence_name};
my $SEQUENCECHECKED     = $params->{sequence_clean};
my $SEQUENCELENGTH      = $params->{sequence_length};

my $html_file = "$TEMPDIR/result.html";

open my $out, ">", $html_file or die "Can't write result page: $!";
select $out;  # Redirect STDOUT to file

# Now all your print statements go into result.html
# print "<!DOCTYPE html>\n";
# print "<html><head><meta charset='UTF-8'><title>Job $job Results</title></head><body>";
print <<'HTML';
<html><head>
  <title>Batch Results</title>
  <link rel="stylesheet" href="/css/results.css">
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
<h2>Results</h2>
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
    our @rna_motif;

    
    print $cgi->h2("Here are the results for JOB ID: $job with sequence name: ". CGI::escapeHTML($SEQNAMECHECKED));

    # running analysis
    &analysis;
}



##calls all the new subrotines
sub analysis {
    chdir $TEMPDIR;

    print "<pre>";
    &TRANS if $do_TRANS;
		
	&IRE if $do_IRE;

    if (length $SEQUENCECHECKED <= $maxcoloredseqlen) {
        &createfolding;
        &checkstems;
    } else {
        print "<br><b>Length:</b>        $SEQUENCELENGTH";
        print "     *some information is only available up to $MAXFOLDINGLEN nt\n" if ($SEQUENCELENGTH > $MAXFOLDINGLEN);
        print "<br><b>Origin:</b>        $dnarna<br>";
    }

    &ARE;
    

    print "<br><b>Catalytic RNA:</b><br>";
    &smsite;

    &tRNA if $do_trna;
    &RNAMOTIF if $do_rnamotif;
    &microRNA if $do_mirna;

    # --- Capture transcript models ---
    my %transcripts;
    if ($do_augustus) {
        %transcripts = %{ AUGUSTUS() };   # Returns hashref
    } else {
        %transcripts = %{ CPC2() };       # Also returns hashref
        $cpc = "TRUE"; 
    }

    foreach my $tid (keys %transcripts) {
        my $model = $transcripts{$tid};

        my ($utr5_ref, $utr3_ref, $utrprintout_ref, $utr_coords_ref) = predict_utrs(
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
        }

        # Optionally store results back into the model
        $model->{utr5} = $utr5_ref;
        $model->{utr3} = $utr3_ref;
        $model->{utr_coords} = $utr_coords_ref;
    }

    normalize_transcript_features(
        \%transcripts,
        \@exons,
        \@utr,
        \@polyasignal,
        \@mirnatarget,
        \@structured_regions
    );

    miRNAtarget(\%transcripts) if $do_mirnatarget;

    # Optional: enable if scan_rbs supports multi-transcript
    # scan_rbs(\%transcripts);

    &csfce;
    &protA1bisite;

    &createfoldingpicture;

    &location_table;

    if (length $SEQUENCECHECKED <= $MAXFOLDINGLEN) {
        &stemggpairs;
    }

    &drawcoloredsequence;

    print "</pre>";
}

sub TRANS{
	#this checkes the ciona-consensus !		
			@transcionareturnvalues=RNASERVER::TRANS2::ciona($SEQUENCECHECKED); ## Ist das der CIONA ??? Auf jeden Fall aber SCHISOTSOMA
			print "<b>Trans-Splicing:</b><br>";
			#print "<b> <big>Putative trans-splicing Schistosoma-consensus</b></big> search:";
			if (@transcionareturnvalues==1) {
				#print "No hit detected<br>";
				print " Schistosoma: none<br>";
			}
			else {
				$hits=pop @transcionareturnvalues;
				for ($count=0;$count<$hits;$count++){
					print " Schistosoma: <br>";
					print "  Position:   $transcionareturnvalues[$count*6+0]<br>";
					$transcionareturnvalues[$count*6+4]=uc ($transcionareturnvalues[$count*6+4]);
					$transcionareturnvalues[$count*6+2]=uc ($transcionareturnvalues[$count*6+2]);
					print "  Stem1:      $transcionareturnvalues[$count*6+4]<br>";
					print "  Structure:  $transcionareturnvalues[$count*6+5]<br>";
					print "  Energy:     $transcionareturnvalues[$count*6+3]<br>";
					print "  Sm-Site:    $transcionareturnvalues[$count*6+2] at pos: $transcionareturnvalues[$count*6+1]<br>";
					@transsplicing=(@transsplicing,$transcionareturnvalues[$count*6+0]-35,($transcionareturnvalues[$count*6+1]+length $transcionareturnvalues[$count*6+2])-1); #for the formatted output of sequence
				}
			}
			#this checkes the c. elegans consensus !!
			@transcelegansvalues=RNASERVER::TRANS2::celegans($SEQUENCECHECKED);	
			if (@transcelegansvalues==1){
				#print "No hit detected<br>";
				print " C. elegans:  none<br>";
			}
			else {
				$hits=pop @transcelegansvalues;
				for ($count=0;$count<$hits;$count++){
					$transcelegansvalues[$count*10+5]=uc($transcelegansvalues[$count*10+5]);
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
					
					
					
					@transsplicing=(@transsplicing,$transcelegansvalues[$count*10+0]-21,($transcelegansvalues[$count*10+9]+(length $transcelegansvalues[$count*10+5])+(length $transcelegansvalues[$count*10+6])));
				}
			}
		
}

sub IRE{
	#So jetzt machen wir den gleichen Spass wie mit trans
			@irereturnvalues=RNASERVER::IRE::suboptimalfindire($SEQUENCECHECKED);
			$irelineprintout=0;
			if (@irereturnvalues>1){
				########new #####
				
				$posprintout=0;
				for ($count=0;$count<=@irereturnvalues-1;$count=$count+7){
				    print "<br><b>Iron-resp Ele.:</b><br>" if ($irelineprintout==0);
				    $irelineprintout=1;
				    
				    if ($posprintout!=$irereturnvalues[$count+0]) {
				    @ire=(@ire,$irereturnvalues[$count+0]-16,$irereturnvalues[$count+0]+22) if ($posprintout!=$irereturnvalues[$count+0]);
				    print " Position:     $irereturnvalues[$count+0]<br>" if ($posprintout!=$irereturnvalues[$count+0]);
				   		
				    print " Sequence:     $irereturnvalues[$count+2]<br>" if ($posprintout!=$irereturnvalues[$count+0]); 
			        }
				    print " Structure:    $irereturnvalues[$count+3]";
			    
				    printf ("  %2f kcal/mol ",$irereturnvalues[$count+4]);
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
        	print "<br><br><b>CstF:</b>          start      mismatch<br>" if ($putativeCVfound==0);
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
		    	print "<b>StemGGpair:</b>    start     -     end<br>" if ($stemggleadlineprinted==0);
			$stemggleadlineprinted=1;
			
			printf (" Hit:          %-5d     -  %5d<br>",$anfang-2,$ende-2);
                        @stemggpairs=(@stemggpairs,$anfang-4,$ende); #these are pointing here
								#	((.((    )).))
			$stemggpairfound=1;			#       ^            ^
								#       |            |
		    }
		}
                pos($str)=pos($str)-2;
            }        
        }
        pos($str)=$anfang-2;
    }
    print "<b>StemGGpair:</b>    none<br>" if ($stemggpairfound==0); 
    @sequ=();
    $str='';
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
                	print "<br><b>Pr.A1 bin.site:</b>start  -  mismatch  - seq<br>" if ($prota1startingline==0);
			$prota1startingline=1;
			my $pa1seq=substr($SEQUENCECHECKED,$wh,38);;			
			printf ("               %-5d        %-2d       %s<br>",$wh,$pa1mismatches,$pa1seq);
			
        	}
	}
	print "<br><br><b>Pr.A1 bin.site:</b>none<br>" if ($prota1startingline==0);
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


sub tRNA {
	#Looking for tRNAs using tRNAscan-SE

	$answertrnascan=`$TRNASCANFOLDER/tRNAscan-SE -Q -y -f $TEMPDIR/$job.trnascanout $TEMPDIR/$job.seq`;
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
		print "<b>tRNAscan Results:</b>         none<br>";
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
		
		
		print "<br><b>Length:</b>        $SEQUENCELENGTH";
	        print "     *some information is only available up to $MAXFOLDINGLEN nt\n" if ($SEQUENCELENGTH >$MAXFOLDINGLEN); 
	        print "<br><b>Origin:</b>        $dnarna<br>";
}


# sub checkstems {
# 		#evaluate the structure
# 		$fldklauf=0; #Klammern auf pro stem
# 		$fldstemsauf=0;
# 		$fldklzu=0;
# 		$fldstemszu=0;

# 		for ($i=0;$i<=@structure-1;$i++) {
# 		    $fldklauf++  if ($structure[$i] eq '(');
# 		    $fldstemsauf++ if ($structure[$i] eq ')' && $fldklauf>=5);
# 		    $fldklauf=0 if ($structure[$i] eq ')');
		    
# 		    $fldklzu++  if ($structure[$i] eq ')');
# 		    $fldstemszu++ if ($structure[$i] eq '(' && $fldklzu>=5);
# 		    $fldklzu=0 if ($structure[$i] eq '('); 
# 		} 
# 		#Now a correction for Stemzu!
# 		$fldstemszu++ if ($fldklzu>=5);
# 		print "<b>Energy:</b>        $energy kcal/mol<br>";
# 		print "<b>Stems:</b>         $fldstemsauf Stem-Structure/s<br>" if ($fldstemsauf==$fldstemszu);
# 		print "<b>Stems:</b>         $fldstemsauf-$fldstemszu Stem-Structure/s<br>" if ($fldstemsauf<$fldstemszu);
# 		print "<b>Stems:</b>         $fldstemszu-$fldstemsauf Stem-Structure/s<br>" if ($fldstemsauf>$fldstemszu);
# 		if ($fldstemsauf != 0 && $fldstemszu != 0) {
# 			if ((2*@structure)/($fldstemsauf+$fldstemszu)<60) { #The last number determines when the structure is interesting! 60 means at least on Stem in each 60 nt!
# 				print "<br>The sequence seems to contain a lot of secondary structure. If the RNA structure search below<br>does not find a result, it might be interesting to have a closer look at the structures.<br>";
# 				print 'You might find it useful to look in the book <br><br>"RNA Motifs and Regulatory Elements"<br>Thomas Dandekar (Ed.)<br>Published by Springer<br>ISBN 3-540-41701<br>';
# 			}
# 		}
# 		#An message if it is highly structured that it could be an rRNA
# 		if (($fldstemauf+$fldstemzu)/2>=10) {
# 			print "<br>          Highly structured RNA, could this be a ribosomal RNA ?<br>";
# 		}
# 		#Okay, that's it for the stem detection!
# }

sub checkstems {
    # newer version of checkstens as previously it was more simpler looking for >5 brackets
    # counts stems with paiting not just consecutive brackets
    # also keeps stacks while counting 
    my @pairing = ();
    my %visited;
    my @stack;

    my ($fldstemsauf, $fldstemszu) = (0, 0);

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

    # Scanning for full stems with small bulge tolerance
    for (my $i = 0; $i < @structure; $i++) {
        next if $visited{$i};
        next unless defined $pairing[$i];
        next if $pairing[$i] < $i;  # avoid counting the same pair twice

        my $left = $i;
        my $right = $pairing[$i];
        my $count = 1;

        # Store all positions involved to mark as visited later
        my @positions = ($left, $right);

        # keeping 1-nt bulges
        my ($l, $r) = ($left, $right);
        while (1) {
            my $found = 0;
            foreach my $gap_left (0, 1) {
                foreach my $gap_right (0, 1) {
                    my $nl = $l + 1 + $gap_left;
                    my $nr = $r - 1 - $gap_right;
                    next if $nl >= $nr;
                    if (defined $pairing[$nl] && $pairing[$nl] == $nr && !$visited{$nl} && !$visited{$nr}) {
                        $count++;
                        $l = $nl;
                        $r = $nr;
                        push @positions, $nl, $nr;
                        $found = 1;
                        last;
                    }
                }
                last if $found;
            }
            last unless $found;
        }

        # 5 or >5 pairing 
        if ($count >= 5) {
            $fldstemsauf++;
            $fldstemszu++;
            $visited{$_} = 1 for @positions;
        }
    }

    # print findings
    print "<b>Energy:</b> $energy kcal/mol<br>";
    if ($fldstemsauf == $fldstemszu) {
        print "<b>Stems:</b> $fldstemsauf Stem-Structure/s<br>";
    } else {
        print "<b>Stems:</b> $fldstemsauf-$fldstemszu Stem-Structure/s<br>";
    }

    my $avg_spacing = (2 * scalar @structure) / ($fldstemsauf + $fldstemszu);
    my $stem_count = ($fldstemsauf + $fldstemszu) / 2;

    # if many stems are detected
    if ($avg_spacing < 80 && $stem_count >= 20) {
        print "<b>Highly structured RNA — could this be an rRNA or similar?</b><br>";
    }

    # less structure but can be still important 
    elsif ($avg_spacing < 120 && $stem_count >= 10) {
        print "The RNA seems structurally interesting, could be a regulatory RNA.<br>";
        print "You might find it useful to look in the book: <br>RNA Motifs and Regulatory Elements<br>Thomas Dandekar (Ed.)<br>Published by Springer<br>ISBN 3-540-41701<br>";
    }
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
			if ($? != 0) { ##added by AA, error handling
    			return ("RNAfold failed", 0);
			}
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

		##alright we shall undo this change and return after the lines below AA

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
    
	print "<br><b>Pred. Protein<sup>1</sup>:</b>";
    if (defined $cpc){
        print "Protein is predicted from CPC2 output.";
    }
	if (@predprot>1){
		print "<br>";	
		for ($count1=1;$count1<=@predprot;$count1++) {
			print $predprot[$count1-1];
			print "<br>" if ($count1%120==0);
		}
	print "<br>";
	}
	else {
		print " none";
	}
}

# microrna search should be full length and should show potential micrornas with a warning. 
sub RNAMOTIF {
    my $tblout_file = "$TEMPDIR/$job.tblout";  # Table format output
	my $output_file = "$TEMPDIR/$job.out";     # Full verbose output

	my $cmd = "$CMSCAN/cmscan -E 0.001 --tblout $tblout_file -o $output_file $RFAM $TEMPDIR/$job.seq > /dev/null 2>&1";
	system($cmd);

	my $format = "%-12s %-12s %-12s %-6s %-8s %-10s %-20s\n";

	my $found = 0;
    @rna_motif =();
    my @results;
    my $found = 0;
    my $i = 0;

    print "<link rel='stylesheet' href='/css/fornac.css'>\n";
    print "<script src='/js/d3.v3.min.js'></script>\n";
    print "<script src='/js/fornac.js'></script>\n";

	print "<pre>\n";
	print "<b>RNA motif search:</b><br>\n";

	
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

		my $row = sprintf($format, $match, "$family_link     ", $from, $to, $score, $e_value, $description);
        # adding forna visual
       
        my $div_id = "rna_ss_$i";
        my $forna_html = "";

        $forna_html .= "<button onclick=\"toggleStructure$i()\">View Structure</button>\n";
        $forna_html .= "<div id='$div_id' style='width: 650px; height: 650px; display: none; margin-top: 10px;'></div>\n";
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
        $forna_html .= "        drawBackground: false\n";
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

        push @results, $row . $forna_html;

		$found = 1;
        $i++;

        print "DEBUG: MOTIF: @rna_motif";
	}
	
	close $fh_tbl;

	
	if ($found) {
		# Print header only if results exist
		printf $format, "Match", "Family", "From", "To", "Score", "E-Value", "Description";
		print "-" x 80, "\n";  # Simple separator

		# Print stored results
		print @results;

	} else {
		print "No motif recognized\n";
    }
	print "</pre>\n";
}

sub CPC2 {

	my $cpc_input = "$TEMPDIR/$job.dna";
	my $cpc_output = "$TEMPDIR/$job.cpc2";

	my $RUN_CPC="python3 $CPC/CPC2.py --ORF -r -i $cpc_input -o $cpc_output";
	my $exit_code = system($RUN_CPC);
	if ($exit_code != 0) {
    	print "CPC2 execution failed with exit code: $exit_code\n";
	}

	print "<b>Checking coding potential:</b>\n";

	open(my $fh_cpc2, "<", "$cpc_output.txt") or die "Cannot open CPC2 result $cpc_output: $!";
	my @results;
	my $found = 0;

	my $format = "%-10s %-18s %-15s %-10s %-10s %-10s %-15s %-10s %-10s\n";

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
		push @results, sprintf($format, $id, $transcript_length, $peptide_length, $fickett_score, $pI, $orf_integrity, $orf_start, $coding_probability, $label);

        $found = 1;
    }
    close $fh_cpc2;

    if ($found) {
        printf $format, "ID", "Transcript Length", "Peptide Length", "Fickett", "pI", "ORF", "ORF Start", "Coding Prob.", "Label";
        print "-" x 110, "\n";
        print @results;
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
        &scan_structured_regions;
    }
    return \%transcripts;
}

sub microRNA {

		my $mirbase_output = "$TEMPDIR/$job.mirtbl";
		my $mirbase_out    = "$TEMPDIR/$job.mir";

		my $mirna_search = "$HMMER/nhmmer --rna --watson -Z 3.73 --tblout $mirbase_output -o $mirbase_out $TEMPDIR/$job.seq $MIRBASE/hairpin.fa";

		system($mirna_search);

		my $format = "%-18s %-6s %-6s %-10s %-8s %-15s %-40s\n";

		my @results;

		open my $fh_tbl, '<', $mirbase_output or die "Cannot open miRNA result: $!";
		while (my $line = <$fh_tbl>) {
			next if $line =~ /^#/;
			chomp $line;
			my @columns = split(/\s+/, $line);
			next unless @columns >= 17;
			my $desc_full = join(" ", @columns[15..$#columns]);

			my ($match, $from, $to, $e_value, $score) =
    			($columns[0], $columns[7], $columns[8], $columns[12], $columns[13]);

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
			print "<b>miRNA search:</b><br>\n";
			printf $format, "Match", "From", "To", "E-Value", "Score", "Accession", "Description";
			print "-" x 120, "\n";

			# Split into human vs others (e.g., match starts with hsa- or description contains Homo sapiens)
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
				printf $format, $hit->{match}, $hit->{from}, $hit->{to}, $hit->{e_value}, $hit->{score}, $link, $hit->{description};
			}

			print "<br><b>Total microRNA hits found:</b> $total<br>\n";
			print "Due to the number of hits, the sequence likely contains microRNA(s)<br>\n";
		} else {
			print "No regions matching a mircroRNA was found.<br></br>\n";
		}

		print "<br></br>";


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

  
    print "<b>miRanda target prediction (3' UTR only):</b><br>\n";
    printf "%-18s %-6s %-6s %-10s %-6s\n", "miRNA", "From", "To", "Energy", "Score";
    print "-" x 60 . "\n";

    my @regions;
    foreach my $line (@filtered) {
        my (undef, $mirna, $score, $energy, $start, $end) = split /\t/, $line;
        printf "%-18s %-6s %-6s %-10.2f %-6.1f\n", $mirna, $start, $end, $energy, $score;
        push @regions, [$start, $end];
    }

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
    foreach my $tid (sort keys %transcripts) {
        my $model = $transcripts{$tid};
        
        my ($gene_id) = $tid =~ /^(g\d+)\./;
        print "<b>Transcript $i</b> (ID: $tid";
        print ", gene: $gene_id" if $gene_id;
        print ")<br>";
        $i++;

        foreach my $exon (@{ $model->{exons} }) {
            my ($start, $end) = @$exon;
            my @segments = split_exon_by_cds($start, $end, @{ $model->{cds} });
            foreach my $seg (@segments) {
                my ($seg_start, $seg_end, $type) = @$seg;
                printf("Exon: %-6d - %6d     %s<br>", $seg_start, $seg_end, $type);
            }
        }

        $model->{protein} =~ s/\s//g if exists $model->{protein};

        print "<br>";
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
                    last;
                }
            }
            if ($utr_seq =~ /A{10,}/i) {
                my $tail_pos = $-[0] + $start;
                print "         PolyA tail near $tail_pos<br>";
            }
        }
    }

    print "<i>No valid UTRs inferred.</i><br>" unless @utr;

    return (\@new5primeutr, \@new3primeutr, \@utrprintout, \@utr);
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
            push @polyasignal, $signal_pos;
            print "PolyA signal ($motif) detected at $signal_pos<br>";
            last;
        }
    }

    if ($sequence =~ /A{10,}/g) {
        my $tail_candidate = pos($sequence);
        if ($tail_candidate >= $cds_end && $tail_candidate - $cds_end < 200) {
            $tail_pos = $tail_candidate;
            push @polyatail, $tail_pos;
            print "PolyA tail detected near $tail_pos<br>";
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

    if ($printout) {
        print "<br>";
        print "<b>Structured regions detected:</b><br>";
        printf "%-12s %-6s %-6s %-10s %-10s<br>", "Region", "From", "To", "Stems", "Energy";
        print "-" x 60 . "<br>";
        foreach my $r (@high_struct_regions) {
            printf "%-12s %-6d %-6d %-10s %-10s<br>", 
                "Region", $r->[0], $r->[1], $r->[2], $r->[3];
        }
        print "<br>** Highly structured regions found. Consider tRNA, rRNA, or ncRNA elements.<br>";
    } else {
        print "<i>No regions with significant RNA structure detected.</i><br>";
    }
}

sub normalize_transcript_features {
    my ($transcripts_ref, $exons_ref, $utr_ref, $polyasignal_ref, $mirnatarget_ref, $structured_ref) = @_;

    @$exons_ref = ();
    @$utr_ref = ();
    @$polyasignal_ref = ();
    @$mirnatarget_ref = ();
    @$structured_ref = ();

    foreach my $tid (keys %$transcripts_ref) {
        my $model = $transcripts_ref->{$tid};

        if (defined $model->{cds_start} && defined $model->{cds_end}) {
            push @$exons_ref, $model->{cds_start}, $model->{cds_end};
        }

        push @$utr_ref, @{ $model->{utr5} } if ref $model->{utr5} eq 'ARRAY';
        push @$utr_ref, @{ $model->{utr3} } if ref $model->{utr3} eq 'ARRAY';

        push @$polyasignal_ref, @{ $model->{polyAcoords} } if ref $model->{polyAcoords} eq 'ARRAY';
        push @$mirnatarget_ref, @{ $model->{mirnatargets} } if ref $model->{mirnatargets} eq 'ARRAY';
        push @$structured_ref, @{ $model->{structured_regions} } if ref $model->{structured_regions} eq 'ARRAY';
    }
    print "DEBUG 1: UTRs: @utr_ref<br>";
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
        print "<h3>RNA Structure Visualization:</h3>";
        
        print "<div style='width: 700px;'>\n";
        print "  <div id='rna_ss' style='width: 700px; height: 700px;'></div>\n";
        print "  <p style='font-size: 0.9em; color: gray; text-align: center;'> Drag to pan, scroll to zoom</p>\n";
        print "<b>Download Folding As: </b>\n";
        print "<a href='$svg_url' target='_blank'><button>SVG File</button></a>";
        print "<a href='$ps_url' target='_blank'><button>PS File</button></a>\n";
        print "</div>";

        # Include required scripts
        print "<link rel='stylesheet' href='/css/fornac.css'>\n";
        print "<script src='/js/d3.v3.min.js'></script>\n";
        print "<script src='/js/fornac.js'></script>";

        my $color_text = "";

        # Process RNA motifs (red)
        
        for (my $j = 0; $j < @utr; $j += 2) {
            my ($from, $to) = @utr[$j, $j + 1];
            for (my $i = $from; $i <= $to; $i++) {
                $color_text .= "$i:lightblue ";
            }
        }

        # Process exons (green)
        for (my $j = 0; $j < @exons; $j += 2) {
            my ($from, $to) = @exons[$j, $j + 1];
            for (my $i = $from; $i <= $to; $i++) {
                $color_text .= "$i:lightgreen ";
            }
        }
        
        # process motif
        for (my $j = 0; $j < @rna_motif; $j += 2) {
            my ($from, $to) = @rna_motif[$j, $j + 1];
            for (my $i = $from; $i <= $to; $i++) {
                $color_text .= "$i:red ";
            }
        }

        $color_text =~ s/\s+$//;  # Trim trailing space

        # Inject JavaScript block to visualize RNA
        print "<script>\n";
        print "  window.onload = function () {\n";
        print "    var container = new fornac.FornaContainer(\"#rna_ss\", { animation: false, labelInterval: 50, allowPanningAndZooming: true, drawBackground: false});\n";
        print "    var options = {\n";
        print "      structure: '$structure',\n";
        print "      sequence: '$sequence'\n";
        print "    };\n";
        print "    container.addRNA(options.structure, options);\n";

        print "    var colorText = \"$color_text\";\n";
        print "    container.addCustomColorsText(colorText);\n";
        print "  };";
        print "</script>";

    } else {
        print "<br><b>Maximum folding limit reached</b><br>";
    }
}

sub location_table {

    sub format_flat_ranges {
        my @data = @_;
        my @formatted;
        for (my $i = 0; $i < @data; $i += 2) {
            last if $i + 1 > $#data;
            push @formatted, "$data[$i] to $data[$i+1]";
        }
        return join(" ", @formatted);
    }

    # Begin styling and table
    print "<style>
        table, th, td {
            border: 1px solid black;
            border-collapse: collapse;
            padding: 8px;
        }
        th {
            background-color: #f2f2f2;
            text-align: left;
        }
    </style>";

    print "<h2>Locations of the Detected Structures:</h2>\n";
    print "<table>\n";

    # Only print each row if the array is not empty
    if (@rna_motif) {
        my $motif_str = format_flat_ranges(@rna_motif);
        print "<tr><th>RNA Motifs</th><td>$motif_str</td></tr>\n";
    }

    if (@exons) {
        my $exons_str = format_flat_ranges(@exons);
        print "<tr><th>Exons</th><td>$exons_str</td></tr>\n";
    }

    if (@utr) {
        my $utr_str = format_flat_ranges(@utr);
        print "<tr><th>UTRs</th><td>$utr_str</td></tr>\n";
    }

    print "</table>\n";
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
    # print "<br>DEBUG: EXONS: @exons<br>";
    # print "DEBUG: UTRs: @utr<br>";
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

    #print $cgi->h3("");
    print $cgi->h1("Summary");
    
    #Print out the predicted protein! reformat for a nicer looking on the page	
    &predprotein;


    #Now create the html output!!!!!!!!!!!!

	if ($SEQUENCELENGTH <= $maxcoloredseqlen) {
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
        
        
        if ($c % 120 == 0) {   
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
        if ($c % 120 == 0) {
            print "     ".$c."<br>";
        }
        if ($c %10 == 0) {
            print "";
        }
    }
    #print;
    @seq=();
    print "<br>Output not possible:<br>@nooutputpossible<br>" if (@nooutputpossible>0);
	
	&drawcoloredstructure;

	print "<br><br><b>Legend:</b><br><b>BOLD</b> marks EXONS; <u>UNDERLINED</u> marks IRE or TRANSSPLICNG hits; <i>ITALIC</i> marks putative UTRs; <font color=red>RED </font> marks SMSITES or snRNP-binding motifs;";
    print "<br><font color=blue>BLUE</font> marks AU-rich regions; <font color=fuchsia>FUCHSIA</font> with <font color=green>GREEN</font> marks Stemm-GG-Pairs; <font color=lime>LIME</font> marks PolyA-signal; <font color=purple>PURPLE</font> marks a Promotor<br>";
	
    }
	else {
		print "<b>Sequence limit reached for creating colored sequence</br>";
	}
 
    #print "<br><br>%change_fontcolor<br><br>";
}

sub drawcoloredstructure {
		
	my @seq=@structure;
    normalize_transcript_features(); 

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


    # print '<font face="monospace">';

    ########### ATTENTION !!!!! ####################
    ### Old problem the index of array starts at 0 !!! #####
    ### So for the next lines to be correct we will add #####
    ### an x to the beginning of the array !!!!!!!!!! #########
    ### this won't be printed out, but then the ##########
    ### nomenclature is okay !!!!#####################
    @seq=('x',@seq);
    ####
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

    

    #Now create the html output!!!!!!!!!!!!
	print "<br></br>";
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
        
        
        if ($c % 120 == 0) {   
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


    for ($c=1;$c<=@seq-1;$c++){
        print "$seq[$c]";
        if ($c % 120 == 0) {
            print "     ".$c."<br>";
        }
        if ($c %10 == 0) {
            print "";
        }
    }
    #print;
    @seq=();
    print "<br>Output not possible:<br>@nooutputpossible<br>" if (@nooutputpossible>0);
}

print "<br>*Pr.A1.bin.site = Protein A1 binding site<br>";

write_file("$TEMPDIR/result.txt", "done\n");

print "</main></body></html>";
close $out;