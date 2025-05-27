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


$debug=0;

#Dir-Localisations
$TRNASCANFOLDER=abs_path('../bin/tRNAscan-SE/bin');
$VIENNARNAFOLDDIR=abs_path('../bin/ViennaRNA-2.7.0/bin'); #pointing to the RNAfold dir
$CMSCAN=abs_path('../bin/tRNAscan-SE/infernal-1.1.5/src'); #CMSCAN directory alrady present in transcan
$RFAM=abs_path('../databases/rfam/Rfam.cm'); #RFAM database
$CPC=abs_path('../bin/cpc2/CPC2_standalone-1.0.1/bin'); #Coding potential calculator 2 location
$HMMER=abs_path('../bin/hmmer/bin'); #hmmer location 
$MIRBASE=abs_path('../databases/mirbase'); #MirBASE database
$MIRANDA=abs_path('../bin/miranda/bin'); #miRanda 
$INTARNA=abs_path('../bin/IntaRNA/bin'); #path to INTARNA
$AUGUSTUS=abs_path('../bin/Augustus/bin/augustus'); #augugtus 
$MAXFOLDINGLEN=5000;
$MAXFOLDINGLENUTR=5000;
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

$SEQUENCESTART = read_file("$TEMPDIR/input.txt");
$SEQUENCESTART =~ /(.*)/s;
$SEQUENCESTARTCHECKED = $1;
$SEQUENCE = $SEQUENCESTARTCHECKED;

my $params_json = read_file("$TEMPDIR/params.json");
my $params = decode_json($params_json);

# Now extract them like form params:
my $do_rnamotif = $params->{RNAmotif};
my $do_augustus = $params->{run_coding};
my $do_mirna    = $params->{mirna};
my $do_trna     = $params->{trna};
my $do_IRE      = $params->{IRE};
my $do_TRANS    = $params->{TRANS};
my $species     = $params->{species};
my $do_mirnatarget  = $params->{mirna_target};
my $dnarna          = $params->{dnarna};
my $SEQNAMECHECKED   = $params->{sequence_name};
my $SEQUENCECHECKED  = $params->{sequence_clean};
my $SEQUENCELENGTH = $params->{sequence_length};


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
      <h1>RNA Analyzer 2.0</h1>
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
	print SEQPIC ">$job\n$SEQUENCECHECKED\n"; #shall create a fasta-format sequence file !!!
	close SEQPIC;


    ####### Initializing some variables for the colored output ########
    @exons=();@transsplicing=();@ire=();@smsite=();@aurichregion=();
    @stemggpairs=();@polyasignal=();@utr=();@cds=();
    print $cgi->h2("Here are the results for JOB ID: $job with sequence name: ". CGI::escapeHTML($SEQNAMECHECKED));

    &mainrun;

    # running analysis
    &analysis;

    &drawcoloredsequence;
}


sub mainrun {
	print "<pre>";
	
		if ($do_TRANS){
			&TRANS;
		}
		
		if ($do_IRE){
			&IRE;
		}
	
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



# sub oppositestrand {
# 	$SEQUENCECHECKED=reverse ($SEQUENCECHECKED);
# 	$SEQUENCECHECKED=uc($SEQUENCECHECKED);
# 	$SEQUENCECHECKED=~s/C/g/g;
# 	$SEQUENCECHECKED=~s/G/c/g;
# 	$SEQUENCECHECKED=~s/A/u/g;
# 	$SEQUENCECHECKED=~s/U/a/g;
# 	$SEQUENCECHECKED=lc($SEQUENCECHECKED);
# }

##################################################

##calls all the new subrotines
sub analysis {
  
	chdir $TEMPDIR;	

	#Create a picture!!!!
	if (length $SEQUENCECHECKED<=$maxcoloredseqlen) { #we should not fold sequences larger than this !!!
		
		
		##############################
		&createfolding;
		##############################
		
		##############################
		&checkstems; #this has to be run soon after &createpicture (?? has it ??)
		##############################
		
	}
	else {
		print "<br><b>Length:</b>        $SEQUENCELENGTH";
		print "     *some information is only available up to $MAXFOLDINGLEN nt\n" if ($SEQUENCELENGTH >$MAXFOLDINGLEN); 
		print "<br><b>Origin:</b>        $dnarna<br>";
	
	}
	
	&ARE; #search the ARE

    &createfoldingpicture;

	print "<br><b>Catalytic RNA:</b><br>";

	&smsite;

	if (do_trna){
	&tRNA; #search for tRNA
    }

	if ($do_rnamotif){
	&RNAMOTIF; #motif search		
	}
	 
	if ($do_mirna){
	&microRNA; #microrna search suing miRbase
	}

    if ($do_mirnatarget){
	&miRNAtarget; #microRNA target but scan takes too long
    }
    
    if ($do_augustus){
	&AUGUSTUS; #gene prediction replacment for genscan
	}
	else{
		&CPC2; #coding potential
	}
	

	&csfce; #Subroutine containing the long search program for those sequences
	#WRONG POSITION!!!!!!!! MUST BE OUTSIDE THE 500 nt barrier!!!!
	
	########################################################
	# Looking for the protein A1 binding site C9E
	&protA1bisite;
	if (length $SEQUENCECHECKED<=$MAXFOLDINGLEN) { 
		&stemggpairs;
	        
	}
   print"</pre>";

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


sub createfoldingpicture{
     # Run RNAfold with input and output specified
		my $seq_file = "$TEMPDIR/$job.seq";
    	my $svg_file = "$TEMPDIR/${job}_ss.svg";
		my $ps_file = "$TEMPDIR/${job}_ss.ps";
        my $svg_url = "/tmp/jobs/job_$job/${job}_ss.svg";
        my $ps_url = "/tmp/jobs/job_$job/${job}_ss.ps";

	if ($SEQUENCELENGTH <= $MAXFOLDINGLEN) {
		

		system("$VIENNARNAFOLDDIR/RNAplot --infile=$TEMPDIR/$job.foldout -f svg --filename-full"); #clicable image so it opens in new tab
			
		##need to print the image it somewhere else## AA##
		#print "<img src='$svg_path' width='800' height='650' alt='RNA Structure'>";
		print "<br><b>Folding:</br></b>";
		print "<br><img src= '$svg_url' width='650' height='650' alt='RNA Structure'</br>";
		print "\n<b>Download As: </b>";
		print "<a href='$svg_url' target='_blank'><button>SVG File</button></a>";
		print "<a href='$ps_url' target='_blank'><button>PS File</button></a>\n";
	}
	else {
		print "<br></br>";
		print "<b>Maximum folding limit reached</br>";
	}
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
	if (@predprot>1){
		print "<br>";	
		for ($count1=1;$count1<=@predprot;$count1++) {
			print $predprot[$count1-1];
			print "<br>" if ($count1%120==0);
		}
	print "<br>";
	$predprotforAnDom=join('',@predprot);
 #	&runAnDom;
	}
	else {
		print " none";
	}
}

##### AA in this subroutine I need to remove all the dependency from older genscan file
##### .genscan file can be replaced by .augugtus ( basically DNA sequence)
#okay, now creating the colored RNAStructure for rev2
sub drawcoloredstructure {
		
	#my @seq=split ('',$SEQUENCECHECKED);
	my @seq=@structure;

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

# microrna search should be full length and should show potential micrornas with a warning. 
sub RNAMOTIF {
    my $tblout_file = "$TEMPDIR/$job.tblout";  # Table format output
	my $output_file = "$TEMPDIR/$job.out";     # Full verbose output

	# Run cmscan, saving outputs to files and suppressing screen output
	my $cmd = "$CMSCAN/cmscan -E 0.001 --tblout $tblout_file -o $output_file $RFAM $TEMPDIR/$job.seq > /dev/null 2>&1";
	system($cmd);

	my $format = "%-12s %-12s %-6s %-6s %-8s %-10s %-20s\n";

	my $found = 0;

	print "<pre>\n";
	print "<b>RNA motif search:</b><br>\n";

	# Read and parse the tabular output (--tblout)
	open my $fh_tbl, '<', $tblout_file or die "Cannot open RNAmotif scan file: $!";
	while (my $line = <$fh_tbl>) {
		next if $line =~ /^#/;  # Skip comments
		chomp $line;
		my @columns = split(/\s+/, $line, 18);
		next unless @columns >= 16;

		my ($match, $family, $from, $to, $score, $e_value, $description) = ($columns[0], $columns[1], $columns[7], $columns[8], $columns[14], $columns[15], $columns[17]);

		my $family_link = "<a href=\"https://rfam.org/family/$family\" target=\"_blank\">$family</a>";

		push @results, sprintf($format, $match, $family_link, $from, $to, $score, $e_value, $description);
		$found = 1;
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

	my $cpc_input = "$TEMPDIR/$job.seq";
	my $cpc_output = "$TEMPDIR/$job.cpc2";

	my $RUN_CPC="python3 $CPC/CPC2.py -i $cpc_input -o $cpc_output";
	my $exit_code = system($RUN_CPC);
	if ($exit_code != 0) {
    	print "CPC2 execution failed with exit code: $exit_code\n";
	}

	print "<b>Checking coding potential:</b>\n";

	open(my $fh_cpc2, "<", "$cpc_output.txt") or die "Cannot open CPC2 result $cpc_output: $!";
	my @results;
	my $found = 0;

	my $format = "%-10s %-18s %-15s %-10s %-10s %-10s %-15s %-10s\n";

	# Read file line by line
	while (my $line = <$fh_cpc2>) {
		next if $line =~ /^#/;  # Skip comment/header lines
		chomp $line;

		my @columns = split(/\t/, $line);  # Split by tab
		next unless @columns >= 7;  # Ensure enough columns exist

		# Extract relevant fields
		my ($id, $transcript_length, $peptide_length, $fickett_score, $pI, $orf_integrity, $coding_probability, $label) = 
			($columns[0], $columns[1], $columns[2], $columns[3], $columns[4], $columns[5], $columns[6], $columns[7]);

		# Format numeric values with 2 decimal places
		$fickett_score = sprintf("%.2f", $fickett_score);
		$pI = sprintf("%.2f", $pI);
		$coding_probability = sprintf("%.6f", $coding_probability);  # Keep precision for probability

		# Store formatted row
		push @results, sprintf($format, $id, $transcript_length, $peptide_length, $fickett_score, $pI, $orf_integrity, $coding_probability, $label);

		$found = 1;
	}
	close $fh_cpc2;

	# Print results in table format
	if ($found) {
		printf $format, "ID", "Transcript Length", "Peptide Length", "Fickett", "pI", "ORF", "Coding Prob.", "Label";
		print "-" x 100, "\n";  # Separator
		print @results;
	} else {
		print "No result in CPC2 output\n";
	}


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

## micRNA scan using miRanda
## slow, takes over 2 minutes to scan
## implemented INTARna for this purpose and then removing overlaping regions so we can identify regions
sub miRNAtarget {
    my $mitar_input = "$TEMPDIR/$job.seq";  # Target sequence FASTA (input)
    my $miranda_out = "$TEMPDIR/$job.miranda.tsv";  # Parsed output from wrapper
    my $mirna_db = "$MIRBASE/mature.fa";
    my $raw_out = "$TEMPDIR/$job.miranda.out";

    # Run the wrapper
    my $cmd = "python3 $MIRANDA/miranda_wrapper.py --parsed_out $miranda_out --miranda_bin $MIRANDA/miranda --tmpdir $TEMPDIR $mirna_db $mitar_input $raw_out";

    my $exit_code = system($cmd);
    if ($exit_code != 0) {
        print "miRanda wrapper execution failed with exit code: $exit_code\n";
        return ();
    }

    # Read and parse output TSV
    open(my $fh, "<", $miranda_out) or die "Can't open $miranda_out: $!";
    my $header = <$fh>;  # skip header

    my @lines = <$fh>;
    chomp @lines;
    close($fh);

    # Sort and filter overlapping query regions by energy (lower = better)
    my @sorted = sort {
        (split /\t/, $a)[4] <=> (split /\t/, $b)[4] ||  # Query_Start
        (split /\t/, $a)[5] <=> (split /\t/, $b)[5] ||  # Query_End
        (split /\t/, $a)[3] <=> (split /\t/, $b)[3]     # Energy
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

    # Output HTML table
    print "<b>miRanda target prediction:</b><br>\n";
    printf "%-18s %-6s %-6s %-10s %-10s\n", "miRNA", "From", "To", "Energy", "Query";
    print "-" x 60 . "\n";

    my @regions;
    foreach my $line (@filtered) {
        my ($query, $mirna, $score, $energy, $start, $end) = split /\t/, $line;
        printf "%-18s %-6s %-6s %-10.2f %-10s\n", $mirna, $start, $end, $energy, $query;
        push @regions, [$start, $end];
    }

    return @regions;
}





########################
##augustus replacing the old genscan; need to check the whole sub for errors
sub AUGUSTUS {
    (my $dna_sequence = $SEQUENCECHECKED) =~ tr/uU/tT/;
    my ($species) = @_;
    $species ||= "human";
    my %utr_supported_species = map { $_ => 1 } qw(human);  # what is this doing?

    my $output_gff = "$TEMPDIR/$job.augustus";
    my $input_dna = "$TEMPDIR/$job.dna.fa";

    @predprot = ();
    @exons = ();
    my $found_exon = 0;
    my @exon_lines = ();
    my $protein_seq = '';
    my $capturing_protein = 0;
    my @protein_lines = ();

    my ($tss_pos, $tts_pos, $cds_start, $cds_end, $strand);
    my @cds_ranges = ();
    my @exon_ranges = ();

    # Write input sequence
    open(my $fh, '>', $input_dna) or die "Cannot write the input file for augustus: $!";
    print $fh ">$job\n$dna_sequence\n";
    close($fh);

    my $augustus_cmd = "$AUGUSTUS --softmasking=0 --protein=on --UTR=on --species=$species $input_dna > $output_gff 2>&1";
    system($augustus_cmd) == 0 or die "Failed to run AUGUSTUS: $!";

    # Parse Augustus output
    open(my $GFF, '<', $output_gff) or die "Can't open AUGUSTUS result: $!";
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
                next;
            }
            push @protein_lines, $line;
            next;
        }

        next if $line =~ /^#/;
        my @fields = split("\t", $line);
        next unless @fields >= 9;

        my ($seqid, $source, $type, $start, $end, $score, $this_strand, $phase, $attributes) = @fields;

        $strand = $this_strand if defined $this_strand and $type eq "gene";

        if ($type eq "tss") {
            $tss_pos = $start;
        } elsif ($type eq "tts") {
            $tts_pos = $end;
        } elsif ($type eq "CDS") {
            push @cds_ranges, [$start, $end];
            $cds_start = $start unless defined $cds_start;
            $cds_end = $end;
        } elsif ($type eq "exon") {
            push @exon_ranges, [$start, $end];
            $found_exon = 1;
        }
    }
    close $GFF;

    my @cds = map { @$_ } @cds_ranges;
    
    # Function to split exon regions
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

        if ($cursor <= $exon_end) {
            push @segments, [$cursor, $exon_end, 'noncoding'];
        }

        return @segments;
    }

    # Final exon processing
    if ($found_exon) {
        print "<b>Exons<sup>1</sup>:</b>        start  -   end       type<br>";
        foreach my $exon (@exon_ranges) {
            my ($start, $end) = @$exon;
            my @segments = split_exon_by_cds($start, $end, @cds_ranges);
            foreach my $seg (@segments) {
                my ($seg_start, $seg_end, $type) = @$seg;
                push @exons, ($seg_start, $seg_end);
                push @exon_lines, sprintf("<br> Exon:         %-6d - %6d     %s<br>", $seg_start, $seg_end, $type);
            }
        }
        print "$_<br>" for @exon_lines;

        $protein_seq = join('', @protein_lines);
        $protein_seq =~ s/\s//g;
        @predprot = split('', $protein_seq);

        infer_and_print_UTRs($tss_pos, $tts_pos, $cds_start, $cds_end, $strand, length($dna_sequence));

        if (@new3primeutr) {
            print "<br><i>Checking predicted 3' UTR for polyA signals...</i><br>";
            for (my $i = 0; $i <= $#new3primeutr; $i += 2) {
                my $utr_start = $new3primeutr[$i];
                my $utr_end   = $new3primeutr[$i + 1];
                scanPolyAinPredictedUTR($SEQUENCECHECKED, $utr_start, $utr_end);
            }
        } else {
            print "<i>No 3' UTR predicted by AUGUSTUS. Attempting to infer from polyA signal...</i><br>";
            refineUTRwithPolyA($SEQUENCECHECKED, $cds_end, $strand, length($SEQUENCECHECKED));
        }
    } else {
        print "<b>Exons<sup>1</sup>:</b>        none<br>";
        &scan_structured_regions();  # fallback RNA scan
    }

    print "<br>";
    return (@exons, @predprot, @cds);  # Return both exon coords and predicted protein
}


sub infer_and_print_UTRs {
    my ($tss_pos, $tts_pos, $cds_start, $cds_end, $strand, $seq_length) = @_;

    our @new5primeutr = ();
    our @new3primeutr = ();
    our @utrprintout  = ();
    our @utr          = ();

    # UTR prediction from augustus found!
    if (defined $tss_pos || defined $tts_pos) {
        print "<i>UTR site(s) detected from AUGUSTUS output.</i><br>";

        if (defined $tss_pos && defined $cds_start && $tss_pos < $cds_start) {
            push @new5primeutr, $tss_pos, $cds_start - 1;
            push @utrprintout, 5, $tss_pos, $cds_start - 1;
            # print " Predicted 5' UTR: $tss_pos - " . ($cds_start - 1) . "<br>";
        } else {
            # print " No 5' UTR predicted (CDS starts at or before TSS)<br>";
        }

        if (defined $tts_pos && defined $cds_end && $cds_end < $tts_pos) {
            push @new3primeutr, $cds_end + 1, $tts_pos;
            push @utrprintout, 3, $cds_end + 1, $tts_pos;
            # print " Predicted 3' UTR: " . ($cds_end + 1) . " - $tts_pos<br>";
        } else {
            # print " No 3' UTR predicted (CDS ends at or after TTS)<br>";
        }

	# no UTR prediction from augustus

    } else {
        print "<i>No UTR sites detected. Inferring UTRs from CDS and strand...</i><br>";
        if ($strand eq '+') {
            if (defined $cds_start && $cds_start > 1) {
                push @new5primeutr, 1, $cds_start - 1;
                push @utrprintout, 5, 1, $cds_start - 1;
                print " Inferred 5' UTR: 1 - " . ($cds_start - 1) . "<br>";
            }
            if (defined $cds_end && $cds_end < $seq_length) {
                push @new3primeutr, $cds_end + 1, $seq_length;
                push @utrprintout, 3, $cds_end + 1, $seq_length;
                print " Inferred 3' UTR: " . ($cds_end + 1) . " - $seq_length<br>";
            }
        } elsif ($strand eq '-') {
            if (defined $cds_end && $cds_end < $seq_length) {
                push @new5primeutr, $cds_end + 1, $seq_length;
                push @utrprintout, 5, $cds_end + 1, $seq_length;
                print " Inferred 5' UTR: $seq_length - " . ($cds_end + 1) . "<br>";
            }
            if (defined $cds_start && $cds_start > 1) {
                push @new3primeutr, 1, $cds_start - 1;
                push @utrprintout, 3, 1, $cds_start - 1;
                print " Inferred 3' UTR: 1 - " . ($cds_start - 1) . "<br>";
            }
        }
    }

    @utr = (@utr, @new5primeutr, @new3primeutr);
    @utr = sort { $a <=> $b } @utr;

    ### Optional: stability folding like calcUTR ###
    print "<b>UTR:</b>           start  -   end   -  stems - energy<br>";
    foreach my $i (0 .. $#utrprintout / 3) {
        my $type  = $utrprintout[$i*3];
        my $start = $utrprintout[$i*3+1];
        my $end   = $utrprintout[$i*3+2];

        printf (" %d'            %-6d - %6d", $type, $start, $end);

        my $utr_seq = substr($SEQUENCECHECKED, $start, $end - $start + 1);
        my @returnout = &checkstemsonly($utr_seq, 1);
        print "       $returnout[0]" if ($returnout[0] == $returnout[1]);
        print "       $returnout[0]-$returnout[1]" if ($returnout[0] != $returnout[1]);
        print "     $returnout[2]<br>" if ($returnout[2] != 1);
    }

    # Return same structure as calcUTR expected
    return (\@new5primeutr, \@new3primeutr, \@utrprintout, \@utr);
}

sub scanPolyAinPredictedUTR {
    my ($sequence, $utr_start, $utr_end) = @_;

    # Ensure DNA format
    $sequence =~ tr/uU/tT/;

    # Extract UTR sequence
    my $utr_seq = substr($sequence, $utr_start - 1, $utr_end - $utr_start + 1);

    my @motifs = qw(AATAAA ATTAAA TATAAA AAGAAA AATATA AGTAAA);
    my $signal_found = 0;
    my $tail_found = 0;

    # Scan for polyA signal motifs
    foreach my $motif (@motifs) {
        if ($utr_seq =~ /$motif/i) {
            my $rel_pos = $-[0];  # relative to UTR
            my $abs_pos = $utr_start + $rel_pos;
            print "PolyA signal $motif detected at position $abs_pos<br>";
            $signal_found = 1;
            last;
        }
    }

    # Scan for polyA tail
    if ($utr_seq =~ /A{10,}/i) {
        my $tail_pos = $-[0] + $utr_start;
        print "PolyA tail detected near position $tail_pos<br>";
        $tail_found = 1;
    }

    unless ($signal_found || $tail_found) {
        print "No polyA signal or tail detected in 3' UTR ($utr_start - $utr_end)<br>";
    }
}


sub refineUTRwithPolyA {
    
	my ($sequence, $cds_end, $strand, $seq_length) = @_;
	$sequence =~ tr/uU/tT/;

    our @new3primeutr = ();
    our @utrprintout  = ();
    our @utr          = ();
    our @polyasignal  = ();
    our @polyatail    = ();

    my @motifs = qw(AATAAA ATTAAA TATAAA AAGAAA AATATA AGTAAA);
    my $window_start = $cds_end + 1;
    my $window_end = $seq_length;
    my $utr_seq = substr($sequence, $window_start - 1, $window_end - $window_start + 1);

    my $signal_pos = -1;
    my $signal_motif = "";
    my $tail_pos = -1;

    # Scan for signal in last ~200bp
    foreach my $motif (@motifs) {
        if ($utr_seq =~ /$motif/i) {
            $signal_pos = $-[0] + $window_start;
            $signal_motif = $motif;
            push @polyasignal, $signal_pos;
            print "PolyA signal ($motif) detected at $signal_pos<br>";
            last;
        }
    }

    # Scan for A-rich tail
    if ($sequence =~ /A{10,}/g) {
        my $tail_candidate = pos($sequence);
        if ($tail_candidate >= $cds_end && $tail_candidate - $cds_end < 200) {
            $tail_pos = $tail_candidate;
            push @polyatail, $tail_pos;
            print "PolyA tail detected near $tail_pos<br>";
        }
    }

    # If either was found, define 3' UTR
    if ($signal_pos > 0 || $tail_pos > 0) {
        my $utr_start = $cds_end + 1;
        my $utr_end = $tail_pos > 0 ? $tail_pos : ($signal_pos + 20);  # buffer

        push @new3primeutr, $utr_start, $utr_end;
        push @utrprintout, 3, $utr_start, $utr_end;
        push @utr, $utr_start, $utr_end;

        print "Inferred 3' UTR based on polyA: $utr_start - $utr_end<br>";
    } else {
        print "No strong polyA signal/tail detected in 3' region.<br>";
    }

    return (\@new3primeutr, \@polyasignal, \@polyatail);
}

sub scan_structured_regions {
    my $count = 0;
    my $len = length $SEQUENCECHECKED;
    my $printout = 0;
    my $last = 0;

    print "<b>Region:        From     To    Stems   Energy Remark</b><br>";

    while ($count <= $len) {
        my ($query, @answer);
        if ($count + 150 < $len) {
            $query = substr($SEQUENCECHECKED, $count, 150);
            printf("<br> Pos:       %6d - %6d", $count, $count + 150);
        } else {
            $query = substr($SEQUENCECHECKED, $count, $len - $count);
            printf("<br> Pos:       %6d - %6d", $count, $len);
            $last = 1;
        }

        @answer = &checkstemsonly($query, 1);
        print "      $answer[0]" if $answer[0] == $answer[1];
        print "    $answer[0]-$answer[1]" if $answer[0] != $answer[1];
        print "    $answer[2]";

        if ($answer[0] >= 3 && $answer[1] >= 3) {
            print "   **<br>";
            $printout = 1;
        } else {
            print "<br>";
        }

        last if $last;
        $count += 100;
    }

    if ($printout) {
        print "** Three or more stem loops in this region! This is a highly structured region. <br>Please check whether tRNA, rRNA or another highly structured RNA is encoded here!<br>";
    } else {
        print "<br>";
    }
}



print "<br>*Pr.A1.bin.site = Protein A1 binding site<br>";



write_file("$TEMPDIR/result.txt", "done\n");

print "</main></body></html>";
close $out;