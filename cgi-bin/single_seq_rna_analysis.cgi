#!/usr/bin/perl -I /var/www/rnaanalyzer/cgi-bin/

use CGI;
use RNASERVER::TRANS2;
use RNASERVER::IRE;
use Bio::Tools::Genscan;
use Cwd;

$rv=CGI::new();

$debug=0;

#Dir-Localisations
$TEMPPICSDIR='/var/www/rnaanalyzer/session/';
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
print '<body text="#000000" bgcolor="#C0C0C0">';
print '<font face="monospace">';

$errorsyetprintout=0; 
#indicates that an error message has not yet been print out! Will be set to 1 if that has happened!

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

	for ($seqcountfs=1;$seqcountfs<@FASTASEQ;$seqcountfs++){
		$SEQUENCE=$FASTASEQ[$seqcountfs];
		$SEQUENCE='>'.$SEQUENCE;
		$job=&jobnumber; 
                #creates a unique number for job! herewith can manage the deal with the picture!!!!
		&startproggi;
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
	@exons=();@transsplicing=();@ire=();@smsite=();@aurichregion=();
        @stemggpairs=();@polyasignal=();@utr=();@promotor=();

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
#print "<hr>References:<br><sup>1</sup>Burge, C. and Karlin, S. (1997) Prediction of complete gene structures in human genomic DNA. J. Mol. Biol. 268, 78-94 ";
#print "<br><sup>2</sup>Lowe, T.M. and Eddy, S.R. (1997) tRNAscan-SE: A program for improved detection of transfer RNA genes in genomic sequence, Nucl. Acids Res., 25, 955-964.<br>";

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
			}
			else {
				$hits=pop @transcionareturnvalues;
				for ($count=0;$count<$hits;$count++){
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
				}
			}
			#this checkes the c. elegans consensus !!
			@transcelegansvalues=RNASERVER::TRANS2::celegans($SEQUENCECHECKED);	
			#print "<b> <big>Putative trans-splicing C.elegans-consensus</b></big> search:";
			if (@transcelegansvalues==1){
				#print "No hit detected<br>";
				print " C. elegans:  none<br>";
			}
			else {
				$hits=pop @transcelegansvalues;
				for ($count=0;$count<$hits;$count++){
					$transcelegansvalues[$count*10+5]=uc($transcelegansvalues[$count*10+5]);
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

#############################################################
#	$SEQUENCE=~/^(>[ ]?([A-Za-z0-9:\|]+))/;
#	$removefastatag=$1;
#	$SEQNAMECHECKED=$2 if ($2 ne 'NOFASTA');
#	#$SEQUENCE=~s/$removefastatag// if ($2 ne 'NOFASTA');
# commented by liang, instead by next part
# because old program here is buggy, can not get rid of seqence name completely

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
# commented out by liang
#	open (SEQPIC,">$TEMPDIR/$job.seq"); #don't know if this works!!	
#	print SEQPIC ">$job\n$SEQUENCECHECKED\n"; #shall create a fasta-format sequence file !!!
#	close SEQPIC;
  
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
#		print '<img src="'."/cgi-bin/get_result_picture.pl?"."$job"."_"."$md5sum".'"><br>';
		print '<img src=./makepic.pl?seq='.$SEQUENCECHECKED.'></img>';
	        
	}
   print"</pre>";
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
##############ARE_old is peters script, contains bugs###################
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
				printf (" Pos:       %6d - %6d",$count,$count+150);
			}
			else {
				$query=substr($SEQUENCECHECKED,$count,$len-$count);
	           		printf (" Pos:       %6d - %6d",$count,$len);	
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
#		@answerout=`$VIENNAFOLDDIR/Fold $job.seq`;
		@answerout=`echo "$SEQUENCECHECKED" | $VIENNARNAFOLDDIR/RNAfold -X`;
# -X no stream stdout, requires the modified RNAfold (by liang)

                
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
#		$md5in=`md5sum $TEMPDIR/$job.seq`;
#		$md5in=~/([0-9abcdef]{32})/;
#		$md5sum=$1;
		
		$actualdir=cwd();
		chdir "$TEMPDIR";
		
		#Create the picture !!
# comment out by liang
# because I will use a seperated cgi to complete it
#		$answerconvert=system("convert "."$job".'_ss.ps'.' -crop 0x0'." $job"."_"."$md5sum".".jpg");
		#print $rv->p("Answerconvert: $answerconvert");  #3_ss.ps
#	         $answermove=system("mv $job.jpg $TEMPPICSDIR/$job.jpg");
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
#				print "<br>The sequence seems to contain a lot of secondary structure. If the RNA structure search below<br>does not find a result, it might be interesting to have a closer look at the structures.<br>";
#				print 'You might find it useful to look in the book <br><br>"RNA Motifs and Regulatory Elements"<br>Thomas Dandekar (Ed.)<br>Published by Springer<br>ISBN 3-540-41701<br>';
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

    print '<a target="_blank" href="./export_result_to_andom.cgi?'.$predprotforAnDom;
    print '">Analyze the Predicted Protein with AnDom</a><br>';
return;


# commented out by liang 
# because original is too verbose, can be simplified
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


