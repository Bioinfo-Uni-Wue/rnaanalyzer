package RNASERVER::IRE;

use strict;
use lib "/usr/local/lib/perl5/site_perl/5.18.2";
use lib "/mnt/c/Users/ama55id/Nextcloud/RNA_analyzer/rnaanalyzer/bin/ViennaRNA-2.6.4/interfaces/Perl"; 
use RNA;
$RNA::noLonelyPairs=1;	#important for correct use of the RNA-module!

sub findire {

    my $SEQUENCE=$_[0];
    my $hits=0;
    my $findire; ##ATTENTION THIS DOES NOT WORK, BUT NEVERMIND! WE'LL CHANGE THIS LATER
    my $structure;my $min_en;my $string; my $length;my $motif2;my @structure;my @string;my $guterTreffer=0;my $motif1;my $lnt0;my $lnt1;my $lnt2;my $lnt3;my $lnt4;my $lnt5;
    my $j;my $run1;my $anzguteTreffer;my @returnvalues=();my $evalpoints;
    #print $findire->p("Here are the results for your query: $SEQNAME");
    my @gen=split("",$SEQUENCE);
    my $genlen=@gen-1;  
    ##################################################
    ##################################################
    for ($j=0;$j<$genlen-10;$j++) {
        #$guterTreffer=0;
        $anzguteTreffer=0;
        $motif1=0;
        $lnt0=uc($gen[$j]);$lnt1=uc($gen[$j+1]);$lnt2=uc($gen[$j+2]);$lnt3=uc($gen[$j+3]);$lnt4=uc($gen[$j+4]);$lnt5=uc($gen[$j+5]);
        $motif1=$motif1+2  if ((($gen[$j] eq 'c') && ($gen[$j+4] eq 'g'))||(($gen[$j] eq 'u') && ($gen[$j+4] eq 'a'))||(($gen[$j] eq 'g') && ($gen[$j+4] eq 'c')));
        $motif1++          if ($gen[$j+1] ne 'g');
        $motif1++          if (($gen[$j+2] eq 'a') || ($gen[$j+2] eq 'g'));
        $motif1++          if ($gen[$j+5] ne 'g');

        #Continuing when we have more than 4 points in the loop !!
        if ($motif1>4) {		
            $motif2=0;
                    ########## Vorbereitung fuer RNAFold #######################
                    #$string=substr($gensq,$j-25,44);
            $string=substr($SEQUENCE,$j-20,40);
                    ########## Achtung ab jetzt RNAFold Teil ###################
            $string = uc($string);
            $length = length($string);
                    #printf("length = %d\n", $length); #if ($istty);                 
            $structure=$string; # wierd way to allocate space     UUGCUUUCCAACUUCAGCUACAGUGUUAGCUAAGUUUGGA
            $min_en = RNA::fold($string, $structure);
          
                    ###########################################UUGCUUUCCAACUUCAGCUACAGUGUUAGCUAAGUUUGGA#################
                    ############# Jetzt wird das Motif beschrieben #############
            @structure=split ("",$structure);
            @string=split ("",$string);
            for ($run1=6;$run1<=$length-7;$run1++){
                if (($string[$run1] eq $lnt0) && ($string[$run1+1] eq $lnt1) && ($string[$run1+2] eq $lnt2) && ($string[$run1+3] eq $lnt3) && ($string[$run1+4] eq $lnt4) && ($string[$run1+5] eq $lnt5)) {
                    if (($structure[$run1] eq '.') && ($structure[$run1+1] eq '.') && ($structure[$run1+2] eq '.') && ($structure[$run1+3] eq '.') && ($structure[$run1+4] eq '.') && ($structure[$run1+5] eq '.' && ($structure[$run1+6] eq ')'))) {
                                #Jetzt haben wir schon mal den Loop
                        if ($structure[$run1-1] eq '(' && $structure[$run1-2] eq '(' && $structure[$run1-3] eq '(' && $structure[$run1-4] eq '(' && $structure[$run1-5] eq '(' &&$structure[$run1-6] eq '.' && $string[$run1-6] eq 'C') {
                            #Or doing it more loosely!
                                $evalpoints=0;
                                $evalpoints++ if ($structure[$run1-1] eq '(');
                                $evalpoints++ if ($structure[$run1-2] eq '(');
                                $evalpoints++ if ($structure[$run1-3] eq '(');
                                $evalpoints++ if ($structure[$run1-4] eq '(');
                                $evalpoints++ if ($structure[$run1-5] eq '(');
                                $evalpoints++ if ($structure[$run1-6] eq '.'); #So we could now say that if we have 5 points, this is good enough!!!! #But we'll have to implement this further !!!
                                   ######## Dann sind schon einige Bedingungen fuer einen Treffer gegeben !! #############
                                    ### Wobei dies auch echt harte und keine wirklichen loose Kriterien sind !! ev noch 1 mismatch einbauen !!
                                    #Aber vorest bleibt es mal zum Webserver-testen so!
                            $guterTreffer=1;
                            $anzguteTreffer++;
                            push (@returnvalues,$string,$structure,$min_en,$j);
                            if ($min_en>-6) {
                                push (@returnvalues,'weak');
                            }
                            else {
                                push (@returnvalues,1);
                            }
                            $hits++;
                        }
                    }
                }
            }
        }
        
    }
    print;
    if ($guterTreffer==1) {
            push (@returnvalues,$hits);            
            return @returnvalues;            
            #print $findire->p("Possible hit detected with free energy $min_en");
        }
        else {
            return 0;            
            #print $findire->p("Results for your $SEQNAME");
            #print $findire->p("No hit detected");

        }


}

sub suboptimalfindire {

	my $gensq=$_[0];


#Novel version 4.072002 written July 2002 for the experiments of mayka sanchez (embl)

#Here at first we use a strict consensus of 



	my @subopthits=();
	##### Here are the main values!! #################
	my $loopdown=17;
	my $loopup=22;

	#cutoff values!
	my $upperstemcutoff=4;
	my $lowerstemcutoff=5;
	my $upperstemcutoffbadhit=3;
	my $lowerstemcutoffbadhit=3;
	my $hitbeforecdscutoff=200;
	my $hitaftercdscutoff=2000;
	#energycutoff=??;
	my $pathtornasuboptandrnafold='/var/www/rnaanalyzer/bin/ViennaRNA-1.5/Progs/RNAsubopt'; #Please locate your version of RNAfold / RNAsubopt

	###########################################

	$gensq="nnnnnnnnnn".$gensq."nnnnnnnnnn"; #now we append some n's in order to prevent the passing hits due to too long excicion and not uninitialized errors! But remember later !!
	
	#####################################
	#####Now starting the novel coding#########
	#####################################
	#####################################	
	#while ($gensq=~m/c[acgu]{5}([c|u]agu[g|a][a|c|u])/g){    
	while ($gensq=~m/c[acgu]{5}(([c|u]agu[g|a][a|c|u])|(ccgagcu)|(cugggc)|(ccgcgc)|(gcgccg)|(gagucg)|(gagagu))/g){
		#        bulgedc       knownloop             options from henderson et al JBC 1996
	    my $looppresent=0;
	    my $bulgedcpresent=0;
	    my $upperstempaired=0;
	    my $lowerstempaired=0;
	    my $passhit=0;
	    my @rnafoldingenergy=();
	    my @rnafoldingstruct=();
	    #@subopthits=();
	    #end of initialization
	    
	    my $looppos=(pos $gensq) - (length $1)+1; #explicitely write +1 for correction 
	    #########################################################
	    ############### Starting subopt #############################
	    #########################################################    
	    my $foldstring=substr($gensq,$looppos-$loopdown,$loopdown+$loopup); #15+20 means the formerly 15 + further 20
	    
	    my @rnafoldinganswer=`echo $foldstring | $pathtornasuboptandrnafold`;
	     #directs to the last hit in @rnafoldinganswer
	    for (my $rnafoldingstructcount=1;$rnafoldingstructcount<=@rnafoldinganswer-1;$rnafoldingstructcount++){	    
		$rnafoldinganswer[$rnafoldingstructcount]=~/([.()]+)[ ]+\(?([-0-9.]+)\)?/;
		push @rnafoldingstruct,$1;
		push @rnafoldingenergy,$2;
	    }
	    my @rnafoldingenergymax=sort {$a<=>$b} @rnafoldingenergy;
	    
	    my $rnafoldinganswerlines=@rnafoldingstruct-1;
	    for (my $suboptcount=0; $suboptcount<=$rnafoldinganswerlines; $suboptcount++){
		$looppresent=0;
		$bulgedcpresent=0;
		$upperstempaired=0;
		$lowerstempaired=0;
		#$foldstructure=$foldstring; # wierd way to allocate space       
		#$min_en = RNA::fold($foldstring, $foldstructure);
	        $passhit=0;	
		my @foldstring=split ('',$foldstring);
		my @foldstructure=split('',$rnafoldingstruct[$suboptcount]);
		
		$looppresent=1 if ($foldstructure[$loopdown-2] eq '(' && $foldstructure[$loopdown-1] eq '.' && $foldstructure[$loopdown] eq '.' &&$foldstructure[$loopdown+1] eq '.' &&$foldstructure[$loopdown+2] eq '.' &&$foldstructure[$loopdown+3] eq '.' &&$foldstructure[$loopdown+4] eq '.' && $foldstructure[$loopdown+5] eq ')');
		$bulgedcpresent=1 if ($foldstructure[$loopdown-7] eq '.');
		
		for (my $upperstemcount=1;$upperstemcount<=5;$upperstemcount++){ #in the end we will force here 4nt to be paired!
		    $upperstempaired++ if ($foldstructure[$loopdown-7+$upperstemcount] eq '('); #$loopdown-7 points to the bulged c and then we're going up the stem
		}
		
		for (my $lowerstemcount=1;$lowerstemcount<=11;$lowerstemcount++){  #we have no idea yet, how many will have to be paired
		    $lowerstempaired++ if ($foldstructure[$loopdown-7-$lowerstemcount] eq '('); #$loopdown-7 points to the bulged c and then we're going up the stem
		}
		#### Now the evaluation of the hit! ####
		
		if ($looppresent==1 && $bulgedcpresent==1) { #this is absolutely nessecary, now more detailed!!!
			if ($upperstempaired>=$upperstemcutoff && $lowerstempaired>=$lowerstemcutoff){
				$passhit=1; #indicating we have a very good hit !!!
			}
			elsif ($upperstempaired>=$upperstemcutoffbadhit && $lowerstempaired>=$lowerstemcutoffbadhit){
				$passhit=2; #indicating we have a bad hit !!!
			}
		}
		if ($passhit>0){
			
			push @subopthits,$looppos-10,$passhit,$foldstring,$rnafoldingstruct[$suboptcount],$rnafoldingenergy[$suboptcount],$upperstempaired,$lowerstempaired;
		#giving back the result: [0]=Looppos, [1]=Quality 1=good 2=bad, [2] hitstructure, [3] energy, [4]upperstempaired, [5]=lowerstempaired
		}		
	    }
	}
	return @subopthits;
	    #########################################################
	    ############### Ending  subopt #############################
	    #########################################################
}
