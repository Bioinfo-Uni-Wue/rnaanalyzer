package RNASERVER::TRANS2;

use Cwd qw(abs_path);
use strict;
use lib '/mnt/c/Users/ama55id/Nextcloud/RNA_analyzer/rnaanalyzer/bin/ViennaRNA-2.7.0/interfaces/perl';
use RNA;
use RNASERVER::COMMON;

#Attention !! Often recivied an error "Process terminated with exit code 0" But only got this error while debugging with Komodo 1.2 (beta)
# In normal use, this shouldn't happen. Maybe some return code ist missing, but once got "Indentifier too long" ?????
# Wohl eher der normale return-code, nur warum ist der 0 ?? Naja wohl egal. weil nirgends 1 zurueckgegeben wird

sub celegans {
    $RNA::noLonelyPairs=1;
    #print "\nDeb0";
    my @returnvalues=();
    my $gensq=$_[0];

    #we will nowinitializte the variables!
    my $newstem1klzu; my $gensqstore;my @gen;my $nt;my $genlen;my $idflag;my $Oldmotifsum;my $leader;my $seqscan1;my $motifsum;my $a1;my $a2;my $crit1;my $crit2;my $crit3;my $seqscan2; my $gguakorrektpresent;
    my $stem1test;my $string;my $length;my $structure;my $min_en;my $stem1string; my $stem1structure;
    my @stem1stringarray; my @stem1structurearray; my $lengthstem1; my $test; my $matchc; my $seqscan3;my $seqscan3ersatz;
    my $loop1nt; my $ucount; my $acount; my $seqscan4; my $Gremark; my $trans1; my $trans2; my $b1; my $b2;
    my $b3; my $b4; my $b5; my $seqscan6; my $seqscan7; my @struct; my $klammer2; my $punkt2;
    my $fold2pass; my $stem2ntbefore; my $stem2ntafter; my $str2; my $count;
    my $structurekorrigiert; my $klammer2korrigiert; my $str2korrigiert; my @structkorrigiert; my $punkt2korrigiert;my $gepaart2;
    my $gepaart2grenzwert; my $stem2auf; my $stem2zu; my $stem2reg; my $structkorrigiert2;
    my $stem2bulges; my $stem2bulgednt; my $stem2partonelength; my @structkorrigiert2;
    my $stem2bulges2; my $stem2bulgednt2; my $stem2partonelength2; my $stem2bulgedntgesamt; my $stem2korrigiertlength;
    my $knockout; my $klammerauf; my $twostem2present;my $alternativefoldpass2;my $twostem2intermediate;
    my $twostem2loopnt; my $stem2stems; my $klammeraufmax; my $stem2string; my $sumenergy; my $sumenergy2;
    my $seqscan12; my $seqscan13; my $energy2; my $strechC; my $strechD;
    my $point12;my $point13;my $stem3test;my $stem2structure;my @structstem3;my $klammer3;my $seqscan14;my $sumsum;
    my $leaderw;my $OldGremark;my $oldmotifsum; my $newhit;my $smsite; my $leaderseq;my $hits=0;

    $gensqstore=$gensq;              
    @gen=split ("",$gensq);          
    #########################################################
    #########################################################
    #Das davor haben wir schon !!

    #$nt=$nt+@gen;
    $genlen=@gen-1;                   
    $idflag=0;
    #Jetzt suchen wir das transsplicing motif
    $Oldmotifsum=0;$leader=0;
    for ($seqscan1=1;$seqscan1<$genlen;$seqscan1++) { 
        if ($gen[$seqscan1] eq 'g' && $gen[$seqscan1-1] eq 'g' && $gen[$seqscan1+1] eq 'u' && $gen[$seqscan1+2] eq 'a') {
	    #print"\nDeb GGUA";
		$newhit=0;        	#This is 0 when a hit has not jet been print out. the idea is that if a hit with that ggua has been prnit out, no further will be print out too because they are mainly the same !!    
            $motifsum=0;$a1=0;$a2=0;                # logo, Variablen auf 0 setzen
            $a1=1 if ($gen[$seqscan1+1] eq 'u');
            $a2=1 if ($gen[$seqscan1+2] eq 'a');    #prueft die ggau sequenz siehe consensus
            $crit1=0;$crit2=0;$crit3=0;
            if ($gen[$seqscan1+3]=~/u|c/){     #checkt Kriterium3
                $crit3=1;
            }
            if ($gen[$seqscan1+4] eq 'g'){    #checkt Kriterium2
                $crit2=1;
            }
            if ($crit3 && $crit2 && ($gen[$seqscan1+5] eq 'u')){   #und Kriterium 1, Bewertungsreihenfolge noch unklar
                $crit1=1;
            }
	    #print "\nDebPoint1";
            for ($seqscan2=$seqscan1-2;$seqscan2>=$seqscan1-11;$seqscan2=$seqscan2-1) {
                if ((($gen[$seqscan2] eq 'u') ||($gen[$seqscan2] eq 'c'))  && (($gen[$seqscan2-1] eq 'u') || ($gen[$seqscan2-1] eq 'c'))) { #die sind die Basen die sich mit dem GG paaren sollen!!!!
                    ###################### Versuch die stem1 und den loop per rnafold zu bestimmen #############################
                    $gguakorrektpresent=0;
		    #$stem1test=substr($gensq,($seqscan2-16+1),($seqscan1+16-($seqscan2-16)));
                    $stem1test=substr($gensq,($seqscan2-18+1),($seqscan1+16-($seqscan2-18)));
		            $string=$stem1test;
                    $string = uc($string);
                    
                    $length = length($string);
                    #printf("length = %d\n", $length); #if ($istty);
                    # $structure=$string;
                      # Save the original nucleotide sequence
                    ($structure, $min_en) = RNA::fold($string);
                    
                    $stem1string=$string;
                    $stem1structure=$structure;
                
		    #print "\n$stem1string\n$stem1structure";
                    ########## Versuch das ggua-Motif anhand des fold-outputs besser zu evaluieren #########
                    @stem1stringarray=split ("",$stem1string);
                    @stem1structurearray=split ("",$stem1structure);
		    $newstem1klzu=0;
		    foreach $test (@stem1structurearray) {
		    	$newstem1klzu++ if ($test eq ')');
		    }
                    $lengthstem1=@stem1structurearray-1;
                    for ($test=0;$test<=$lengthstem1;$test++){
                        if ($stem1stringarray[$test] eq 'G' && $stem1stringarray[$test+1] eq 'G' && $stem1stringarray[$test+2] eq 'U' && $stem1stringarray[$test+3] eq 'A'){
                            if ($stem1structurearray[$test] eq ')' && $stem1structurearray[$test+1] eq ')' && ($stem1structurearray[$test+2] eq ')' || $stem1structurearray[$test+3] eq ')') && $newstem1klzu>5){# && $stem1structurearray[$test-1] eq '.' && $stem1structurearray[$test-2] eq '.' && $stem1structurearray[$test-3] eq '.') {
                                $gguakorrektpresent=1;
				#print "\nGGUAkorrektpresent1 mit $newstem1klzu Klzu";
				
                            }
			    
                        }
                    }
                            
                        
                    #printf("\n minimum free energy = %6.2f kcal/mol\n", $min_en);
                    ###################### Ende des stem1   ####################################################################
                    #$stem1string=$string;$stem1structure=$structure;
                    $matchc=0;
                    #Das Ziel der naechsten Schleife ist es, den folgenden "horizontalen" paired zu finden !!!
                    #genauer den Teil mit 3 of 4 must at least pair !!!!!!!
                    for ($seqscan3=$seqscan2-2;$seqscan3>=$seqscan2-5;$seqscan3=$seqscan3-1) {
                        $seqscan3ersatz=$seqscan3;
                        

                        #for ($bulgestem1=$seqscan2-3;$bulgestem1>=$seqscan2-5;$bulgestem1--){
                            #if ($bulgestem1==$seqscan3){
                             #   $seqscan3ersatz++;
                            #}
                            $matchc++ if (($gen[$seqscan3ersatz] eq 'a' && ($gen[$seqscan1+$seqscan2-$seqscan3-1] eq 'u')));
                            $matchc++ if (($gen[$seqscan3ersatz] eq 'u' && ($gen[$seqscan1+$seqscan2-$seqscan3-1]=~/a|g/)));
                            $matchc++ if (($gen[$seqscan3ersatz] eq 'g' && ($gen[$seqscan1+$seqscan2-$seqscan3-1]=~/c|u/)));
                            $matchc++ if (($gen[$seqscan3ersatz] eq 'c' && ($gen[$seqscan1+$seqscan2-$seqscan3-1] eq 'g')));
                       # }
                    }
                    
                    $loop1nt=0;
		    #print "\nMatchc=$matchc";
                    if ($gguakorrektpresent==1) {    #wenn also ein gepaarter Stem vorliegt		    
                    #Dann suchen wir den loop davor und schauen ob ja auch brav mind 3 U vorkommen !!!
                        $ucount=0;$acount=0;
			#print "\nCountUs from $seqscan1 -2 bis $seqscan2 +1";
                        

                        for ($seqscan4=$seqscan1-2;$seqscan4>=$seqscan2+1;$seqscan4=$seqscan4-1) {
                            $ucount++ if ($gen[$seqscan4] eq 'u');
                            $acount++ if ($gen[$seqscan4] eq 'a');
                            $loop1nt++;
                            
                        }
                        

			#print "\nUcount=$ucount";
                        if ($ucount>=1){ #&& ($seqscan1-$seqscan2-$acount-$ucount<6)) 
                            #print "\nUcount passed";
				$Gremark=$seqscan1;
                            $trans1=$seqscan1+4;
                            $trans2=$seqscan2-5;      #zeigt auf den Anfang der non-obligatory-features
                            #Testen ob der Stamm verlngert werden kann!
                            

                            while ($Gremark-$trans1+$matchc+3>0 && 0<$trans1 && $trans1<$genlen && 1<$trans2 && $trans2<$genlen) {
                                $trans1++;$trans2=$trans2-1;
                                $matchc++ if (($gen[$trans1] eq 'a') && ($gen[$trans2] eq 'u'));
                                $matchc++ if (($gen[$trans1] eq 'c') && ($gen[$trans2] eq 'g'));
                                $matchc++ if (($gen[$trans1] eq 'g') && ($gen[$trans2]=~/c|u/));
                                $matchc++ if (($gen[$trans1] eq 'u') && ($gen[$trans2]=~/a|g/));
                                

                            }
                            #Testen ob das G/U U/A U/A U motif am Ende des splicing stems ist
                            $b1=0;$b2=0;$b3=0;$b4=0;$b5=0;
                            $b1++ if ($gen[$trans2]=~/g|u/);
                            $b2++ if ($gen[$trans2+1] eq 'u');
                            $b3++ if ($gen[$trans2+2]=~/u|a/);
                            $b4++ if ($gen[$trans2+3]=~/u|a/);
                            $b5++ if ($gen[$trans2+5] eq 'u');
                            $motifsum=$a1+$a2+$b1+$b2+$b3+$b4+$b5;
                           
                            #Jetzt wird gecheckt ob eine Sm-Site downstream in der Naehe ist.
			    #my $newupperlimit;
			   #$newupperlimit=80;
			   #$newupperlimit=
                            for ($seqscan6=$trans1+9;$seqscan6<=$trans1+80;$seqscan6++) {
                                #$smsitepresent=0;
                                #$alternativsmcheck=substr($gensq,($trans1+9),$trans1+80);
                                #if ($alternativsmcheck=~/(aacugg)|(aaccug)|(aauu[u]+cuuu[u]+ga)|(aauuuuuu[u]+auaa)/){
                                    #$smsitepresent=1;
                               # }
                                if ($gen[$seqscan6]=~/a|g/) {
                                    if ($gen[$seqscan6+1]=~/a|g/) {       #u exeptionally allowed
                                        $seqscan7=$seqscan6+2;
                                        while ($gen[$seqscan7] eq 'u') {
                                            $seqscan7++;
                                        }                              #diese Konstruktion erlaubt auch ein anderes nt
                                        $seqscan7++;
                                        
                                        while ($gen[$seqscan7] eq 'u') {
                                            $seqscan7++;
                                        }
                                        if ((($seqscan7-$seqscan6)>5)){ #| $smsitepresent==1) {
                                            if (($gen[$seqscan7]=~/a|g/) && ($gen[$seqscan7+1]=~/a|g/)||($gen[$seqscan7] eq 'c' && $gen[$seqscan7+1] eq 'a')) {
						    $smsite='';
                            
						    
						    #for ($xfgf=1;$xfgf<10;$xfgf++){}
						    
						    
						    #for ($newcount=$seqscan6;$newcount<=$seqscan7+1;$newcount++){
							    #	    $smsite=$smsite.$gen[$newcount];
							    # }
						    #Is there a loop between Sm-Site and transsplicing loop? Check!
                                                #print "\n-------------------------------\nSequenzausgabe: Hier drinnen sollte ein Stem sein:\n";
                                                #for ($printout=$trans1;$printout<=$seqscan6;$printout++){
                                                    #print "$gen[$printout]";
                                                #}
                                              
					     	$smsite=substr($gensq,($seqscan6),($seqscan7-$seqscan6+2));	
						################### Anpassung an RNA-Fold ######################
                                                $string=substr($gensq,($trans1+1),($seqscan6-$trans1));
                                                ################### Teil aus RNA-Fold ##########################						    
                                                $string = uc($string);
                                                $length = length($string);
                                                #printf("length = %d\n", $length); #if ($istty);						   						  
                                                $structure=$string; # wierd way to allocate space						 
                                                ($structure, $min_en) = RNA::fold($string);
                                                
                                                ######### Diese Zeile startet den Versuch, zu pruefen, ob das GGUA wirklich gepaart ist, und ob davor ein loop ist #####
                                             
                                                    
                                                #print "$string\n$structure";						  
                                                #printf("\n minimum free energy = %6.2f kcal/mol\n", $min_en);						   
                                                ################### Ende RNA-Fold ##############################
                                                ################### Entscheidungskriterien fuer foldoutput #####
                                                @struct=split ("",$structure);
                                                $klammer2=0;
                                                $punkt2=0;
                                                $fold2pass=0;
                                                $stem2ntbefore=0;
                                                $stem2ntafter=0;
                                                foreach $str2 (@struct) {
                                                    if ($str2 eq '(' || $str2 eq ')'){
                                                        $klammer2++;
                                                    }
                                                    else {
                                                        $punkt2++;
                                                    }
                                                }
                                                #print "\nGepaarte: $klammer1 und ungepaarte: $punkt1";
                                               # if ($klammer1>$punkt1){
                                                    #print "\nStamm von fold akzeptiert";
                                                    #$fold1pass=1;
                                                #}
                                                ######### Maxmimal 7nt vorher und nachher zulassen !! #########
                                                for ($count=0;$count<=9;$count++){
                                                    if ($struct[$count] eq '('){
                                                        $stem2ntbefore=1;
                                                    }
                                                }
                                                for ($count=-7;$count<=-1;$count++){
                                                    if ($struct[$count] eq ')'){
                                                        $stem2ntafter=1;
                                                    }
                                                }
                                                ######### Korrigierte Stemwerte ###############################
                                                #$structure=~/(\(\))/;
                                                $structure=~/\.*(\(.*\))/;
                                                $structurekorrigiert=$1;
                                                @structkorrigiert=split("",$structurekorrigiert);
                                                $klammer2korrigiert=0;
                                                $punkt2korrigiert=0;
                                                foreach $str2korrigiert (@structkorrigiert) {
                                                    if ($str2korrigiert eq '(' || $str2korrigiert eq ')'){
                                                        $klammer2korrigiert++;
                                                    }
                                                    else{
                                                        $punkt2korrigiert++;
                                                    }
                                                }
                                                if ($klammer2korrigiert!=0 && $punkt2korrigiert!=0) {
                                                    $gepaart2=($klammer2korrigiert/($klammer2korrigiert+$punkt2korrigiert));
                                                    $gepaart2grenzwert=0.00266*($klammer2korrigiert+$punkt2korrigiert)+0.5174;
                                                    if ($gepaart2>=$gepaart2grenzwert) {
                                                        $fold2pass=1;						   
                                                        #print "\nGepaart: $gepaart2 Grenze: $gepaart2grenzwert";
                                                    }  
                                                }
                                                ################### Testen ob der Stem aus "richtige" Propotionen hat ###########
                                                $stem2auf=0;
                                                $stem2zu=0;
                                                $stem2reg=1;
                                                foreach $test (@struct){
                                                    if ($test eq '('){
                                                        $stem2auf=1;
                                                    }
                                                    if ($test eq ')'){
                                                        $stem2zu=1;
                                                    }
                                                    if (($test eq '(')&& $stem2zu==1){
                                                        $stem2reg=0;    # hier bedeutet 1, dass stem2 nicht aus einem "richtigen" Stamm besteht, ev. aber aus 2 kurzen St.
                                                    }
                                                }
                                                if ($stem2reg==0){
                                                    #print "\nDieser Stamm ist nicht regulaer !!!!";
                                                    $fold2pass=0;
                                                }
                                                ################## bestimmte Anzahl von Bulges zulassen ####################
                                                @structkorrigiert2=@structkorrigiert;
                                                while ($structkorrigiert2[-1] eq ')' || $structkorrigiert2[-1] eq '.'){
                                                    pop @structkorrigiert2;
                                                }				
                                                $stem2bulges=0;
                                                $stem2bulgednt=0;
                                                $test=0;
                                                $stem2partonelength=@structkorrigiert2;
                                                $stem2partonelength--;
                                                for ($test=0;$test<=$stem2partonelength;$test++){
                                                    if ($structkorrigiert2[$test] eq '('){
                                                        #$test++;                            #noch wird hier der Loop als Bulge angesehen!!
                                                    }
                                                    else{
                                                        if ($structkorrigiert2[$test] eq '.'){
                                                            $stem2bulges++;
                                                            $stem2bulgednt++;
                                                            $test++;
                                                            while ($structkorrigiert[$test] eq '.') {
                                                                $stem2bulgednt++;
                                                                $test++;
                                                            }
                                                        } 	
                                                    }
                                                }
                                                #print "\nIm aufsteigenden Stem2 $stem2bulges mit insgesamt $stem2bulgednt";
                                                ################## Das gleiche jetzt fuer den anderen Teil des Stems #############
                                                @structkorrigiert2=@structkorrigiert;
                                                while ($structkorrigiert2[0] eq '(' || $structkorrigiert2[0] eq '.'){
                                                    shift @structkorrigiert2;
                                                }
                                                $stem2bulges2=0;
                                                $stem2bulgednt2=0;
                                                $stem2partonelength2=@structkorrigiert2;
                                                $stem2partonelength2--;
                                                for ($test=0;$test<=$stem2partonelength2;$test++){
                                                    if ($structkorrigiert2[$test] eq '.'){
                                                        $stem2bulges2++;
                                                        $stem2bulgednt2++;
                                                        $test++;
                                                        while ($structkorrigiert2[$test] eq '.'){
                                                            $stem2bulgednt2++;
                                                            $test++;
                                                        }
                                                    }
                                                }
                                                #print"\nFuer den absteigenden Teil findet man: $stem2bulges2 mit $stem2bulgednt2 bulged nt";
                                                $stem2bulgedntgesamt=$stem2bulgednt+$stem2bulgednt2;
                                                ################## Bulged Grenzwerte ##########################
                                                $stem2korrigiertlength=@structkorrigiert;
                                                if (($stem2korrigiertlength <=20 && ($stem2bulgedntgesamt)>=3) || ($stem2korrigiertlength>20 && $stem2korrigiertlength<=40 && $stem2bulgedntgesamt >=7) || ($stem2korrigiertlength>40 && $stem2bulgedntgesamt>=19) ){
                                                    $knockout=0;
                                                    #print "\nKnockout !!!!!!!!";
                                                }
                                                else {
                                                    $knockout=1;
                                                }
                                                ################## Knockout-Zeile #############################
                                                if ($stem2reg==0 || $klammer2korrigiert<4 || $gepaart2<=$gepaart2grenzwert || $stem2ntbefore==0 || $stem2ntafter==0 || $knockout==0) {
                                                    $fold2pass=0;
                                                    #print "\nFold2pass nicht passiert";
                                                    #print "\nKlammern: $klammer2korrigiert";
                                                }
                                                if ($fold2pass==1){
                                                    #print "\nFold2pass ok";
                                                }
                                             
                                                ################### Ende dieser Entscheidungskriterien ########
                                                ################### Dieser Abschnitt soll pruefen, ob es eventuell 2 stems im Abschnitt gibt ##############
                                                #nur wenn foldpass negativ
                                                $klammerauf=0;
                                                $twostem2present=0;
                                                $test=0;
                                                $alternativefoldpass2=0;
                                                $twostem2intermediate=1;
                                                $twostem2loopnt=0;
                                                for ($stem2stems=1;$stem2stems<=2;$stem2stems++){
                                                    while ($structkorrigiert[$test] eq '.'){
                                                        $test++;
                                                        $twostem2intermediate++;
                                                    }
                                                    while ($structkorrigiert[$test] eq '('){
                                                        $klammerauf++;
                                                        $test++;
                                                    }
                                                    if ($structkorrigiert[$test] eq '.'){ #ev. hier besser while und dann if und dann wieder while
                                                        $test++;
                                                        while ($structkorrigiert[$test] eq '('){
                                                            $test++;
                                                            $klammerauf++;          #Diese Konstrukition zaehlt die auf-Klammern und laesst einen Bulge zu!!
                                                        }
                                                    }
                                                    $twostem2loopnt++;
                                                    while ($structkorrigiert[$test] eq '.') {
                                                        $test++;
                                                        $twostem2loopnt++;
                                                    }
                                                    $klammeraufmax=$klammerauf;
                                                    while ($structkorrigiert[$test] eq ')'){
                                                        $test++;
                                                        $klammerauf--;
                                                    }
                                                    if ($structkorrigiert[$test] eq '.'){
                                                        $test++;
                                                    }
                                                    while ($structkorrigiert[$test] eq ')'){
                                                        $test++;
                                                        $klammerauf--;
                                                    }
                                                    if ($klammerauf==0) {
                                                        $twostem2present++;
                                                    }
                                                }
                                                #print "\ntwostem2loopnt=$twostem2loopnt";
                                                if ($twostem2present==2 && $klammeraufmax>=4 && $twostem2loopnt<=13){# && $twostem2intermediate <=6) {
                                                    $alternativefoldpass2=1;
                                                }
                                                #print "\nTest auf 2 Stems war positiv, falls eine 1 folgt : $alternativefoldpass2";	      					   
                                                ################### Ende dieser Prfung #######################
                                                $stem2string=$string; $stem2structure=$structure;
                                                
                                                #print "\n***\n";
                   ################# Ab hier Suchroutine weg !!! ##################################################
                                                $sumenergy=0;						  
                                                if ($fold2pass==1 || $alternativefoldpass2==1){
						#print "\nFoldpassed!!!!!!!<<<<<>>>>><<<<<>>><<<<<<<<<>>>>>>>"; 
						
						##########################################################################
						
							################ Schauen ob es wirklich nt-Paare gibt ###############
							$point12=$seqscan7+2; #sollte das Ende der smsite sein!
							$point13=$seqscan7+42;
							$stem3test=substr($gensq,($point12),($point13-$point12+1));
							$string=$stem3test;
							$string = uc($string);
							$structure=$string;
							($structure, $min_en) = RNA::fold($string);							
							#print"\nStem3: $string\n       $structure";
							@structstem3=split ("",$structure);
							$klammer3=0;
							foreach $test (@structstem3){
							    if ($test eq '(') {
								$klammer3++;
							    }
							}
                            

							#print "\nKlammer3:$klammer3";
							####################################################################

							    #Abfrage ob Treffer besser ist, als der vorher, wenn nicht, keine Ausgabe !!!!!
							   
								if ($klammer3>1 && $newhit==0) {
								    $newhit=1;
                                    

								    #$stem3test=substr($gensq,($point12),($point13-$point12+1));
								    #$string=$stem3test;
								    #$string = uc($string);
								    #$length = length($string);
								    #printf("length = %d\n", $length); #if ($istty);
								    #$structure=$string;
								    #$min_en = RNA::fold($string, $structure);
								    
								    $leaderseq='';
								    




								    push (@returnvalues,$seqscan1,$stem1string,$stem1structure,$stem2string,$stem2structure,$smsite,$string,$structure); #ev noch leader !!;
								    if ($leaderseq eq '') {
									    push (@returnvalues,0);
								    }
								    else {
									push (@returnvalues,$leaderseq);
								    }
								    push (@returnvalues,$seqscan6);#angefuegt um auch die Pos der smsite zu uebermitteln, notwendig fuer den graphischen output
								    $hits++;
								    #$gensq=$gensqstore;
								    #$OldGremark=$Gremark;
								    #$Oldmotifsum=$motifsum;
								}
							    
							
						        
                                                }
                                            }                                  
                                        }                                                     
                                    }                                                                      
                                }                                                                                             
                            }                                                                          
                        }                                                                                                                
                    }                                                                                               
                }                                                                                                      
            }                                                                                                
        }                                                                                                                              
    }
    push (@returnvalues,$hits);
    return @returnvalues;
}


sub ciona {#seems to work for the moment !!!!!!

    #This subroutine takes a sequence as first parameter and the consensusstrength as second parameter.
    #It will return an array containing the results with 12 strings:
    # $arr[0]=smsite at Position: $arr[1],$arr[2] distance smsite ggua,$arr[3]=gguposition,$arr[4] is the mfe for str, $arr[5]=firstclosed,$arr[6]=lastclosed,$arr[7]=lengthofstemclose
    #,$arr[8]=closedonesofstem1,$arr[9]=stem1loopstart,$arr[10]=firstopen,$arr[11]=openonesofstem1,$arr[12]=string,$arr[13]=structure);
    # the last file of array contains 0 for nothing found, 1 for hit, 2 for ggua+smsite but no hit


    $RNA::noLonelyPairs=1;	#important for correct use of the RNA-module!
	    #This subroutine will take the parameters in the following form: ($gensq, $consensusstrength) and it will return an array of results !
    my @returnvalues=();
    my $gensq=$_[0];
    my $consensusstrength=$_[1];
    $gensq=Trans::common::correctgensq($gensq);
    #print "\nCorrected gensq: $gensq\n";
    my $hit=0; #knows if it is a hit and it will be returned! For debugging reasons only
    while ($gensq=~m/((ggua[aguc]{2,50}?)([a|g|c][a|g]u*[^u]?uuu+[a|g][a|g]))/g){ #non-greedy version #Beautiful line and really works well !!! 
	#print "SITE!\n";
	my $len=length ($2);	#$1 contains seq with smsite, $2 contains the seq without the smsite !!! Now checking $2
	my $poswithsmsite=pos($gensq);	####ACHTUNG !! Hier liegt eine Stolperfalle !!!! obwohl wir jetzt nur $2 auslesen, steht pos !! immer noch auf dem Ende der smsite !!!!!
	my $poswithoutsmsite=$poswithsmsite-length ($3);
	my $gguaposition=$poswithsmsite-length($1);
	#print "\nggua-Position=$gguaposition\n";
	#Please comment-out one of the next two lines, depending if you want with or without smsite !!
	#$pos=$poswithsmsite;
	my $pos=$poswithoutsmsite;
    
	#Having extracted this sequence, we should now extract the part to be folded !!!
	my $unten1=35;	#It seems as if 35 is the best, considering the amount of nts in the examples !!!
	$unten1=$pos-$len if ($pos-$len<35); #now we have extracted the part of the ggua till end of smsite and the upstreamnts as in $unten1
	my $checkseq=substr($gensq,$pos-$len-$unten1,$unten1+$len);
	my $string=$checkseq;
	my $structure=$string;   #obligatory! weird way to allocate space, but absolutely needed !!!!
	my $min_en;
    ($structure, $min_en) = RNA::fold($string);
	my @struct=split ("",$structure);
	#print "folded seq: $string\nfolded str: $structure\n";
	my @stringarr=split ("",$string);
	my $ggua=$pos-$len+1; # This shall now point on the ggua, but I'm not sure, if it really points directly to the ggua
	#
	my$gguafold=$unten1; #This shall point to the ggua in the folded structure and it really does? But why? Doesn't matter !!
       
	my $letztezu=0;
	my $count1=$gguafold;
	#One thing all the testseq have is a ggua at the splice site, followed by an ag, this might be by chance, but never mind, we'll use it
	my $gguaag=0;
	$gguaag=1 if ($stringarr[$gguafold+4] eq 'a' && $stringarr[$gguafold+5] eq 'g');
	my $stem1klammerzu=0;		# Bitte noch die Anzahl der Klammern ermitteln und vielleicht noch die Energie ermitteln
	while (defined $struct[$count1] && $struct[$count1]=~/\)|\./) {
	    $letztezu=$count1 if ($struct[$count1] eq ')');
	    $stem1klammerzu++ if ($struct[$count1] eq ')');
	    $count1++;				#Diese Schleife ermittelt jetzt das Ende des Stems vom ggua aus gesehen ! Die naechste Schleife sollte den Start finden !!
	}
	#
	$count1=$gguafold-1;
	my $erstezu=0;
	while (defined $struct[$count1] && $struct[$count1]=~/\)|\./) {
	    $erstezu=$count1 if ($struct[$count1] eq ')');
	    $stem1klammerzu++ if ($struct[$count1] eq ')');
	    $count1--;				#Diese Schleife ermittelt jetzt das Ende des Stems vom ggua aus gesehen ! Die naechste Schleife sollte den Start finden !!
	}
	#
	# Laenge von stem 1 nach vorne ermitteln
	$count1=$erstezu-1;
	my $stem1loop=0;
	while (defined $struct[$count1] && $struct[$count1] eq '.') {
	    $stem1loop=$count1;
	    $count1--;
	} # Now we have the $stem1loop pointing to the most upstream nt of loop1, so we can calculate the number of nts of the loop via difference of $erstezu stem1loop
	######
	my $stem1klammernauf=0;
	my $ersteauf=0;
	while (defined $struct[$count1] && $struct[$count1]=~/\(|\./) {
	    $ersteauf=$count1 if ($struct[$count1] eq '(');
	    $stem1klammernauf++ if ($struct[$count1] eq '(');
	    $count1--;				#Diese Schleife ermittelt jetzt das Ende des Stems vom ggua aus gesehen ! Die naechste Schleife sollte den Start finden !!
	}
	######
	my $diff=$letztezu-$erstezu+1;
	$diff='0' if ($letztezu==0 || $erstezu==0);
	my $stem1loopabs=$gguaposition-$unten1+$stem1loop;
	#Attention, now checking knock out criteria for output !
	
	if ($erstezu>=15 && $diff>=7 && $stem1klammerzu>=7 && $stem1klammerzu==$stem1klammernauf && $gguaag==1 && $min_en<=-14) {
	    $hit++;	
	    #print "smsite:$3 at pos: $pos len: $len\n";
	    #print "Energie: $min_en Erstezu: $erstezu Letztezu: $letztezu LaengeStemzu: $diff stem1klammerzu: $stem1klammerzu\n";
	    #print "Stem1loopStart= $stem1loop ErsteKlammerAuf $ersteauf Stem1KlammerAufGesamt $stem1klammernauf\n";
	    #print "         1         10        20        30        40        50        60        70        80\n";
	    #print "SeqExtr: $string\n";
	    #print "SeqFold: $structure\n";
	    #Attention for the output: ALL VARIBALES THAT ARE DIRECTLY POINTING TO A NT HAVE AN INDEX THAT IS 1 TO SMALL, DUE
	    #TO THE ARRAY STARTING WITH FIELD 0 !!!!!!!!! YOU SHALL CORRECT THIS HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    #TAKE CARE OF THESE FOR THE FINAL PROGRAM OUTPUT !!!!!!!!!!!!
	    #ATTENTION: WE STILL HAVE A MIXTURE OF VARIABLES POINTING TO THE ABSOLUTE NTS AND SOME ARE POINTING TO THE RELATIVE NTS !!!
	    #YOU OUGHT TO CHANGE THIS !!!!!
	    #push (@returnvalues,$3,$pos,$len,$gguaposition,$min_en,$erstezu,$letztezu,$diff,$stem1klammerzu,$stem1loopabs,$ersteauf,$stem1klammernauf,$string,$structure);
	    push (@returnvalues, $gguaposition+1,$pos+1,$3,$min_en,$string,$structure);
	    #gguapos, with smsite at $pos with $seq($3)

	    #print "I'd locate to start of the loop to: gguaposition-unten1+stem1loop(ev+1): $stem1loopabs"; #shall be okay !!		
	    #print "xxxxxxxxxx\n@returnvalues";
	}
    }
    push (@returnvalues,$hit); #here we write the number of hits so we know how often we will print out the results!
    return @returnvalues;
}

return 1;
