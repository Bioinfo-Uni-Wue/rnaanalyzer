package RNASERVER::IRE;

use strict;
use lib "/mnt/c/Users/ama55id/Nextcloud/RNA_analyzer/rnaanalyzer/bin/ViennaRNA-2.6.4/interfaces/Perl";
use lib "/mnt/c/Users/ama55id/Nextcloud/RNA_analyzer/rnaanalyzer/cgi-bin/";
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
