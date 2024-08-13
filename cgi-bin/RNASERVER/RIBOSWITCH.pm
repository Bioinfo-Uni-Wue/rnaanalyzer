package RNASERVER::RIBOSWITCH;

#Riboswitch.pm

sub findRiboswitch {
    
    $pathtornasuboptandrnafold='/mnt/c/Users/ama55id/Nextcloud/RNA_analyzer/rnaanalyzer/bin/ViennaRNA-2.6.4/src/bin/RNAsubopt -s -e 0.75'; #Please locate your version of RNAfold / RNAsubopt
    $searchsq=$_[0]; #import the Gensq
    $direction=$_[1];
    $consensustype=-1;
    $consensustype=$_[2]; #1=strict 2=general 3=loose
    $consensustype=1 if ($consensustype==-1);
    @returnarray=();
    $overallhits=0;
    #this monster regular expression shall contain the following elements: $1 -> complete seq, $2->Stem P1 5' Part $3->P2 5'Part $4->P2 3'Part $5->P3 5'Part $6->P3 3'Part $7->P3-P1 juction $8->P1 3'Part
    $strictconsensus="([ACGU]{2}([ACGU]{9}U[AU])UA([ACGU]{5,8})[ACGU]AU[ACGU]{1,2}GG([ACGU]{6,8})GU[ACGU][UC]CUAC([ACGU]{5,7})CC[ACGU]{2,3}AA([ACGU]{5,7})(GC|AC|AGC|AU|UC)([ACGU]{12}))";
     $generalconsensus="([ACGU]{2}([ACGU]{9}U[AU])UA([ACGU]{5,8})[ACGU][ACGU]U[ACGU]{1,2}GG([ACGU]{6,8})GU[ACGU][ACGU]CUAC([ACGU]{5,7})CC[ACGU]{2,3}[ACGU][ACGU]([ACGU]{5,7})(GC|AC|AGC|AU|UC)([ACGU]{12}))";
     $looseconsensus="([ACGU]{2}([ACGU]{9}U[AU])UA([ACGU]{5,8})[ACGU][ACGU][ACGU][ACGU]{1,2}[ACGU][ACGU]([ACGU]{6,8})GU[ACGU][ACGU]CUAC([ACGU]{5,7})[ACGU][ACGU][ACGU]{2,3}[ACGU][ACGU]([ACGU]{5,7})(GC|AC|AGC|AU|UC)([ACGU]{12}))";
    $consensusweusenow=$strictconsensus if ($consensustype==1);
    $consensusweusenow=$generalconsensus if ($consensustype==2);
    $consensusweusenow=$looseconsensus if ($consensustype==3);

    
    while ($searchsq=~/$consensusweusenow/g) {
        $overallhits++;
        $poshit=pos ($searchsq)-length($1)+1;
        $poshit=(length $searchsq)-$poshit+1 if ($direction eq 'minus'); #corrects the invertation of the strand
        #print "\n----------------------\n$id$DEcomplete\n$_[1] strand at position $poshit: \n$1\n";
        $hitseq=$1;
        @rnafoldingstruct=();
        @rnafoldingenergy=();
            
        @rnafoldinganswer=`echo $1 | $pathtornasuboptandrnafold`;
        #directs to the last hit in @rnafoldinganswer
        for ($rnafoldingstructcount=1;$rnafoldingstructcount<=@rnafoldinganswer-1;$rnafoldingstructcount++){	    
            $rnafoldinganswer[$rnafoldingstructcount]=~/([.()]+)[ ]+\(?([-0-9.]+)\)?/;
            push @rnafoldingstruct,$1;
            push @rnafoldingenergy,$2;
        }
        
        $suboptstructures=@rnafoldingstruct; #now we know how many subopt structures we will return for that special hit
        
        push @returnarray, $suboptstructures, $direction, $poshit, $hitseq; #the first part of this hit is saved
        
        #@rnafoldingenergymax=sort {$a<=>$b} @rnafoldingenergy;
        #$rnafoldinganswerlines=@rnafoldingstruct-1;
        #print @rnafoldingstruct;
        $csvout=1;
        for ($temp=0;$temp<=@rnafoldingstruct-1;$temp++){
            #print "$rnafoldingstruct[$temp] $rnafoldingenergy[$temp]\n";
            push @returnarray, $rnafoldingstruct[$temp],$rnafoldingenergy[$temp]; #2nd part of hit is saved
            &checkFolding($rnafoldingstruct[$temp]);
        }
        #&checkFolding($rnafoldingstruct[0]);
    }
    push @returnarray,$overallhits;
    return @returnarray;
    
}

sub checkFolding {
    #this routine shall check a certain folding with the question, if it
    #consistent with the rules for a riboswitch
    $folding=$_[0];
    #$folding='..(((.(((((((....((((.........))))..........((((((.......))))))..))))))).))).';
    @fold=split ('',$folding);
    $foldLen=@fold-1;
    #print "\n@fold\n";
    
    $p1=-1;$p2=-1;$p3=-1;$p4=-1;$p5=-1;$p6=-1;$p7=-1;$p8=-1;$p9=-1;$p10=-1;$p11=-1;$p12=-1;
    $kl12=0;$kl56=0; $kl78=0;
    
    #      stem1        stem2                       stem3
    # ..(((((((((((....((((((.......)))))).........(((((((.......)))))))..)))))))).)))
    # p 1         2    3    4       5    6         7     8       9     10 11         12
    
    for ($var1=0;$var1<=$foldLen;$var1++) {
        if ($fold[$var1] eq '(') {    #okay, scannen die seq und versuchen p1, p4, p5, p6, p7, p8, p9, p12 zuzuordnen
            $p1=$var1 if ($p1==-1); #wenn p1 noch nicht zugeordnet ist, dann muss es dies sein.
            $p4=$var1 if ($p1!=-1 && $p5==-1);
            $p7=$var1 if ($p7==-1 && $p5!=-1);
            $p8=$var1 if ($p7!=-1 && $p9==-1);
        }
        
        if ($fold[$var1] eq ')') {
            $p5=$var1 if ($p5==-1);
            $p6=$var1 if ($p5!=-1 &&$p7==-1);
            $p9=$var1 if ($p9==-1 && $p7!=-1);
            $p12=$var1;
        } 
    }
    
    #versuch, p2 und p3 festzulegen: zuerst klammern zw p5 und p6 ermitteln und dann soviele klammern von p4 zurueck zaehlen
    
    for ($var1=$p5;$var1<=$p6;$var1++){ #zaehlt die Klammern
        $kl56++ if ($fold[$var1] eq ')');
    }
    $klzaehler=0;
    for ($var1=$p4;$var1>=0;$var1--) { #ermittelt p3
        $p3=$var1;
        $klzaehler++ if ($fold[$var1] eq '(');
        last if ($klzaehler==$kl56);
    }
    for ($var1=$p3-1;$var1>=0;$var1--){ #ermittelt p2
        $p2=$var1;
        last if ($fold[$var1] eq '(');
    }    
    
    #versuch, diesesmal nun p10 und p11 festzulegen
    
    for ($var1=$p7;$var1<=$p8;$var1++){ #ermittelt die Klammern zw p7 und p8
        $kl78++ if ($fold[$var1] eq '(');
    }
    $klzaehler=0;
    for ($var1=$p9;$var1<=$foldLen;$var1++) { #ermittelt p10
        $p10=$var1;
        $klzaehler++ if ($fold[$var1] eq ')');
        last if ($klzaehler==$kl78);
    }
    for ($var1=$p10+1;$var1<=$foldLen;$var1++) {
        $p11=$var1;
        last if ($fold[$var1] eq ')');
    }
    
    # Jetzt fehlt nur noch die Anzahl der Klammer zw p1 und p2 = kl12
    for ($var1=$p1;$var1<=$p2;$var1++) {
        $kl12++ if ($fold[$var1] eq '(');
    }
    
    #hitquality
    $hitquality=0;
    if ($kl12>=5 && $kl56>=5 && $kl78>=5 && $p2>=5) { #means this is a good hit
        $hitquality='good  ';
    }
    elsif ($kl12>=3 && $kl56>=3 && $kl78>=3 && $p2>=4) {#middle / bad
        $hitquality='middle';
    }
    else {
        $hitquality='bad   ';
    }
    
    #print OUT "$folding\nStemKlammern: $kl12  $kl56  $kl78\n";
    #print OUT "Erg:\n$p1  $p2  $p3  $p4  $p5\n$p6  $p7  $p8  $p9  $p10\n$p11  $p12\n";
    #print  "KlammernStems: $kl12 $kl56 $kl78 Positionen: $p1 $p2 $p3 $p4 $p5 $p6 $p7 $p8 $p9 $p10 $p11 $p12\n";
    #print  "$acc;$poshit;$direction;$DEshort;$hitseq\n" if ($csvout==1);
    $csvout=0;

    push @returnarray, $hitquality, $kl12, $kl56, $kl78, $p1, $p2, $p3, $p4, $p5, $p6, $p7, $p8,$p9,$p10,$p11,$p12;
}
1
