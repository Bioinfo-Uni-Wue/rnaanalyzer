      AUGUSTUS UTR available?
           |
    Yes  ←───┴────→ No
     |                   |
Scan for polyA     Predict polyA site
in UTR region      in sequence end
     |                   |
 Add polyA info      Define 3' UTR
 to annotation       using polyA signal


##new subroutine:

To detect a polyadenylation signal or polyA tail downstream of the coding sequence (CDS) and use that to define or refine a 3′ UTR when AUGUSTUS fails to predict one.

🧠 Step-by-Step Breakdown
```perl
my $window_start = $cds_end + 1; my $window_end = $seq_length;

yaml
Copy
Edit

- The **3′ UTR should start after the coding sequence ends**.
- We are scanning the region from the **CDS end to the end of the sequence**, assuming it contains potential UTR and polyA-related elements.

---

### ```perl
my $utr_seq = substr($sequence, $window_start - 1, $window_end - $window_start + 1);
We extract the actual subsequence to scan.

-1 accounts for Perl’s 0-based indexing.

```

```perl
my @motifs = qw(AATAAA ATTAAA TATAAA AAGAAA AATATA AGTAAA);

yaml
Copy
Edit

- These are **known polyA signal motifs** (hexamers).
  - `AATAAA` is canonical (~60–70% of polyA sites).
  - The rest are common **variants** used by ~30–40% of genes.

---

### ```perl
foreach my $motif (@motifs) {
    if ($utr_seq =~ /$motif/i) {
        ...
        last;
    }
}
Scans the extracted 3′ region for any of the motifs.

The first match wins — you could make this more flexible later.

```perl
if ($sequence =~ /A{10,}/g) { ... }

yaml
Copy
Edit

- Searches for **long runs of A’s**, a hallmark of the **polyA tail** that’s added to mRNA after transcription.
- We check whether this **tail** appears **reasonably close to the CDS end** (within ~200 nt), to validate it's connected to this transcript.

---

### 🧩 Now the Key Logic: If Either Was Found

```perl
if ($signal_pos > 0 || $tail_pos > 0) {
    my $utr_start = $cds_end + 1;
    my $utr_end = $tail_pos > 0 ? $tail_pos : ($signal_pos + 20);
Start of the UTR is always after the CDS.

If a tail was found, it likely marks the true end of the transcript → use that.

If only a signal is found, add ~20 nt padding, since cleavage sites are typically ~10–30 nt downstream of the signal.

This is based on known biological structure:

🧬 polyA signals (e.g., AATAAA) appear ~10–30 nt upstream of the cleavage site (polyA tail).

💾 Then We Push to Your Annotation Arrays:
perl
Copy
Edit
push @new3primeutr, $utr_start, $utr_end;
push @utrprintout, 3, $utr_start, $utr_end;
push @utr, $utr_start, $utr_end;
This integrates the result exactly like infer_and_print_UTRs(), so the rest of your pipeline sees the same UTR structure.

🧬 Summary of Biological Logic Encoded:
Subroutine Action	Biological Meaning
Scans after CDS end	UTR is always after the coding region
Looks for polyA motifs	Signals that mRNA processing adds polyA tail
Looks for polyA tails (A{10,})	Indicates mature mRNA end
Combines signal + tail or best guess	Simulates actual cleavage/polyadenylation site
Returns 3′ UTR for annotation	Provides meaningful UTR even if AUGUSTUS fails


## how does this help with UTR prediction?

 How This Subroutine Helps with UTR Prediction
🧠 The problem:
AUGUSTUS sometimes doesn’t predict UTRs, especially 3′ UTRs.

That’s bad for biological interpretation and downstream analysis:

You can not define the true end of the transcript.

You miss important regulatory regions (e.g., microRNA binding, stability elements).

🧩 What refineUTRwithPolyA() does:
It gives you a biologically meaningful guess for the 3′ UTR even when AUGUSTUS doesn’t.

Here’s how:

🧬 1. Predicts where the transcript might end
PolyA signals like AATAAA appear near the end of mRNAs, just before the cleavage and polyA tail.

If you see such a motif (and//or a polyA tail), it’s a strong clue that the transcript likely ends nearby.

➕ Result:
You can say:

“No UTR from AUGUSTUS? No problem. Based on polyA signal/tail, UTR likely ends around position 3570.”

🛠️ 2. Defines a UTR region you can annotate
With both start (CDS end + 1) and end (polyA signal or tail), you now have:

text
Copy
Edit
3′ UTR = [cds_end + 1, polyA_site]
➕ Result:
You can now:

Print a UTR in your report

Draw it in visual sequence output

Pass it to structure folding

Add it to GFF or HTML output

🔗 3. Fits into your existing pipeline
The UTR range is stored in:

@new3primeutr

@utrprintout

@utr

These are the same arrays used by AUGUSTUS prediction logic

➕ Result:
Your app doesn’t need to care if the UTR came from AUGUSTUS or polyA logic — it just works.

######## replaced subroutine for calculating UTR
######## also call UTR from the new subroutine after augustus  

sub calcUTR {
	#Trying to calculate 3' and 5' regions! But this is really very basic done
	#This block looks for single-exon predictions, and assumes it's RNA, because multi-exon UTR logic is handled elsewhere.
	@utrprintout=();
	$numberofexons=@exons/2;
	

    if ($numberofexons==1){
	###### RNA 5' and 3' Detection ###########
	############### S T A R T ################



	###New: checking for coding sequence!!! Perhaps this should better be done in UTR, but nevermind

    ##this grabs Single-exon coding boundaries and Putative PolyA signals (used to help detect 3′ UTR end)
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
        #This defines the region from the end of the exon to just after the polyA signal, and marks it as 3′ UTR.
		#oder man sagt, die 3' region geht dann bis zum Schluss, koennte man auch sagen, das soll aber mal okay sein!
	}
	#now 5' UTR works only with the first sngl exon because here we can located the 5'UTR Start to 1 bzw 0, otherwise we can't, because we
	#do not know where ist starts, is it the PolyASignal, PolyATail or elsewhere ??
    #It assumes the 5′ UTR goes from the beginning (position 1) to just before the exon starts.
	if ($singleexonboundaries[0]==$exons[0]){ #the it might be rna and the 5'UTR is 1 to $exons[0]-1
		push @new5primeutr,1,$exons[0]-1;
	}

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