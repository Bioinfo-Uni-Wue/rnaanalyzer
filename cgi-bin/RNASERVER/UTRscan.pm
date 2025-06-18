package RNASERVER::UTRscan;

use Cwd qw(abs_path);
use strict;
use warnings;

sub predict_utrs {
    my (%args) = @_;

    my $seq        = $args{seq};
    my $cds_start  = $args{cds_start};
    my $cds_end    = $args{cds_end};
    my $strand     = $args{strand} || '+';
    my $source     = $args{source} || 'unknown';
    my $tss_pos    = $args{tss};  # optional
    my $tts_pos    = $args{tts};  # optional

    $seq =~ tr/uU/tT/;
    my $seq_length = length($seq);

    print "[DEBUG] Source: $source\n"; 
    print "[DEBUG] CDS start: $cds_start, CDS end: $cds_end, Sequence length: $seq_length\n";
    
    unless (defined $cds_start && defined $cds_end && $cds_start < $cds_end) {
        warn "[UTRscan] Invalid or missing CDS positions: start=$cds_start, end=$cds_end\n";
        return;
    }

    our @new5primeutr = ();
    our @new3primeutr = ();
    our @utrprintout  = ();
    our @utr          = ();


    print "<i>UTR prediction source: $source</i><br>";

    if ($source eq 'augustus' && (defined $tss_pos || defined $tts_pos)) {
        print "<i>Using AUGUSTUS-predicted TSS/TTS boundaries.</i><br>";

        if (defined $tss_pos && defined $cds_start && $tss_pos < $cds_start) {
            push @new5primeutr, $tss_pos, $cds_start - 1;
            push @utrprintout, 5, $tss_pos, $cds_start - 1;
        }
        if (defined $tts_pos && defined $cds_end && $cds_end < $tts_pos) {
            push @new3primeutr, $cds_end + 1, $tts_pos;
            push @utrprintout, 3, $cds_end + 1, $tts_pos;
        }

    } else {
        print "<i>Inferring UTRs from CDS boundaries...</i><br>";

        if ($strand eq '+') {
            if (defined $cds_start && $cds_start > 1) {
                push @new5primeutr, 1, $cds_start - 1;
                push @utrprintout, 5, 1, $cds_start - 1;
            }
            if (defined $cds_end && $cds_end < $seq_length) {
                push @new3primeutr, $cds_end + 1, $seq_length;
                push @utrprintout, 3, $cds_end + 1, $seq_length;
            }
        } elsif ($strand eq '-') {
            if (defined $cds_end && $cds_end < $seq_length) {
                push @new5primeutr, $cds_end + 1, $seq_length;
                push @utrprintout, 5, $cds_end + 1, $seq_length;
            }
            if (defined $cds_start && $cds_start > 1) {
                push @new3primeutr, 1, $cds_start - 1;
                push @utrprintout, 3, 1, $cds_start - 1;
            }
        }
    }

    @utr = sort { $a <=> $b } (@new5primeutr, @new3primeutr);

    print "<b>UTR:</b>           start  -   end   -  stems - energy<br>";
    for (my $i = 0; $i < @utrprintout; $i += 3) {
        my ($type, $start, $end) = @utrprintout[$i, $i+1, $i+2];

        printf (" %d'            %-6d - %6d", $type, $start, $end);

        my $utr_seq = substr($seq, $start - 1, $end - $start + 1);
        my @returnout = &checkstemsonly($utr_seq, 1);
        print "       $returnout[0]" if ($returnout[0] == $returnout[1]);
        print "       $returnout[0]-$returnout[1]" if ($returnout[0] != $returnout[1]);
        print "     $returnout[2]<br>" if ($returnout[2] != 1);

        # Optional motif scan
        if ($type == 5) {
            my $sd = ($utr_seq =~ /AGGAGG/i) ? "SD motif" : "";
            my $kozak = ($utr_seq =~ /gcc[AG]ccATGG/i) ? "Kozak motif" : "";
            print "         $sd $kozak<br>" if $sd || $kozak;
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

    return (\@new5primeutr, \@new3primeutr, \@utrprintout, \@utr);
}

return 1;


