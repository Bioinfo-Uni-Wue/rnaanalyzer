package RNASERVER::IRE;

use strict;
use lib "/rnaanalyzer/bin/ViennaRNA-2.7.0/interfaces/Perl";
use lib "/usr/lib/perl5/5.26.1";
use lib "/rnaanalyzer/cgi-bin/";
use RNA;
use warnings;

$RNA::noLonelyPairs=1;	#important for correct use of the RNA-module!

# uses RNAsubopt from viennaRNA 
# $pathtornasuboptandrnafold='ViennaRNA-2.7.0/bin/RNAsubopt'; #Please locate your version of RNAfold / RNAsubopt
# added more comments for better understanding for future updates


sub suboptimalfindire {
    my $sequence = $_[0];

    # Configuration parameters
    my $loop_downstream = 17;
    my $loop_upstream   = 22;

    # quality assessment
    my $upper_stem_cutoff_good = 4;
    my $lower_stem_cutoff_good = 5;
    my $upper_stem_cutoff_bad  = 3;
    my $lower_stem_cutoff_bad  = 3;
    my $hit_before_cds_cutoff  = 200;
    my $hit_after_cds_cutoff   = 2000;

    # Path to RNAsubopt executable
    my $rnasubopt_path = '/rnaanalyzer/bin/ViennaRNA-2.7.0/bin/RNAsubopt';

    my @subopt_hits = ();

    # prevent boundary errors
    $sequence = "nnnnnnnnnn" . lc($sequence) . "nnnnnnnnnn";

   
    my $evaluate_candidate = sub {
    my ($loop_position, $require_bulged_c, $upper_good, $lower_good, $upper_bad, $lower_bad, $hit_label) = @_;

        my $loop_present       = 0;
        my $bulged_c_present   = 0;
        my $upper_stem_paired  = 0;
        my $lower_stem_paired  = 0;
        my $pass_hit           = 0;
        my @folding_energies   = ();
        my @folding_structures = ();

        # sequence for folding analysis
        my $fold_string = substr(
            $sequence,
            $loop_position - $loop_downstream,
            $loop_downstream + $loop_upstream
        );

        # skip invalid windows
        return if !defined $fold_string || length($fold_string) < ($loop_downstream + $loop_upstream);

        # subopt to get alternate foldings
        my @rna_folding_results = `echo $fold_string | $rnasubopt_path`;

        # Parse folding results
        for my $result_line_index (1 .. @rna_folding_results - 1) {
            if ($rna_folding_results[$result_line_index] =~ /([.()]+)[ ]+\(?([-0-9.]+)\)?/) {
                push @folding_structures, $1;
                push @folding_energies,   $2;
            }
        }

        # adding a limit of max 3 alternate structures
        my $max_structures = 3;
        if (@folding_structures > $max_structures) {
            @folding_structures = @folding_structures[0 .. $max_structures - 1];
            @folding_energies   = @folding_energies[0 .. $max_structures - 1];
        }

        # Analyze each suboptimal structure
        for my $structure_index (0 .. @folding_structures - 1) {

            # Reset variables for each structure analysis
            $loop_present      = 0;
            $bulged_c_present  = 0;
            $upper_stem_paired = 0;
            $lower_stem_paired = 0;
            $pass_hit          = 0;

            my @fold_structure = split('', $folding_structures[$structure_index]);

            # Check for required loop structure: 6 unpaired nucleotides with flanking pairs
            if ($fold_structure[$loop_downstream - 2] eq '(' &&
                $fold_structure[$loop_downstream - 1] eq '.' &&
                $fold_structure[$loop_downstream]     eq '.' &&
                $fold_structure[$loop_downstream + 1] eq '.' &&
                $fold_structure[$loop_downstream + 2] eq '.' &&
                $fold_structure[$loop_downstream + 3] eq '.' &&
                $fold_structure[$loop_downstream + 4] eq '.' &&
                $fold_structure[$loop_downstream + 5] eq ')') {
                $loop_present = 1;
            }


            # Check for bulged cytosine (unpaired nucleotide in stem)
            if ($fold_structure[$loop_downstream - 7] eq '.') {
                $bulged_c_present = 1;
            }

            # Count paired nucleotides in upper stem (above bulged C)
            for my $upper_stem_position (1 .. 5) {
                if ($fold_structure[$loop_downstream - 7 + $upper_stem_position] eq '(') {
                    $upper_stem_paired++;
                }
            }

            # Count paired nucleotides in lower stem (below bulged C)
            for my $lower_stem_position (1 .. 11) {
                if ($fold_structure[$loop_downstream - 7 - $lower_stem_position] eq '(') {
                    $lower_stem_paired++;
                }
            }

            # Evaluate hit quality based on structural requirements
            if ($loop_present == 1 && (!$require_bulged_c || $bulged_c_present == 1)) {
                if ($upper_stem_paired >= $upper_good &&
                    $lower_stem_paired >= $lower_good) {
                    $pass_hit = 1; # Good hit
                }
                elsif ($upper_stem_paired >= $upper_bad &&
                       $lower_stem_paired >= $lower_bad) {
                    $pass_hit = 2; # not a good hit
                }
            }

            # Store hit data if it passes quality filters
            if ($pass_hit > 0) {
                push @subopt_hits,
                    $loop_position - 10,                 # Loop position (corrected for padding)
                    $pass_hit,                           # Quality (1=good, 2=bad)
                    $fold_string,                        # Folded sequence
                    $folding_structures[$structure_index], # Secondary structure
                    $folding_energies[$structure_index],   # Folding energy
                    $upper_stem_paired,                  # Upper stem paired count
                    $lower_stem_paired,                  # Lower stem paired count
                    $hit_label;
            }
        }
    };


    # canonical IRE scan
    while ($sequence=~m/c[acgu]{5}(([c|u]agu[g|a][a|c|u])|(ccgagcu)|(cugggc)|(ccgcgc)|(gcgccg)|(gagucg)|(gagagu))/g){
        my $loop_position = (pos($sequence) - length($1)) + 1;

        $evaluate_candidate->(
            $loop_position,
            1,  # require bulged C
            $upper_stem_cutoff_good,
            $lower_stem_cutoff_good,
            $upper_stem_cutoff_bad,
            $lower_stem_cutoff_bad,
            'canonical'
        );
    }

    # non canonical but verified IRE- PFN2 from Luscieti et al blood 2017
    while ($sequence =~ m/c[acgu]{5}((aaguu[acgu]))/g) {
        my $seed_start = (pos($sequence) - length($1)) + 1;

        # for non-canonical motifs, center around the seed more flexibly
        my $loop_position = $seed_start;

        $evaluate_candidate->(
            $loop_position,
            1,  # do NOT strictly require bulged C
            $upper_stem_cutoff_good,
            $lower_stem_cutoff_good,
            $upper_stem_cutoff_bad,
            $lower_stem_cutoff_bad,
            'non-canonical'
        );
    }

    return @subopt_hits;
}


# now another whole loop for EPAS1 IRE (with an u bulkfhe and a typical loop structure)
sub find_epas1_ire {
    my ($sequence) = @_;

    my $rnasubopt_path = '/rnaanalyzer/bin/ViennaRNA-2.7.0/bin/RNAsubopt';

    my $loop_downstream = 17;
    my $loop_upstream   = 22;
    my $max_structures  = 3;

    my $upper_stem_cutoff_good = 3;
    my $lower_stem_cutoff_good = 4;
    my $upper_stem_cutoff_bad  = 2;
    my $lower_stem_cutoff_bad  = 3;

    my $block_len_good_min = 7;
    my $block_len_good_max = 9;
    my $block_len_bad_min  = 6;
    my $block_len_bad_max  = 10;

    my @epas1_hits;

    # padding like original logic
    my $padded = "nnnnnnnnnn" . $sequence . "nnnnnnnnnn";

    # strict EPAS1-like local motif
    while ($padded =~ /(ugucca(cagag[acgu]))/g) {

        my $full_match    = $1;
        my $loop_seq      = $2;
        my $loop_position = (pos($padded) - length($loop_seq)) + 1;   # 1-based in padded seq
        my $reported_pos  = $loop_position - 10;                      # same correction style

        my $fold_string = substr(
            $padded,
            $loop_position - $loop_downstream,
            $loop_downstream + $loop_upstream
        );

        next if !defined $fold_string;
        next if length($fold_string) < ($loop_downstream + $loop_upstream);

        my @rna_folding_results = `echo $fold_string | $rnasubopt_path`;

        my @folding_structures;
        my @folding_energies;

        for my $i (1 .. $#rna_folding_results) {
            chomp $rna_folding_results[$i];
            if ($rna_folding_results[$i] =~ /([.()]+)\s+\(?([-0-9.]+)\)?/) {
                push @folding_structures, $1;
                push @folding_energies,   $2;
            }
        }

        next unless @folding_structures;

        if (@folding_structures > $max_structures) {
            @folding_structures = @folding_structures[0 .. $max_structures - 1];
            @folding_energies   = @folding_energies[0 .. $max_structures - 1];
        }

        for my $structure_index (0 .. $#folding_structures) {

            my $structure = $folding_structures[$structure_index];
            my $energy    = $folding_energies[$structure_index];
            my @fold_structure = split('', $structure);

            my $loop_start = $loop_downstream - 1;
            my $loop_end   = $loop_start + length($loop_seq) - 1;

            my $loop_unpaired = 0;
            for my $i ($loop_start .. $loop_end) {
                next if $i < 0 || $i > $#fold_structure;
                $loop_unpaired++ if $fold_structure[$i] eq '.';
            }

            # contiguous unpaired block containing the loop
            my $block_left  = $loop_start;
            my $block_right = $loop_end;

            $block_left--  while $block_left  >= 0                && $fold_structure[$block_left]  eq '.';
            $block_right++ while $block_right <= $#fold_structure && $fold_structure[$block_right] eq '.';

            $block_left++;
            $block_right--;

            my $unpaired_block_len = $block_right - $block_left + 1;

            my $left_block_flank_ok  = 0;
            my $right_block_flank_ok = 0;

            $left_block_flank_ok  = 1 if $block_left  - 1 >= 0                && $fold_structure[$block_left  - 1] eq '(';
            $right_block_flank_ok = 1 if $block_right + 1 <= $#fold_structure && $fold_structure[$block_right + 1] eq ')';

            my $upper_stem_paired = 0;
            my $lower_stem_paired = 0;
            my $pass_hit          = 0;

            for my $i (($block_left - 5) .. ($block_left - 1)) {
                next if $i < 0 || $i > $#fold_structure;
                $upper_stem_paired++ if $fold_structure[$i] eq '(';
            }

            for my $i (($block_right + 1) .. ($block_right + 8)) {
                next if $i < 0 || $i > $#fold_structure;
                $lower_stem_paired++ if $fold_structure[$i] eq ')';
            }

            if ($loop_unpaired == length($loop_seq) &&
                $unpaired_block_len >= $block_len_good_min &&
                $unpaired_block_len <= $block_len_good_max &&
                $left_block_flank_ok &&
                $right_block_flank_ok &&
                $upper_stem_paired >= $upper_stem_cutoff_good &&
                $lower_stem_paired >= $lower_stem_cutoff_good) {
                $pass_hit = 1;
            }
            elsif ($loop_unpaired >= length($loop_seq) - 1 &&
                   $unpaired_block_len >= $block_len_bad_min &&
                   $unpaired_block_len <= $block_len_bad_max &&
                   $left_block_flank_ok &&
                   $right_block_flank_ok &&
                   $upper_stem_paired >= $upper_stem_cutoff_bad &&
                   $lower_stem_paired >= $lower_stem_cutoff_bad) {
                $pass_hit = 2;
            }

            if ($pass_hit > 0) {
                push @epas1_hits,
                    $reported_pos,
                    $pass_hit,
                    $fold_string,
                    $structure,
                    $energy,
                    $upper_stem_paired,
                    $lower_stem_paired,
                    'canonical';
            }
        }
    }

    return @epas1_hits;
}