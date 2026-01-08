package RNASERVER::IRE;

use strict;
use lib "$ENV{HOME}/rnaanalyzer/bin/ViennaRNA-2.7.0/interfaces/Perl";
use lib "/usr/lib/perl5/5.26.1";
use lib "$ENV{HOME}/rnaanalyzer/cgi-bin/";
use RNA;
use warnings;

$RNA::noLonelyPairs=1;	#important for correct use of the RNA-module!

# uses RNAsubopt from viennaRNA 
# $pathtornasuboptandrnafold='/storage/srv/bioapps/rnaanalyzer/bin/ViennaRNA-2.7.0/bin/RNAsubopt'; #Please locate your version of RNAfold / RNAsubopt
# added more comments for better understanding for future updates

sub suboptimalfindire {
    my $sequence = $_[0];
    
    # Configuration parameters
    my $loop_downstream = 17;
    my $loop_upstream = 22;
    
    #  quality assessment
    my $upper_stem_cutoff_good = 4;
    my $lower_stem_cutoff_good = 5;
    my $upper_stem_cutoff_bad = 3;
    my $lower_stem_cutoff_bad = 3;
    my $hit_before_cds_cutoff = 200;
    my $hit_after_cds_cutoff = 2000;
    
    # Path to RNAsubopt executable
    my $rnasubopt_path='/storage/srv/bioapps/rnaanalyzer/bin/ViennaRNA-2.7.0/bin/RNAsubopt';
    
    my @subopt_hits = ();
    
    #  prevent boundary errors
    $sequence = "nnnnnnnnnn" . $sequence . "nnnnnnnnnn";
    
    # cannonical sequences or variants from Henderson et al JBC 1996 + newly found from Luscieti et al blood 2017
    while ($sequence =~ m/c[acgu]{5}(([cu]agu[ga][acu])|(aaguu[acgu])|(ccgagcu)|(cugggc)|(ccgcgc)|(gcgccg)|(gagucg)|(gagugu))/g) {
        
        # Initialize variables for this hit
        my $loop_present = 0;
        my $bulged_c_present = 0;
        my $upper_stem_paired = 0;
        my $lower_stem_paired = 0;
        my $pass_hit = 0;
        my @folding_energies = ();
        my @folding_structures = ();
        
        # analyze loop position (+1 for correction)
        my $loop_position = (pos $sequence) - (length $1) + 1;
        
        #  sequence for folding analysis
        my $fold_string = substr($sequence, $loop_position - $loop_downstream, $loop_downstream + $loop_upstream);
        
        # subopt to get alternate foldings
        my @rna_folding_results = `echo $fold_string | $rnasubopt_path`;
        
        # Parse folding results
        for my $result_line_index (1 .. @rna_folding_results - 1) {
            if ($rna_folding_results[$result_line_index] =~ /([.()]+)[ ]+\(?([-0-9.]+)\)?/) {
                push @folding_structures, $1;
                push @folding_energies, $2;
            }
        }
        
	# adding a limit of max 3 alternate structures
        my $max_structures = 3;
        if (@folding_structures > $max_structures) {
            @folding_structures = @folding_structures[0..$max_structures-1];
            @folding_energies = @folding_energies[0..$max_structures-1];
        }

	
        # Analyze each suboptimal structure
        for my $structure_index (0 .. @folding_structures - 1) {
            
            # Reset variables for each structure analysis
            $loop_present = 0;
            $bulged_c_present = 0;
            $upper_stem_paired = 0;
            $lower_stem_paired = 0;
            $pass_hit = 0;
            
            my @fold_sequence = split('', $fold_string);
            my @fold_structure = split('', $folding_structures[$structure_index]);
            
            # Check for required loop structure: 6 unpaired nucleotides with flanking pairs
            if ($fold_structure[$loop_downstream - 2] eq '(' && 
                $fold_structure[$loop_downstream - 1] eq '.' && 
                $fold_structure[$loop_downstream] eq '.' &&
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
            if ($loop_present == 1 && $bulged_c_present == 1) {
                if ($upper_stem_paired >= $upper_stem_cutoff_good && 
                    $lower_stem_paired >= $lower_stem_cutoff_good) {
                    $pass_hit = 1; # Good hit
                }
                elsif ($upper_stem_paired >= $upper_stem_cutoff_bad && 
                       $lower_stem_paired >= $lower_stem_cutoff_bad) {
                    $pass_hit = 2; # not a good hit
                }
            }
            
            # Store hit data if it passes quality filters
            if ($pass_hit > 0) {
                push @subopt_hits, 
                     $loop_position - 10,           # Loop position (corrected for padding)
                     $pass_hit,                     # Quality (1=good, 2=bad)
                     $fold_string,                  # Folded sequence
                     $folding_structures[$structure_index], # Secondary structure
                     $folding_energies[$structure_index],   # Folding energy
                     $upper_stem_paired,            # Upper stem paired count
                     $lower_stem_paired;            # Lower stem paired count
            }
        }
    }
    
    return @subopt_hits;
}
