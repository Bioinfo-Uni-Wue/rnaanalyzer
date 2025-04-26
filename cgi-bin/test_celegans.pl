#!/usr/bin/perl
use strict;
use warnings;

# Adjust this path based on your project structure
use lib ".";  # <-- or adjust as needed
use RNASERVER::TRANS2;

# Optional: uncomment for better error stack tracing
# use Devel::Confess;
 
# === Input test sequence (this one gave you hits in the old webapp) ===
my $sequence = 'aacuaaaacaauuuuugaagaacaguuucuguacuucauugguauguagagacuuccagaaccuaguucugaaauuuuggaccgagcuuucgggcuuuuuuuuuacuuuuuucggguacggcgcgggugugggggggccugugggggguccuaugggggccuauggggggcccagggggggucccgggggggccggccgggcgcccaccccgcccccggcgacccgagcggcgcacgaaaaacgcgaaacgcgugcuuggcgugacucgcgcggcccgugcgcgugcguggugacaaucggcgcgugcugaaugcgacggugcgcgagauguggcggugcggagggucggcccgcgcuuuuuuugaccuccuccaaaaaagccuauaguccgcgcaagagggcggcgcuuuggggacgcgggcuggcccccgugcgccggcccccgcgcgcaggcccgcccccgacccugccccgcccugccccgcccugcccgcgccccggcccucaccccacgccacgccaccccguaucaucgauagacccccgcgcggguacauuucgcgccaaacggguuagagacucuu';
print "\n🔍 [DEBUG] Starting celegans() test...\n";
print "🔬 [DEBUG] Input sequence length: ", length($sequence), "\n";

# === Call the celegans function ===
my @results = RNASERVER::TRANS2::celegans($sequence);

# === Check how many elements returned ===
my $num_results = scalar(@results);
print "\n📦 [DEBUG] Returned $num_results values from celegans()\n";

# === Final position in the array is always hit count ===
my $hits = $results[-1];
print "🎯 [DEBUG] Number of hits: $hits\n";

# === Dump all returned data for inspection ===
if ($hits > 0) {
    print "\n📊 [DEBUG] celegans() result breakdown:\n";
    my $block_size = 10;  # celegans result structure returns 10 values per hit
    for my $i (0 .. $hits - 1) {
        my $offset = $i * $block_size;
        print "\n🧬 Hit #", $i + 1, "\n";
        print "  Position:     $results[$offset+0]\n";
        print "  Stem1:        $results[$offset+1]\n";
        print "  Structure1:   $results[$offset+2]\n";
        print "  Stem2:        $results[$offset+3]\n";
        print "  Structure2:   $results[$offset+4]\n";
        print "  Sm-Site:      $results[$offset+5]\n";
        print "  Stem3:        $results[$offset+6]\n";
        print "  Structure3:   $results[$offset+7]\n";
        print "  Leader:       ", ($results[$offset+8] eq '0' ? "none" : $results[$offset+8]), "\n";
        print "  Sm-Site Pos:  $results[$offset+9]\n";
    }
} else {
    print "\n🚫 [DEBUG] No hits were found by celegans()\n";
}

