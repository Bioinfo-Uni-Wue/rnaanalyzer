#!/usr/bin/perl
use strict;
use warnings;
use File::Slurp;
use Cwd qw(abs_path);

# Set paths
my $rnafold_bin = '../bin/ViennaRNA-2.6.4/src/bin/RNAfold'; 
my $input_file  = './test.seq';
my $output_file = './test.foldout';

# Check everything
die "RNAfold binary not found at $rnafold_bin" unless -x $rnafold_bin;
die "Input file not found: $input_file" unless -e $input_file;

# Absolute paths
$input_file  = abs_path($input_file);
$output_file = abs_path($output_file);

# Construct command
my $cmd = "$rnafold_bin --infile=$input_file --outfile=$output_file";

print "Running RNAfold:\n$cmd\n\n";

# Run command and capture output
my $output = `$cmd 2>&1`;
my $exit_code = $? >> 8;

# Print result
print "Exit code: $exit_code\n";
print "RNAfold output:\n$output\n";

# Show folded output if successful
if ($exit_code == 0 && -e $output_file) {
    print "\n--- RNAfold .foldout contents ---\n";
    print read_file($output_file);
}

