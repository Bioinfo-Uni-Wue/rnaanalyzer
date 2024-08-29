#!/usr/bin/perl


use Cwd;
use strict;
use warnings;


my $TEMPDIR='/mnt/c/Users/ama55id/Nextcloud/RNA_analyzer/rnaanalyzer/tmp';
my $VIENNARNAFOLDDIR='/mnt/c/Users/ama55id/Nextcloud/RNA_analyzer/rnaanalyzer/bin/ViennaRNA-2.6.4/src/bin'; #pointing to the RNAfold dir


sub process_rna {
    my $fold;
    my $svg;
	my $seq_file = "$TEMPDIR/job.seq";
    my $fold_file = "$TEMPDIR/job.fold";
    my $svg_file = "$TEMPDIR/job.svg";

    print "Temporary Directory: $TEMPDIR\n";
    print "Sequence File: $seq_file\n";
    print "Fold File: $fold_file\n";
    print "SVG File: $svg_file\n";

    

    # Run RNAfold with input and output specified
    system("$VIENNARNAFOLDDIR/RNAfold --infile=$seq_file | $VIENNARNAFOLDDIR/RNAplot -o svg --filename-full > $svg_file");
    
}

process_rna()