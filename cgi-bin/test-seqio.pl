# test-seqio.pl
#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

my $seqio = Bio::SeqIO->new(
    -file   => 'test_sequence.fasta',
    -format => 'fasta'
);

while (my $seq = $seqio->next_seq) {
    print "Sequence ID: ", $seq->id, "\n";
    print "Sequence: ", $seq->seq, "\n";
}



while (my $seq = $seqio->next_seq) {
    print "Sequence ID: ", $seq->id, "\n";
    print "Sequence: ", $seq->seq, "\n";
    print "Length: ", $seq->length, "\n";
    print "Subsequence (1-10): ", $seq->subseq(1, 10), "\n";
}

my $out = Bio::SeqIO->new(
    -file   => '>output_sequences.fasta',
    -format => 'fasta'
);

while (my $seq = $seqio->next_seq) {
    $out->write_seq($seq);
}



use Bio::DB::GenBank;

my $gb = Bio::DB::GenBank->new;
my $seq = $gb->get_Seq_by_acc('J01673');  # Example accession number

print "Fetched Sequence ID: ", $seq->display_id, "\n";
print "Description: ", $seq->desc, "\n";


eval {
    while (my $seq = $seqio->next_seq) {
        # Process each sequence
    }
};
if ($@) {
    warn "An error occurred while processing sequences: $@\n";
}
