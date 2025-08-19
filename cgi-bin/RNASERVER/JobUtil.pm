package RNASERVER::JobUtil;
use strict;
use warnings;
use Exporter 'import';
use File::Path qw(make_path);
our @EXPORT_OK = qw(get_next_job_id);

sub get_next_job_id {
    my $dir  = "../tmp";           # Stable and writable
    my $file = "$dir/job.dat";

    # Create directory if it doesn't exist
    unless (-d $dir) {
        make_path($dir) or die "Failed to create $dir: $!";
        chmod 0777, $dir;
    }

    # Create job.dat file if it doesn't exist
    unless (-e $file) {
        open my $new, '>', $file or die "Can't create $file: $!";
        print $new "1000\n";
        close $new;
        chmod 0666, $file;
    }

    open my $fh, "+<", $file or die "Can't open $file: $!";
    flock($fh, 2) or die "Can't lock $file: $!";

    my $id = <$fh> || 1000;
    chomp $id;
    $id =~ s/\D//g;

    seek $fh, 0, 0;
    truncate $fh, 0;
    print $fh $id + 1;
    close $fh;

    return $id;
}

1;
