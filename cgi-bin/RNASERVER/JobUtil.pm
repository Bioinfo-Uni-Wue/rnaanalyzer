package RNASERVER::JobUtil;
use strict;
use warnings;
use Exporter 'import';
use File::Path qw(make_path);
use Digest::SHA qw(hmac_sha256);   # Core/standard module
use MIME::Base64 qw(encode_base64);# Core/standard module
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

    # secret phrase for creating the hash
    my $SECRET = 'rnaanalyzer-bioinfo';

    # HMAC-SHA256(id, secret), then URL-safe Base64, then take first 16 chars
    my $digest = hmac_sha256("$id", $SECRET);
    my $b64    = encode_base64($digest, '');   # no newline
    $b64 =~ tr|+/=|-_|d;                       # URL-safe, strip '=' padding
    my $opaque = substr($b64, 0, 6);          # length 16 (tweakable)

    return $opaque;

}

1;
