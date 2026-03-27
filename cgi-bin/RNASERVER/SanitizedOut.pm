package RNASERVER::SanitizedOut;

use strict;
use warnings;
use Exporter qw(import);
use CGI qw(escapeHTML);
use File::Path qw(make_path);

our @EXPORT_OK = qw(sanitized_output_link create_sanitized_copy);

sub create_sanitized_copy {
    my (%args) = @_;

    my $raw_file = $args{raw_file} or die "Missing raw_file";
    my $out_dir  = $args{out_dir}  or die "Missing out_dir";
    my $out_name = $args{out_name} or die "Missing out_name";

    make_path($out_dir) unless -d $out_dir;

    my $safe_file = "$out_dir/$out_name";

    open my $in,  '<', $raw_file  or die "Cannot open raw file '$raw_file': $!";
    open my $out, '>', $safe_file or die "Cannot write safe file '$safe_file': $!";

    while (my $line = <$in>) {
        next if $line =~ /^\s*(?:Running|Command|Executing|Exec|CMD)\s*[:=]/i;

        next if $line =~ m{
            (?<![A-Za-z0-9._-])
            /
            (?:[A-Za-z0-9._-]+/)+
            [A-Za-z0-9._-]*
        }x;

        $line =~ s/\b(api[_-]?key|token|password|passwd|rnaanalyzer)\b\s*[:=]\s*\S+/$1=[REDACTED]/ig;

        print {$out} $line;
    }

    close $in;
    close $out;

    return $safe_file;
}

sub sanitized_output_link {
    my (%args) = @_;

    my $raw_file = $args{raw_file} or return '';
    my $out_dir  = $args{out_dir}  or return '';
    my $out_name = $args{out_name} or return '';
    my $out_path = $args{out_path} or return '';
    my $label    = $args{label} || 'View Raw Output';

    create_sanitized_copy(
        raw_file => $raw_file,
        out_dir  => $out_dir,
        out_name => $out_name,
    );

    $label = escapeHTML($label);
    my $link = "$out_path/$out_name";

    return qq{<a class="sanitized-out-btn" href="$link" target="_blank" rel="noopener noreferrer">$label</a>};
}

1;