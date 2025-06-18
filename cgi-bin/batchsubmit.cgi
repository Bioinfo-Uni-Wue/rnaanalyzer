#!/usr/bin/perl

use lib "."; 
use lib "./RNASERVER/";
use strict; 
use warnings;
use CGI;
use File::Slurp;
use Bio::SeqIO;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use RNASERVER::JobUtil qw(get_next_job_id);
use JSON;
use CGI qw(:standard escapeHTML);


$CGI::POST_MAX = 1 * 1024 * 1024;  # 1 MB

my $q;

eval {
    $q = CGI->new;
};
if ($@ || !$q) {
    print CGI::header(), "<p style='color:red;'>Upload failed: file size exceeds 2 MB limit.</p>";
    exit;
}

$| = 1;  # Autoflush stdout

print $q->header();

print <<'HTML';
<html><head>
  <title>Batch Results</title>
  <link rel="stylesheet" href="/css/batchresults.css">
</head><body>
<header>
    <a href="http://localhost">    <!--- change after putting on server -->
      <img src="../images/logo.png" alt="RNA Analyzer Logo" class="logo" />
    </a>
    <div class="header-text">
      <h1>RNA Analyzer 2.0</h1>
      <p>Webserver for RNA Sequence Overview</p>
    </div>
    <div class="header-links">
      <a href="../htdocs/about.html" target="_blank">About</a> |
      <a href="../htdocs/contact.html" target="_blank">Contact</a> |
      <a href="https://www.biozentrum.uni-wuerzburg.de/bioinfo" target="_blank">Dandekar Lab</a>
    </div>
  </header>

<div class="results-container">
<h2>Results</h2>
HTML

# Read form inputs
my $input     = $q->param("SEQUENCE") // '';
my $upload_fh = $q->upload("fasta_file");
my $ire       = $q->param("IRE") // '';
my $trans     = $q->param("TRANS") // '';
my $rnamo     = $q->param("RNAmotif") // '';
my $mirna     = $q->param("mirna") // '';
my $trna      = $q->param("trna") // '';
my $coding    = $q->param("run_coding") // '';
my $species   = $q->param("species") // '';
my $mirna_target    =  $q->param("mirna_target") // '';

# Detect refresh based on hidden job_ids
my @job_ids = $q->multi_param('job_ids');
my $is_refresh = scalar @job_ids > 0;

my @fasta_blocks;

if ($input =~ /\S/) {
    @fasta_blocks = read_pasted_fasta($input);
}
elsif (defined $upload_fh) {
    @fasta_blocks = read_uploaded_fasta($upload_fh);
}

if (scalar(@fasta_blocks) > 5) {
    print "<p style='color:orange;'>Error: Not more than 5 sequences are accepted at a given time.</p>";
    exit;
}

sub read_uploaded_fasta {
    my ($fh) = @_;
    die "Undefined filehandle" unless defined $fh;

    my @blocks;
    my $seqio = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');

    while (my $seq = $seqio->next_seq) {
        my $header = $seq->id . ' ' . ($seq->desc // '');
        my $seqstr = $seq->seq;
        push @blocks, ">$header\n$seqstr";
    }
    return @blocks;
}


sub read_pasted_fasta {
    my ($text) = @_;
    $text =~ s/\r//g;

    # Add dummy header if missing
    unless ($text =~ /^>/m) {
        $text = ">Your sequence\n$text";
    }

    my @blocks = split(/^>/m, $text);
    shift @blocks if $blocks[0] !~ /\S/;

    return map { ">" . $_ } @blocks;
}

my $has_valid_input = 0;

# On first submission, launch jobs
unless ($is_refresh) {

    foreach my $block (@fasta_blocks) {
        my $sanitized = sanitize_sequence($block);


        unless ($sanitized->{valid}) {
            next;  # Skip invalid entries
        }

        $has_valid_input = 1;

        my $job_id = get_next_job_id();
        my $job_dir = "../tmp/jobs/job_$job_id";
        mkdir $job_dir;

        # write_file("$job_dir/input.txt", ">$sanitized->{name}\n$sanitized->{cleaned_seq}\n");

        my %params = (
            job_id          => $job_id,
            IRE             => $ire,
            TRANS           => $trans,
            RNAmotif        => $rnamo,
            mirna           => $mirna,
            trna            => $trna,
            run_coding      => $coding,
            species         => $species,
            mirna_target    => $mirna_target,
            dnarna          => $sanitized->{type},
            sequence_name   => $sanitized->{name},
            sequence_clean  => $sanitized->{cleaned_seq},
            sequence_length => $sanitized->{length},
            passsequence    => $sanitized->{valid},
        );

        write_file("$job_dir/params.json", encode_json(\%params));
        system("perl ./webserver_AA.cgi $job_id > /dev/null 2>&1 &");

      push @job_ids, $job_id;
  }
    unless ($has_valid_input) {
        print "<p style='color:red;'>Error: No valid RNA/DNA sequence found in your input.</p>";
        exit;
    }
}


sub sanitize_sequence {
    my ($fasta_block) = @_;
    $fasta_block =~ s/\r//g;  # Remove carriage returns from windows or older systems 
    if (length($fasta_block) > 20_000) {
      die "Input too long.";
    }

    # Extract header and sequence
    my ($header, $seq_body);

    # Clean up any extra whitespace
    $fasta_block =~ s/^\s+//;
    $fasta_block =~ s/\s+$//;

    if ($fasta_block =~ /^>([^\n]*)\n([\s\S]*)/) {
        $header   = $1 || "Your sequence";
        $seq_body = $2 || "";
    } else {
        die "$fasta_block";
        $header   = "Your sequence";
        $seq_body = $fasta_block;
    }

    # SECURITY: sanitize header to remove any shell-injection risk or path traversal
    $header =~ s/[^ a-zA-Z0-9><;+\-_.]/_/g;   # allow only safe characters
    $header =~ s/\.\.//g;                    # prevent "../" path tricks
    $header = substr($header, 0, 20);       # truncate excessively long names

    # Clean input sequence
    $seq_body = lc($seq_body);                   # Lowercase everything
    $seq_body =~ s/[^a-z]//g;                    # Remove non-letters for this analysis

    # Count base characters before sanitizing
    my $letter_count = length($seq_body);
    my $valid_base_count = ($seq_body =~ tr/atgcu//);

    my $valid = 0;
    if ($letter_count >= 10) {
        my $valid_percent = $valid_base_count / $letter_count;
        $valid = 1 if $valid_percent >= 0.9;
    }

    # Only continue if it might be valid
    # Clean and standardize sequence: replace invalids with 'n', convert t to u
    $seq_body =~ s/[^acgtu]/n/g;
    $seq_body =~ tr/t/u/;

    # Final cleaned sequence
    my ($cleaned_seq) = ($seq_body =~ /([acgun]+)/);
    $cleaned_seq ||= '';
    my $length = length($cleaned_seq);

    # Determine RNA or DNA
    my $ucount = ($cleaned_seq =~ tr/u//);
    my $tcount = ($cleaned_seq =~ tr/t//);  # Should be 0, but included for safety
    my $type = 'unknown';
    if ($ucount > $tcount) {
        $type = 'RNA';
    } elsif ($tcount > $ucount) {
        $type = 'DNA';
    }
    $type = 'unknown' if ($ucount > 20 && $tcount > 20);

    return {
        name        => $header,
        cleaned_seq => $cleaned_seq,
        type        => $type,
        length      => length($cleaned_seq),
        valid       => $valid,
    };
}



# Check for completion
my $all_done = 1;
foreach my $jid (@job_ids) {
    unless (-e "../tmp/jobs/job_$jid/result.txt") {
        $all_done = 0;
        last;
    }
}



# Table display
print "<table>";
my $elapsed = $q->param('elapsed_time') || 0;
print "<p>Elapsed time: <span id='timer'>$elapsed</span> seconds</p>" unless $all_done;

print "<tr><th>#</th><th>Job ID</th><th>Sequence Name</th><th>Status</th></tr>";

for (my $i = 0; $i < @job_ids; $i++) {
    my $jid = $job_ids[$i];
    my $result_file = "../tmp/jobs/job_$jid/result.txt";
    my $param_file  = "../tmp/jobs/job_$jid/params.json";

    my $seq_name = "Unknown";
    if (-e $param_file) {
        my $param_json = read_file($param_file);
        my $param_data = eval { decode_json($param_json) };
        $seq_name = $param_data->{sequence_name} if $param_data && $param_data->{sequence_name};
    }
   
    my $status = -e $result_file
        ? "<a class='status-button' href='../tmp/jobs/job_$jid/result.html' target='_blank'>View Result</a>"
        : "Processing";

    print "<tr><td>" . ($i + 1) . "</td><td>$jid</td><td>" . CGI::escapeHTML($seq_name) . "</td><td>$status</td></tr>";
}
print "</table>";

# Completion or waiting message
if ($all_done) {
    print "<p class='done-msg'>All jobs completed.</p>";
    print "<p>Elapsed time: <span id='timer'>$elapsed</span> seconds</p>";
} else {
    print "<p class='waiting-msg'> Waiting for jobs to finish...</p>";
}

# Hidden form to keep job IDs across refreshes
print "<form id='refreshForm' method='POST' action='batchsubmit.cgi'>\n";
foreach my $jid (@job_ids) {
    print "<input type='hidden' name='job_ids' value='$jid'>\n";
}
print "<input type='hidden' id='elapsed_time' name='elapsed_time' value='$elapsed'>\n";
print "</form>\n";

# Auto-refresh script if still processing
unless ($all_done) {
    print <<'JS';
  <script>
  let seconds = parseInt(document.getElementById('timer').textContent);
  const timerEl = document.getElementById('timer');
  const hiddenInput = document.getElementById('elapsed_time');

  function updateTimer() {
    seconds++;
    timerEl.textContent = seconds;
    hiddenInput.value = seconds;
  }

  setInterval(updateTimer, 1000);

  // Auto-refresh
  setTimeout(function() {
    document.getElementById('refreshForm').submit();
  }, 5000);
  </script>
JS
}



print "</div></body></html>";
