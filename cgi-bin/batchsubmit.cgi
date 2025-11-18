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
  <title>RNAanalyzer</title>
  <link rel="stylesheet" href="/css/batchresults.css">
</head><body>
<header>
    <a href="https://rnaanalyzer.bioapps.biozentrum.uni-wuerzburg.de//">    <!--- change after putting on server -->
      <img src="../images/logo.png" alt="RNA Analyzer Logo" class="logo" />
    </a>
    <div class="header-text">
      <h1>RNA Analyzer<sup>3</sup></h1>
      <p>Webserver for RNA Sequence Overview</p>
    </div>
    <div class="header-links">
      <a href="https://rnaanalyzer.bioapps.biozentrum.uni-wuerzburg.de//about.html" target="_blank">Help</a> |
      <a href="https://rnaanalyzer.bioapps.biozentrum.uni-wuerzburg.de//contact.html" target="_blank">Contact</a> |
      <a href="https://www.biozentrum.uni-wuerzburg.de/bioinfo" target="_blank">Dandekar Lab</a>
    </div>
  </header>

<div class="results-container">
<h2>JOBS</h2>
HTML

# code written by liang 2025-09-17 ---------------

use CGI qw(:standard);
use Fcntl ':flock';  # For file locking
use Time::HiRes qw(time);

my @job_ids = $q->multi_param('job_ids');
my $is_refresh = scalar @job_ids > 0;

unless ($is_refresh) {

my $max_visits = 10;          # Maximum visits allowed
my $time_frame = 60;          # Time frame in seconds (e.g., 60 seconds)
my $log_file = "../tmp/ip_log.txt";  # Log file to store IP

# Get the visitor's IP address
my $ip_address = remote_addr();
my %ip_log = ();
if (open(my $log_fh, '<', $log_file)) {
    flock($log_fh, LOCK_SH);
    while (my $line = <$log_fh>) {
        chomp $line;
        my ($ip, $last_time, $visits) = split /:/, $line;
        $ip_log{$ip} = { time => $last_time, visits => $visits };
    }
    close $log_fh;
}

# Get the current time
my $current_time = time();
# Check the current IP's visit log
my $visit_data = $ip_log{$ip_address};
if ($visit_data && $current_time - $visit_data->{time} < $time_frame) {
    # Within the time frame
    if ($visit_data->{visits} >= $max_visits) {
        # Deny access
        print header(-status => '429 Too Many Requests'),
              start_html('Access temporarily denied'),
              h1('Too Many Requests'),
              p('You have exceeded the maximum number of computing tasks. Please try again two minutes later. Thank you very much for understanding.'),
              end_html;
        exit;
    } else {
        $ip_log{$ip_address}->{visits}++;
    }
} else {
    # Reset or initialize the visit count
    $ip_log{$ip_address} = { time => $current_time, visits => 1 };
}

# Save the updated log
if (open(my $log_fh, '>', $log_file)) {
    flock($log_fh, LOCK_EX);
    foreach my $ip (keys %ip_log) {
        my $data = $ip_log{$ip};
        print $log_fh "$ip:$data->{time}:$data->{visits}\n";
    }
    close $log_fh;
}
###########################
use File::Find;
use File::Path qw(remove_tree);
use POSIX qw(strftime);
my $jobdir="../tmp/jobs";
my $time_threshold = time()-30*24*60*60; # Roughly 30 days

sub wanted {
    return unless -d $_;  # Proceed only if it's a directory
    my $dir = $File::Find::name;
    my $mtime = (stat($dir))[9];
    if (defined $mtime && $mtime < $time_threshold) {
        # print p("Removing directory: $dir");
        # Remove the directory and handle errors
        eval { remove_tree($dir) };
        if ($@) {
            print p("Failed to remove $dir: $@");
        }
    }
}
find(\&wanted, $jobdir);
#print p("completed cleanup");
###### finish cleaning up #####
#print "Computing...";
}
# code written by liang 2025-09-17 ---------------


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
my $rbp       = $q->param("RBP") // '';
my $ribo      = $q->param("RIBO") // '';


my @job_ids = $q->multi_param('job_ids');
my $is_refresh = scalar @job_ids > 0;

my @fasta_blocks;
our $formatted_time = $q->param('formatted_time') || 0;
our $count;

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


    # Count all characters in the combined text
    my $all_text = join("\n", @blocks);
    $count = length($all_text);

    return (@blocks, $count);
}


sub read_pasted_fasta {
    my ($text) = @_;
    $text =~ s/\r//g;

    $count = length($text);

    # Add dummy header if missing
    unless ($text =~ /^>/m) {
        $text = ">no_name\n$text";
    }

    my @blocks = split(/^>/m, $text);
    shift @blocks if $blocks[0] !~ /\S/;

    return map { ">" . $_ } @blocks;
}

my $has_valid_input = 0;

#  launch jobs
unless ($is_refresh) {

    my $max_length = 0;

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
            RBP             => $rbp,
	    ribo            => $ribo,
            dnarna          => $sanitized->{type},
            sequence_name   => $sanitized->{name},
            sequence_clean  => $sanitized->{cleaned_seq},
            sequence_length => $sanitized->{length},
            passsequence    => $sanitized->{valid},
        );

        write_file("$job_dir/params.json", encode_json(\%params));
        system("perl ./webserver.cgi $job_id > /dev/null 2>&1 &");

      push @job_ids, $job_id;

      $max_length = $sanitized->{length} if $sanitized->{length} > $max_length;
  }
    unless ($has_valid_input) {
        print "<p style='color:red;'>Error: Invalid characters found in the input sequence. Please check the sequence before submitting.</p>\n";
		print "<button type='button' class='back-button' onclick='history.back()'>Back to submission</button>";               # added button to go back to input
        exit;
    }

    my $estimate_per_250nt = 5;  # estiamted
    my $predicted_time = int(($count / 250) * $estimate_per_250nt);
    $predicted_time = $predicted_time / 60.0;  #for minutes

    #    #mirnatargetscan takes time so adding a min for this 
    #if ($mirna_target) {
    #    $predicted_time += 1;
    #}

    $predicted_time = 1 if $predicted_time < 1;
    if ($predicted_time == int($predicted_time)) {
        $formatted_time = int($predicted_time);
    } else {
        $formatted_time = sprintf("%.1f", $predicted_time);
    }

}


sub sanitize_sequence {
    my ($fasta_block) = @_;
    $fasta_block =~ s/\r//g;  # Remove carriage returns from windows or older systems 
    if (length($fasta_block) > 20_000) {      # how much shall be the input length?
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
        # die "$fasta_block";
        $header   = "Your sequence";
        $seq_body = $fasta_block;
    }

    # sanitize header to remove any shell-injection risk or path traversal
    $header =~ s/[^a-zA-Z0-9>+\-_.]/_/g;   # allow only safe characters
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
    # Clean and standardize sequence: replace invalids with '', convert t to u
    $seq_body =~ s/[^acgtu]//g;
    $seq_body =~ tr/t/u/;

    # Final cleaned sequence
    my ($cleaned_seq) = ($seq_body =~ /([acgu]+)/);
    $cleaned_seq ||= '';
    my $length = length($cleaned_seq);

    # Determine RNA or DNA
    # my $ucount = ($cleaned_seq =~ tr/u//);
    # my $tcount = ($cleaned_seq =~ tr/t//);  # Should be 0, but included for safety
    # my $type = 'unknown';
    # if ($ucount > $tcount) {
    #     $type = 'RNA';
    # } elsif ($tcount > $ucount) {
    #     $type = 'DNA';
    # }
    # $type = 'unknown' if ($ucount > 20 && $tcount > 20);



    return {
        name        => $header,
        cleaned_seq => $cleaned_seq,
        # type        => $type,
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

# if only one sequecen is submitted then forward it to result

if ($all_done && scalar(@job_ids) == 1) {
    my $single_job_id = $job_ids[0];
    print <<JS_REDIRECT;
    <script>
    // Immediately replace current page in browser history
    window.location.replace('../tmp/jobs/job_$single_job_id/result.html');
    </script>
JS_REDIRECT
    exit;
}

# Table display
print "<table>";
my $elapsed = $q->param('elapsed_time') || 0;
print "<div class='timing-info'>\n";
print "<span id='timer' style='display:none'>$elapsed</span>\n" unless $all_done;   # removed print statement to remove showing elapsed time (seconds) stored in $elapsed
print "<div class='estimated'>Max. Estimated Run Time: < $formatted_time minutes</div>\n" unless $all_done;
print "</div>\n";
print "<br>";


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
} else {
    print "<p class='waiting-msg'> Waiting for jobs to finish...</p>";
}

# Hidden form to keep job IDs across refreshes
print "<form id='refreshForm' method='POST' action='batchsubmit.cgi'>\n";
foreach my $jid (@job_ids) {
    print "<input type='hidden' name='job_ids' value='$jid'>\n";
}
print "<input type='hidden' id='elapsed_time' name='elapsed_time' value='$elapsed'>\n";
print "<input type='hidden' name='formatted_time' value='$formatted_time'>\n";
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
