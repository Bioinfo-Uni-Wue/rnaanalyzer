#!/usr/bin/perl

use lib "."; 
use lib "./RNASERVER/";
use strict; 
use warnings;
use CGI;
use File::Slurp;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use RNASERVER::JobUtil qw(get_next_job_id);
use JSON;

my $q = CGI->new;
$| = 1;

# Read form inputs
my $input     = $q->param("SEQUENCE") // '';
my $is_rna    = $q->param("ORIGIN") // '';
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

print $q->header();

print <<'HTML';
<html><head>
  <title>Batch Results</title>
  <link rel="stylesheet" href="/css/batchresults.css">
</head><body>
<header>
    <a href="localhost">    <!--- change after putting on server -->
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


# On first submission, launch jobs
unless ($is_refresh) {
    my @sequences = split(/^>/m, $input);
    shift @sequences if $sequences[0] !~ /\S/;

    foreach my $seq (@sequences) {
        next unless $seq =~ /\S/;

        my $job_id = get_next_job_id();
        my $job_dir = "../tmp/jobs/job_$job_id";
        mkdir $job_dir;

        write_file("$job_dir/input.txt", ">$seq");

        my %params = (
            job_id     => $job_id,
            ORIGIN     => $is_rna,
            IRE        => $ire,
            TRANS      => $trans,
            RNAmotif   => $rnamo,
            mirna      => $mirna,
            trna       => $trna,
            run_coding => $coding,
            species    => $species,
            mirna_target => $mirna_target,
        );

        write_file("$job_dir/params.json", encode_json(\%params));
        system("perl ./webserver_AA.cgi $job_id > /dev/null 2>&1 &");

        push @job_ids, $job_id;
    }
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

print "<tr><th>#</th><th>Job ID</th><th>Status</th></tr>";

for (my $i = 0; $i < @job_ids; $i++) {
    my $jid = $job_ids[$i];
    my $result_file = "../tmp/jobs/job_$jid/result.txt";
   
       my $status = -e $result_file
        ? "<a class='status-button' href='../tmp/jobs/job_$jid/result.html' target='_blank'>View Result</a>"
        : "Processing";

    print "<tr><td>" . ($i + 1) . "</td><td>$jid</td><td>$status</td></tr>";
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
