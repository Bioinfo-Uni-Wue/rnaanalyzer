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
  <style>
    body {
      font-family: Arial, sans-serif;
      background: #f4f6f8;
      padding: 20px;
    }
    .results-container {
      background: #fff;
      padding: 20px;
      border-radius: 10px;
      box-shadow: 0 0 10px rgba(0,0,0,0.1);
      max-width: 700px;
      margin: auto;
    }
    h2 {
      color: #003c7f;
    }
    table {
      border-collapse: collapse;
      width: 100%;
      margin-top: 15px;
    }
    th, td {
      padding: 12px 15px;
      text-align: center;
    }
    th {
      background-color: #34495e;
      color: white;
    }
    tr:nth-child(even) {
      background-color: #f2f2f2;
    }
    .status-button {
      padding: 6px 10px;
      background-color: #47b92f;
      color: white;
      border: none;
      border-radius: 5px;
      cursor: pointer;
      text-decoration: none;
    }
    .status-button:hover {
      background-color: #2980b9;
    }
    .done-msg {
      margin-top: 20px;
      color: green;
      font-weight: bold;
    }
    .waiting-msg {
      margin-top: 20px;
      color: orange;
    }
  </style>
</head><body>
<div class="results-container">
<h2> Batch Results</h2>
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
print "<tr><th>#</th><th>Job ID</th><th>Status</th></tr>";

for (my $i = 0; $i < @job_ids; $i++) {
    my $jid = $job_ids[$i];
    my $result_file = "../tmp/jobs/job_$jid/result.txt";

    my $status = -e $result_file
        ? "<a class='status-button' href='../tmp/jobs/job_$jid/result.html' target='_blank'>View Result</a>"
        : " Still processing";

    print "<tr><td>" . ($i + 1) . "</td><td>$jid</td><td>$status</td></tr>";
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
print "</form>\n";

# Auto-refresh script if still processing
unless ($all_done) {
    print <<'JS';
<script>
  setTimeout(function() {
    document.getElementById('refreshForm').submit();
  }, 5000);
</script>
JS
}

print "</div></body></html>";
