#!/usr/bin/perl

use strict;
use warnings;
use CGI;

my $q = CGI->new;

my $job_id = $q->param('job_id');

$job_id =~ s/\s+//g;

my $results_dir = "../tmp/jobs/job_$job_id";
my $result_file = "$results_dir/result.html";

if (-e $result_file) {
    print $q->redirect($result_file);
    exit;
} else {
    print $q->header();
    print "<p>Job number not found. Please check and try again.</p>";
}
