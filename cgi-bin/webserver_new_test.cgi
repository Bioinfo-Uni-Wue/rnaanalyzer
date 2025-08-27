#!/usr/bin/perl

use strict;
use warnings;
use CGI;

# Create a new CGI object
my $cgi = CGI->new;

# Print the HTTP header and start the HTML output
print $cgi->header('text/html');
print "<html><head><title>RNA Analyzer Test</title></head><body>";

# Retrieve the parameters from the form
my $seqname = $cgi->param('SEQNAME') || 'No name provided';
my $sequence = $cgi->param('SEQUENCE') || 'No sequence provided';

# Display the received parameters
print "<h1>RNA Analyzer Test Results</h1>";
print "<p><strong>Sequence Name:</strong> $seqname</p>";
print "<p><strong>Sequence:</strong></p>";
print "<pre>$sequence</pre>";

# End the HTML
print "</body></html>";

