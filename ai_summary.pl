#!/usr/bin/env perl
use strict;
use warnings;
use utf8;
use open qw(:std :encoding(UTF-8));
use JSON::PP qw(encode_json decode_json);
use LWP::UserAgent;

# ---- config (override via env or argv[0]) ----
my $facts_file  = $ARGV[0] // $ENV{FACTS_FILE} // 'facts.txt';
my $model       = $ENV{MODEL}      // 'gemma3:1b';
my $ollama_url  = $ENV{OLLAMA_URL} // 'http://127.0.0.1:11434';
my $temperature = defined $ENV{TEMP} ? 0 + $ENV{TEMP} : 0.2;   # deterministic by default
my $num_predict = defined $ENV{NUM_PREDICT} ? 0 + $ENV{NUM_PREDICT} : 120;
my $timeout_s   = defined $ENV{TIMEOUT_S} ? 0 + $ENV{TIMEOUT_S} : 30;
my $top_p       = defined $ENV{TOP_P} ? 0 + $ENV{TOP_P} : 0.7;
my $repeat_pen  = defined $ENV{REPEAT_PENALTY} ? 0 + $ENV{REPEAT_PENALTY} : 1.05;

# ---- read facts ----
my $facts = '';
if (-e $facts_file) {
  local $/;
  open my $fh, '<:encoding(UTF-8)', $facts_file or do { print "AI summary not available\n"; exit 0; };
  $facts = <$fh>;
  close $fh;
}
if (!defined $facts || $facts !~ /\S/) { print "AI summary not available\n"; exit 0; }

# ---- helper: fallback deterministic render (no LLM) ----
sub int_to_words {
  my ($n) = @_;
  my %w = (1=>'One',2=>'Two',3=>'Three',4=>'Four',5=>'Five',6=>'Six',7=>'Seven',8=>'Eight',9=>'Nine',10=>'Ten');
  return $w{$n} // $n;
}
# sub fallback_render {
#   my ($facts) = @_;
#   my ($pos) = $facts =~ /position\s*(\d+)/i;
#   my ($rbp) = $facts =~ /(\d+)\s*RNA-?binding proteins/i;
#   $rbp //= $facts =~ /(\d+)\s*RBPs?/i ? $1 : undef;

#   my $line1 = "The RNA is highly structured and likely regulatory.";
#   my $line2 = defined $pos
#     ? "An IRE motif was detected at position $pos by Rfam."
#     : "An IRE motif was detected by Rfam.";
#   my $line3 = ($facts =~ /coding\s*potential\s*:\s*noncoding/i) ? "The RNA is noncoding." : "The RNA is noncoding.";
#   my $line4 = defined $rbp
#     ? sprintf("%s RNA-binding proteins were detected.", int_to_words($rbp))
#     : "RNA-binding proteins were detected.";

#   return join("\n", $line1, $line2, $line3, $line4);
# }

# ---- strict prompt (maps facts -> fixed sentences; bans hedging) ----
my $prompt = <<"PROMPT";
Summarize the facts without adding any information. Should have a little flow and little context.


Rules:
- Do not infer new biology; no new terms.
- Do not change numbers or positions.
- When the coding potenial is equal to non-coding then just say that the RNA is non-coding in nature.
- Do not use "coding potential".
- Do not use hedging words: suggests, may, might, could, appears, potentially, indicates, possibly.
- expand the facts a little.

FACTS:
$facts
PROMPT

# ---- call Ollama /api/generate (stream=false) ----
my $ua = LWP::UserAgent->new(timeout => $timeout_s, agent => 'ai-summary-deterministic/1.0');
my $res = $ua->post(
  "$ollama_url/api/generate",
  'Content-Type' => 'application/json',
  Content => encode_json({
    model   => $model,
    prompt  => $prompt,
    stream  => JSON::PP::false,
    options => {
      temperature     => $temperature,
      top_p           => $top_p,
      repeat_penalty  => $repeat_pen,
      num_predict     => $num_predict
    }
  })
);

# ---- handle response ----
my $out = '';
if ($res->is_success) {
  my $json = eval { decode_json($res->decoded_content) } || {};
  $out  = $json->{response} // '';
  $out =~ s/^\s+|\s+$//g;
}

# ---- validator to catch drift; fallback if needed ----
# sub validate_summary {
#   my ($facts, $text) = @_;
#   return "empty" unless $text && $text =~ /\S/;

#   my @lines = grep { /\S/ } split /\n/, $text;
#   return "lines" unless @lines == 4;

#   # no hedging
#   my $hedge = qr/\b(suggests?|may|might|could|appears?|potential(?:ly)?|indicat(?:e|es|ing))\b/i;
#   return "hedge" if $text =~ $hedge;

#   # line1 must state highly structured + likely regulatory
#   return "l1" unless $lines[0] =~ /highly structured.*likely regulatory/i;

#   # enforce noncoding on line3
#   return "l3" unless $lines[2] =~ /\bnoncoding\b/i;

#   # extract expectations from facts
#   my ($pos) = $facts =~ /position\s*(\d+)/i;
#   my ($rbp) = $facts =~ /(\d+)\s*RNA-?binding proteins/i;
#   $rbp //= $facts =~ /(\d+)\s*RBPs?/i ? $1 : undef;

#   # line2 must contain IRE motif, position (if present), and Rfam
#   return "l2_rfam" unless $lines[1] =~ /\bIRE motif\b/i && $lines[1] =~ /\bRfam\b/i;
#   if (defined $pos) {
#     return "l2_pos" unless $lines[1] =~ /\bposition\s*$pos\b/i;
#   }

#   # line4 must reflect RBP count (if present)
#   if (defined $rbp) {
#     my %w = (1=>'one',2=>'two',3=>'three',4=>'four',5=>'five',6=>'six',7=>'seven',8=>'eight',9=>'nine',10=>'ten');
#     my $want_word = $w{$rbp} // '';
#     return "l4_rbp" unless ($lines[3] =~ /\b\Q$rbp\E\b/i) || ($want_word && $lines[3] =~ /\b$want_word\b/i);
#   }

#   return ""; # OK
# }

# my $err = validate_summary($facts, $out);
# if ($err ne "") {
#   # fallback to deterministic local render
#   $out = fallback_render($facts);
# }

print(($out =~ /\S/) ? "$out\n" : "AI summary not available\n");
