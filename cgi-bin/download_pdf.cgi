#!/usr/bin/perl
use strict;
use warnings;
use CGI qw(param);
use File::Temp qw(tempfile);
use FindBin qw($Bin);

# -------- CONFIG --------

my $BASE_HREF = "http://localhost/";      # or "/"

my $NODE      = "../bin/node/bin/node";     # <-- adjust to your actual node
my $NODE_PATH = "../bin/node/node_modules";
my $TIMEOUTMS = 60000;
# ------------------------

my $job = param('job') // '';
$job =~ s/[^0-9A-Za-z_-]//g;

if (!$job) {
  print "Status: 400\r\nContent-Type: text/plain\r\n\r\nMissing job parameter.\n";
  exit;
}

my $html_path = "../tmp/jobs/job_$job/result.html";
if (!-f $html_path) {
  print "Status: 404\r\nContent-Type: text/plain\r\n\r\nHTML not found: $html_path\n";
  exit;
}

my $TMPDIR = "../tmp/jobs/job_$job/";

# Ensure tmpdir exists
if (!-d $TMPDIR) {
  print "Status: 500\r\nContent-Type: text/plain\r\n\r\nTMPDIR missing: $TMPDIR\n";
  exit;
}

my ($js_fh, $js_path)   = tempfile("render_pdf_XXXX", SUFFIX => ".js",  DIR => $TMPDIR);
my ($pdf_fh, $pdf_path) = tempfile("results_${job}_XXXX", SUFFIX => ".pdf", DIR => $TMPDIR);
close $pdf_fh;

my $renderer_js = <<'JS';
"use strict";
const fs = require("fs");
const puppeteer = require("puppeteer");

function getArg(name, def=null) {
  const i = process.argv.indexOf(name);
  if (i !== -1 && i + 1 < process.argv.length) return process.argv[i+1];
  return def;
}
function injectBaseHref(html, baseHref) {
  if (!baseHref) return html;
  if (/<base\s/i.test(html)) return html;
  const baseTag = `<base href="${baseHref}">`;
  if (/<head[^>]*>/i.test(html)) return html.replace(/<head[^>]*>/i, m => m + "\n" + baseTag + "\n");
  return baseTag + "\n" + html;
}

(async () => {
  const htmlFile  = getArg("--html");
  const outFile   = getArg("--out");
  const baseHref  = getArg("--base", "");
  const timeoutMs = parseInt(getArg("--timeout", "60000"), 10);

  if (!htmlFile || !outFile) process.exit(2);

  const html = fs.readFileSync(htmlFile, "utf8");
  const html2 = injectBaseHref(html, baseHref);

  const browser = await puppeteer.launch({
    executablePath: process.env.CHROME_PATH,
    args: ["--no-sandbox", "--disable-setuid-sandbox"],
  });

  try {
    const page = await browser.newPage();
    await page.setViewport({ width: 1280, height: 800 });

    await page.setContent(html2, { waitUntil: "domcontentloaded", timeout: timeoutMs });
    await page.waitForFunction(() => document.readyState === "complete", { timeout: timeoutMs });

    await page.emulateMediaType("print");

    // Measure full content height (after render)
    const fullHeight = await page.evaluate(() => {
      const doc = document.documentElement;
      const body = document.body;
      return Math.max(
        doc.scrollHeight, doc.offsetHeight, doc.clientHeight,
        body ? body.scrollHeight : 0,
        body ? body.offsetHeight : 0
      );
    });

    const pdf = await page.pdf({
      width: "300mm",            // or "250mm" etc.
      height: `${fullHeight}px`, // single tall page
      printBackground: true,
      margin: { top: "0mm", right: "0mm", bottom: "0mm", left: "0mm" },
      pageRanges: "1",
      preferCSSPageSize: false,
    });

    fs.writeFileSync(outFile, pdf);
  } catch (e) {
    console.error(e && e.stack ? e.stack : String(e));
    process.exit(1);
  } finally {
    await browser.close();
  }
})();
JS

print $js_fh $renderer_js;
close $js_fh;

local $ENV{NODE_PATH} = $NODE_PATH;
local $ENV{CHROME_PATH} = "path/to/chromium"; # add path to chromium before running


# capture output
my $cmd = qq{"$NODE" "$js_path" --html "$html_path" --out "$pdf_path" --base "$BASE_HREF" --timeout "$TIMEOUTMS" 2>&1};
my $out = `$cmd`;
my $rc  = $?;  # exit status

if ($rc != 0 || !-s $pdf_path) {
  unlink $js_path  if -e $js_path;
  unlink $pdf_path if -e $pdf_path;
  print "Status: 500\r\nContent-Type: text/plain\r\n\r\nPDF generation failed (rc=$rc)\n$out\n";
  exit;
}

open my $pdf, "<:raw", $pdf_path or do {
  unlink $js_path  if -e $js_path;
  unlink $pdf_path if -e $pdf_path;
  print "Status: 500\r\nContent-Type: text/plain\r\n\r\nCould not open PDF.\n";
  exit;
};

print "Content-Type: application/pdf\r\n";
print "Content-Disposition: attachment; filename=\"results_$job.pdf\"\r\n\r\n";
while (read($pdf, my $buf, 8192)) { print $buf; }
close $pdf;

unlink $js_path  if -e $js_path;
unlink $pdf_path if -e $pdf_path;

