# Following NAR guideline

Requies the reults to be bookmarked. current implementation does not allow this.
this is done now.

# A new html landing page which is more modern and has sufficient information

Also implementated

# Adding security checks to submitted form:

1. Standardized FASTA block format
    read_uploaded_fasta() and read_pasted_fasta() cna now take input as file or text and then run sanitize on them

2. Fixed sequence parsing logic
    sanitize_sequence() correctly handles inputs with or without FASTA headers.

3. Strips leading/trailing whitespace
    Removes blank lines and space noise before parsing headers and sequences.

4. Sanitizes headers
    Removes unsafe characters from sequence names (shell-safe).
    Truncates headers to a max of 20 characters.
    Blocks path traversal (e.g., ../).

5. Sanitizes sequence body
    Removes all non-letter characters.
    Converts everything to lowercase.
    Converts t → u to standardize to RNA.

6. Validates biological content
    Requires cleaned sequence to be at least 10 bases long.
    Requires at least 90% of characters to be valid (a, c, g, t, u).

7.  Blocks invalid sequences
    Skips job creation if a sequence fails the validity check.

8. Rejects overlong sequences
    Automatically dies if input exceeds 20,000 characters (anti-abuse safeguard).

9. Detects and blocks empty submissions
    Displays error if both the text area and file upload are empty.

10. Client-side validation with JavaScript
    Prevents form submission if no sequence is pasted and no file is uploaded.

11. Safe HTML output
    Uses escapeHTML() when displaying user inputs like sequence names to prevent XSS and broken UI.

12. Limit number of sequences to be 5
13. Added 1 mb uplaod limit for fasta file (given ou curent implementation, bigger input is not required)

# added new routine for UTR prediction

this can predict UTR not only from augustus but also from CPC2 if there is coding transcript.
also done now

# adding rbs-finder
    check https://github.com/deprekate/rbs-finder/blob/master/rbs_finder.pl
    if it can be added to above additionof UTR prediction, would be nice

# get a nucleotide composition overview

# check if kraken 2 can be added for transcript claddification 


# Making a liscense folder and also a link to every liscene we use

# converting to python

with the dealine being pushed, we shall focus on python conversion now.

# Batch Submit is working but is quite limited with perl cgi

continue with python and flask (gives more robustness to job handling for users)

# layout for the output result file 

this is curently being done by Yash