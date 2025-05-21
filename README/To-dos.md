# Following NAR guideline

Requies the reults to be bookmarked. current implementation does not allow this.
this is done now.

# A new html landing page which is more modern and has sufficient information

Also implementated

# Adding security checks to submitted form:

    1. Validate All Inputs
    Ensure required fields are present, correctly formatted, and within safe length limits.

    2. Escape Output to Prevent XSS
    Always HTML-escape user data before displaying it in pages.

    3. Sanitize for File and Shell Safety
    Remove dangerous characters from input used in filenames or system calls.

    4. Limit Abuse
    Enforce limits on sequence size, job count, and request rate to prevent overuse.

    5. Suppress Sensitive Errors
    Avoid exposing internal paths or raw input in error messages; log safely.

# converting to python

with the dealine being pushed, we shall focus on python conversion now.

# Batch Submit is working but is quite limited with perl cgi

continue with python and flask (gives more robustness to job handling for users)

# layout for the output result file 

this is curently being done by Yash