
# Installation RNA Analyzer

## Source

/var/www/rnaanalyzer

## Nginx config in

/etc/nginx/sites-available/rnaanalyzer.bioinfo-wuerz.de

## Dependencies to install in Ubuntu

```bash
sudo apt update
sudo apt install fcgiwrap perl
sudo apt install cpanminus
sudo apt install spawn-fcgi # for adding this to distro service
```

## Perl packages

- CGI
- Bio::Tools::Genscan (part of the Bioperl suite)
- Bio::SeqIO (part of the Bioperl suite)
- File::Temp (core module)
- File::Basename (core module)
- RNASERVER::TRANS2 (custom module)
- RNASERVER::IRE (custom module)

### Install Perl Packages

```bash
cpan CGI
cpan Bio::Perl
cpanm Bio::Perl
```

### Example Perl Script to Verify Modules

```perl
use strict;
use warnings;
use CGI;
use Bio::Tools::Genscan;
use Bio::SeqIO;
use File::Temp;
use File::Basename;

print "All modules loaded successfully.\n";
```

Add custom module paths if necessary:

```perl
use lib '/path/to/your/modules';
```

## Software to install

### Genscan

Example of extracting and installing Genscan:

```bash
tar -xzf genscan.tar.gz
cd genscan
make
```

### tRNAscan-SE

If available via package manager:

```bash
sudo apt install tRNAscan-SE
```

Or install from source:

```bash
git clone https://github.com/lowelab/tRNAscan-SE.git
cd tRNAscan-SE
make
```

### ViennaRNA Package

Install via package manager:

```bash
sudo apt install viennarna
```

Or download and install from source:

```bash
wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.14.tar.gz
tar -xzf ViennaRNA-2.4.14.tar.gz
cd ViennaRNA-2.4.14
./configure
make
sudo make install
```
