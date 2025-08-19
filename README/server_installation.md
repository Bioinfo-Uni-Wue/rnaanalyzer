## Dependencies Installation (Ubuntu)

Run the following commands to install necessary dependencies:

```bash
sudo apt update
sudo apt install fcgiwrap perl cpanminus spawn-fcgi python
```

# some more libraries if already not installed 
```
sudo apt install autoconf automake libtool
sudo apt-get install libmpfr4
sudo apt-get install libmpfr-dev
sudo apt install perl-doc  # important to check perl module information
```

### Installing Perl Packages

To install the necessary Perl modules, use the following commands:

```bash
cpan CGI
cpan Bio::Perl
cpanm Bio::Perl  #if fcapnm is available
cpan Bio::Tools::Genscan (Part of BioPerl)
cpan Bio::SeqIO (Part of BioPerl)
cpan File::Temp (Core)
cpan File::Basename (Core)
RNASERVER::TRANS2 (Custom) # no installation needed
RNASERVER::IRE (Custom)   # no installation needed
sudo cpan install JSON 
sudo cpan install File::Slurp  
```

### Verifying Perl Module Installation

Test the installation of the Perl modules with the following script:

```perl
use strict;
use warnings;
use CGI;
use Bio::Tools::Genscan;
use Bio::SeqIO;
use File::Temp;
use File::Basename;

print "All modules loaded successfully.";

```

## important install biopython systemwide

```bash
sudo pip3 install biopython
sudo pip3 install pandas
```

## Software Installation

### tRNAscan-SE

Clone the source and install:

```bash
mkdir /tmp/tRNAscan-SE 
## Installing tRNAscan-SE:
git clone https://github.com/UCSC-LoweLab/tRNAscan-SE.git /tmp/tRNAscan-SE 
cd /tmp/tRNAscan-SE
autoreconf -fi 
./configure --prefix=/storage/srv/bioapps/rnaanalyzer/bin/tRNAscan-SE 
make 
mkdir -p /storage/srv/bioapps/rnaanalyzer/bin/tRNAscan-SE/bin
cp tRNAscan-SE /storage/srv/bioapps/rnaanalyzer/bin/tRNAscan-SE/bin/
cp tRNAscan-SE.conf /storage/srv/bioapps/rnaanalyzer/bin/tRNAscan-SE/bin/
cp -R lib/ /storage/srv/bioapps/rnaanalyzer/bin/tRNAscan-SE/

# Set PERL5LIB so Perl can locate tRNAscan-SE's Perl modules
export PERL5LIB="/storage/srv/bioapps/rnaanalyzer/bin/tRNAscan-SE/lib:$PERL5LIB"
```

### Infernal (Required for tRNAscan-SE)

Install Infernal:

```bash
cd /storage/srv/bioapps/rnaanalyzer/bin
wget http://eddylab.org/software/infernal/infernal.tar.gz
tar -zxf infernal.tar.gz
cd infernal-1.1.5
./configure --prefix=/storage/srv/bioapps/rnaanalyzer/bin/tRNAscan-SE/
make
```

Create symbolic links:

```bash
cd /storage/srv/bioapps/rnaanalyzer/bin/tRNAscan-SE/bin
ln -s ../infernal-1.1.5/src/cmsearch cmsearch
ln -s ../infernal-1.1.5/src/cmscan cmscan
ln -s ../infernal-1.1.5/src/cmstat cmstat
```

### ViennaRNA Package

```bash
cd /storage/srv/bioapps/rnaanalyzer/bin/
wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_7_x/ViennaRNA-2.7.0.tar.gz
tar -xzf ViennaRNA-2.7.0.tar.gz
cd ViennaRNA-2.7.0
./configure  
make
make install
```

## RFAM database (uses infernal from tRNAscan folder)

```
mkdir /storage/srv/bioapps/rnaanalyzer/databases/
mkdir /storage/srv/bioapps/rnaanalyzer/databases/rfam
cd ../databases/rfam
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz
cmpress Rfam.cm
```

## CPC2_standalone CPC2 (requires python3 and bipython)

```
mkdir /storage/srv/bioapps/rnaanalyzer/bin/cpc2
wget https://github.com/gao-lab/CPC2_standalone/archive/refs/tags/v1.0.1.tar.gz
gunzip v1.0.1.tar.gz | tar -xvf
cd libs/libsvm/
gzip -dc libsvm-3.18.tar.gz | tar xf -
cd libsvm-3.18/
make clean && make
```
## augustus

# requires lots of libraries
```
git clone https://github.com/Gaius-Augustus/Augustus.git
cd Augutus
make augutus 
```

OR
system-wide install

```
sudo apt install augustus
```

# installing InstaRNA

```
https://github.com/BackofenLab/IntaRNA.git

./configure --prefix=/path/to/IntaRNA --with-vrna=/path/to/ViennaRNA-2.7.0 --disable-pkg-config

make
make install
```
