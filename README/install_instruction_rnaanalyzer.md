# RNA Analyzer Installation Guide

## Source Directory

All files for RNA Analyzer will reside in:

```
/var/www/rnaanalyzer
```

## Download Project

git pull https://github.com/Department-of-Bioinformatics/rnaanalyzer.git

## Adjust rights and permissions

### Installing Nginx

sudo apt install nginx
sudo service nginx start
sudo service nginx status

## Nginx Configuration

Nginx configuration is located at:
nginx.conf in rnaalyer folder
or/and copied to

```
/etc/nginx/sites-available/rnaanalyzer.bioinfo-wuerz.de
```

1. Copy the Nginx configuration file to the Nginx `sites-available` directory.
2. Adjust the necessary script lines according to your server setup.

## Dependencies Installation (Ubuntu)

Run the following commands to install necessary dependencies:

```bash
sudo apt update
sudo apt install fcgiwrap perl cpanminus spawn-fcgi
```

### Installing Perl Packages

To install the necessary Perl modules, use the following commands:

```bash
cpan CGI
cpan Bio::Perl
cpanm Bio::Perl  #if fcapnm is available 
```

## Configuring fcgiwrap

### Step 1: Open Port 9000

To verify that port 9000 is open and being used, run:

```bash
sudo iptables -L | grep 9000
nc -zv <your_server_ip> 9000
sudo ss -tuln | grep 9000
```

### Step 2: Edit fcgiwrap Configuration

Edit the `fcgiwrap.service` file:

```bash
sudo nano /etc/systemd/system/fcgiwrap.service
```

Update the service file to listen on port 9000:

```ini
[Unit]
Description=Simple CGI Server
After=nss-user-lookup.target

[Service]
ExecStart=/usr/sbin/fcgiwrap -s tcp:0.0.0.0:9000
User=www-data
Group=www-data
```

## Nginx Configuration for fcgiwrap

To configure Nginx for fcgiwrap, edit the Nginx configuration file:

```bash
sudo nano /etc/nginx/sites-available/rnaanalyzer.bioinfo-wuerz.de
```

Append the necessary configurations to enable CGI processing via fcgiwrap.

## Perl Modules

Ensure the following Perl packages are installed:

```
    cpan CGI (Core)
    cpan Bio::Tools::Genscan** (Part of BioPerl)
    cpan Bio::SeqIO** (Part of BioPerl)
    cpan File::Temp** (Core)
    cpan File::Basename** (Core)
    RNASERVER::TRANS2** (Custom) # no installation needed
    RNASERVER::IRE** (Custom)   # no installation needed
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

print "All modules loaded successfully.
";
```

If custom modules are in a non-standard location, specify the path as follows:

```perl
use lib '/path/to/your/modules';
```

## Software Installation

### Genscan

Extract and install Genscan:

```bash
cd /var/www/rnaanalyzer/bin 
    wget -O /var/www/rnaanalyzer/bin/genscan.zip "https://www.coreunitrdm.biozentrum.uni-wuerzburg.de/index.php/s/ZspAjbN4nqmpTwd/download" 
unzip /var/www/rnaanalyzer/bin/genscan.zip -d /var/www/rnaanalyzer/bin
rm /var/www/rnaanalyzer/bin/genscan.zip
chmod a+x /var/www/rnaanalyzer/bin/genscanlinux
chmod a+r /var/www/rnaanalyzer/bin/genscanlinux/*.smat
```

### tRNAscan-SE

Clone the source and install:

```bash
mkdir /tmp/tRNAscan-SE 
## Installing tRNAscan-SE:
git clone https://github.com/UCSC-LoweLab/tRNAscan-SE.git /tmp/tRNAscan-SE 
cd /tmp/tRNAscan-SE
autoreconf -fi 
./configure --prefix=/var/www/rnaanalyzer/bin/tRNAscan-SE 
make 
mkdir -p /var/www/rnaanalyzer/bin/tRNAscan-SE/bin
cp tRNAscan-SE /var/www/rnaanalyzer/bin/tRNAscan-SE/bin/
cp tRNAscan-SE.conf /var/www/rnaanalyzer/bin/tRNAscan-SE/bin/
cp -R lib/ /var/www/rnaanalyzer/bin/tRNAscan-SE/

# Set PERL5LIB so Perl can locate tRNAscan-SE's Perl modules
export PERL5LIB="/var/www/rnaanalyzer/bin/tRNAscan-SE/lib:$PERL5LIB"
```

### Infernal (Required for tRNAscan-SE)

Install Infernal:

```bash
cd /var/www/rnaanalyzer/bin
wget http://eddylab.org/software/infernal/infernal.tar.gz
tar zxf infernal.tar.gz
cd infernal-1.1.5
./configure --prefix=/var/www/rnaanalyzer/bin/tRNAscan-SE/
make
```

Create symbolic links:

```bash
cd /var/www/rnaanalyzer/bin/tRNAscan-SE/bin
ln -s ../infernal-1.1.5/src/cmsearch cmsearch
ln -s ../infernal-1.1.5/src/cmscan cmscan
ln -s ../infernal-1.1.5/src/cmstat cmstat
```

### ViennaRNA Package

```bash
cd /var/www/rnaanalyzer/bin/
wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_7_x/ViennaRNA-2.7.0.tar.gz
tar -xzf ViennaRNA-2.7.0.tar.gz
cd ViennaRNA-2.7.0
./configure
make
sudo make install
```

#### Set Permissions for RNA Analyzer Directory

Set ownership and permissions for the RNA Analyzer files:

```bash
sudo chown -R www-data:www-data /var/www/rnaanalyzer
sudo find /var/www/rnaanalyzer -type d -exec chmod 755 {} ;
sudo find /var/www/rnaanalyzer -type f -exec chmod 644 {} ;
```

#### Adjust Paths in RNA Analyzer Scripts

Run the following commands to update paths in the RNA Analyzer scripts:

1. Adjust the base path:

```bash
find /var/www/rnaanalyzer/cgi-bin -type f \( -name "*.cgi" -or -name "*.pl" \) -exec sed -i 's|/storage/srv/bioapps/rnaanalyzer|/var/www/rnaanalyzer|g' {} +
```

1. Update ViennaRNA paths:

```bash
find /var/www/rnaanalyzer/cgi-bin -type f \( -name "*.cgi" -or -name "*.pl" \) -exec sed -i 's|/var/www/rnaanalyzer/bin/ViennaRNA-1.5/Progs|/var/www/rnaanalyzer/bin/ViennaRNA-2.4.18/src/bin/|g' {} +
find /var/www/rnaanalyzer/cgi-bin -type f \( -name "*.cgi" -or -name "*.pl" \) -exec sed -i 's|/var/www/rnaanalyzer/bin/ViennaRNA-1.5/Utils|/var/www/rnaanalyzer/bin/ViennaRNA-2.4.18/src/Utils|g' {} +
```

## RFAM

```
mkdir /var/www/rnaanalyzer/databases/
mkdir /var/www/rnaanalyzer/databases/rfam
cd ../databases/rfam
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz
cmpress Rfam.cm
```

## CPC2_standalone CPC2 (requires python3 and bipython)

```
mkdir /var/www/rnaanalyzer/bin/cpc2
wget https://github.com/gao-lab/CPC2_standalone/archive/refs/tags/v1.0.1.tar.gz
gunzip v1.0.1.tar.gz | tar -xvf
cd libs/libsvm/
gzip -dc libsvm-3.18.tar.gz | tar xf -
cd libsvm-3.18/
make clean && make
```

# check python dependencies

## important install biopython systemwide

sudo pip3 install biopython

## Adjust SVG File Path

In the file `/var/www/rnaanalyzer/cgi-bin/webserver_AA.cgi`, update the path for SVG image generation:

> in file /var/www/rnaanalyzer/cgi-bin/webserver_AA.cgi  replace print "<br><img src='$TEMPDIR/"."$job"."\_"."ss".".svg' width='452' height='650' alt='RNA Structure'<br>"; with print "<br><img src='/tmp/"."$job"."\_"."ss".".svg' width='452' height='650' alt='RNA Structure'<br>";

```
sed -i 's/$TEMPDIR///tmp//g' /var/www/rnaanalyzer/cgi-bin/webserver_AA.cgi
```

## add miRNA target scan

```
sudo pip3 install pandas
```
miRanda good tool used globally but really old. but no other good replacement. checked RRNAhybrid as well was much slower. Other known tools for this are species specific like TARGETSCAN. 
Local installion of miRanda was not possible due to old libraries. Instead conda was used to install and then miRanda binary was extracted which works well. 
A python wrapper was made to execute miranda and also parse the output. WOrks well but time is an issue. A full scan takes ~2 minutes which gives a timeout on CGI. Need to think ideas about this.

## augustus

# requires lots of libraries
```
git clone https://github.com/Gaius-Augustus/Augustus.git
cd Augutus
make augutus 
```

might require boost libraries and sql

```
sudo apt install libboost-all-dev
sudo apt install libmysql++-dev
sudo apt install libsqlite3-dev
sudo apt install libgsl-dev
sudo apt install liblpsolve55-dev
```
OR
system-wide install

```
sudo apt install augustus
```
## new commit ##

# some more libraries if already not installed 
```
sudo apt install autoconf automake libtool
sudo apt-get install libmpfr4
sudo apt-get install libmpfr-dev
sudo apt install perl-doc  # important to check perl module information
```

# perl modules
```
sudo cpan install JSON # for reading inputs
sudo cpan install File::Slurp  #for creating and storing correct job folder
```

# updating ViennaRNA to 2.7.0
```
wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_7_x/ViennaRNA-2.7.0.tar.gz
tar -xzf ViennaRNA-2.7.0.tar.gz
cd ViennaRNA-2.7.0

./configure --prefix=/rnaanalyzer/bin/ViennaRNA-2.7.0

make # takes a bit of time needs patience
make install # installs the correct libraries in right location. last time it was not so right!
```

# Installing HMMER3 

```
wget http://eddylab.org/software/hmmer/hmmer.tar.gz
tar zxf hmmer.tar.gz
cd hmmer-3.4
./configure --prefix /your/install/path   # replace /your/install/path with what you want, obv 
make
make check                                # optional: run automated tests
make install                              # optional: install HMMER programs, man pages
```

# adding protein binding motif scan
# Installing FIMO

```
wget https://meme-suite.org/meme/meme-software/5.5.8/meme-5.5.8.tar.gz
tar zxf meme-5.5.8.tar.gz
cd meme-5.5.8
./configure --prefix=$HOME/meme --enable-build-libxml2 --enable-build-libxslt
make
make test
make install
```

# getting rbpdb pfm file 

```
http://rbpdb.ccbr.utoronto.ca/downloads/PFMDir.zip

unzip PFMDir.zip
```

# creating a .meme file using PFMS from RBPDB using custom Python script 

```
change_name_PFM.py
```
