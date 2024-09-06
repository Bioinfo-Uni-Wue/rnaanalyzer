
# RNA Analyzer Installation Guide

## Source Directory

All files for RNA Analyzer will reside in:
```
/var/www/rnaanalyzer
```

## Nginx Configuration

Nginx configuration is located at:
nginx.conf in rnaalyer folder
or/and copied to 
```
/etc/nginx/sites-available/rnaanalyzer.bioinfo-wuerz.de
```

### Installing Nginx

1. Copy the Nginx configuration file to the Nginx `sites-available` directory.
2. Adjust the necessary script lines according to your server setup.

## Dependencies Installation (Ubuntu)

Run the following commands to install necessary dependencies:
```bash
sudo apt update
sudo apt install fcgiwrap perl cpanminus spawn-fcgi
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
- **CGI** (Core)
- **Bio::Tools::Genscan** (Part of BioPerl)
- **Bio::SeqIO** (Part of BioPerl)
- **File::Temp** (Core)
- **File::Basename** (Core)
- **RNASERVER::TRANS2** (Custom)
- **RNASERVER::IRE** (Custom)

### Installing Perl Packages

To install the necessary Perl modules, use the following commands:
```bash
cpan CGI
cpan Bio::Perl
cpanm Bio::Perl
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
tar -xzf genscan.tar.gz
cd genscan
make
```

### tRNAscan-SE

Clone the source and install:
```bash
git clone https://github.com/UCSC-LoweLab/tRNAscan-SE.git
cd tRNAscan-SE/
./configure --prefix=/path/rnaanalyzer/bin/tRNAscan-SE/
automake --add-missing
sudo make
cp tRNAscan-SE bin/
cp tRNAscan-SE.conf bin/
```

### Infernal (Required for tRNAscan-SE)

Install Infernal:
```bash
wget http://eddylab.org/software/infernal/infernal.tar.gz
tar zxf infernal.tar.gz
cd infernal-1.1.5
./configure --prefix=/path/rnaanalyzer/bin/tRNAscan-SE/
make
```

Create symbolic links:
```bash
cd /path/rnaanalyzer/bin/tRNAscan-SE/bin
ln -s ../infernal-1.1.5/src/cmsearch cmsearch
ln -s ../infernal-1.1.5/src/cmscan cmscan
ln -s ../infernal-1.1.5/src/cmstat cmstat
```

### ViennaRNA Package

Install via package manager:
```bash
sudo apt install viennarna
```

If not available, install from source:
```bash
wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_6_x/ViennaRNA-2.4.14.tar.gz
tar -xzf ViennaRNA-2.4.14.tar.gz
cd ViennaRNA-2.4.14
./configure
make
sudo make install
```

## Set Permissions for RNA Analyzer Directory

Set ownership and permissions for the RNA Analyzer files:
```bash
sudo chown -R www-data:www-data /var/www/rnaanalyzer
sudo find /var/www/rnaanalyzer -type d -exec chmod 755 {} ;
sudo find /var/www/rnaanalyzer -type f -exec chmod 644 {} ;
```

## Adjust Paths in RNA Analyzer Scripts

Run the following commands to update paths in the RNA Analyzer scripts:

1. Adjust the base path:
```bash
find /var/www/rnaanalyzer/cgi-bin -type f \( -name "*.cgi" -or -name "*.pl" \) -exec sed -i 's|/storage/srv/bioapps/rnaanalyzer|/var/www/rnaanalyzer|g' {} +
```

2. Update ViennaRNA paths:
```bash
find /var/www/rnaanalyzer/cgi-bin -type f \( -name "*.cgi" -or -name "*.pl" \) -exec sed -i 's|/var/www/rnaanalyzer/bin/ViennaRNA-1.5/Progs|/var/www/rnaanalyzer/bin/ViennaRNA-2.4.18/src/bin/|g' {} +
find /var/www/rnaanalyzer/cgi-bin -type f \( -name "*.cgi" -or -name "*.pl" \) -exec sed -i 's|/var/www/rnaanalyzer/bin/ViennaRNA-1.5/Utils|/var/www/rnaanalyzer/bin/ViennaRNA-2.4.18/src/Utils|g' {} +
```

## Adjust SVG File Path

In the file `rnaanalyzer/cgi-bin/webserver_AA.cgi`, update the path for SVG image generation:
```perl
print "<br><img src='/tmp/"."$job"."_"."ss".".svg' width='452' height='650' alt='RNA Structure'<br>";
```
