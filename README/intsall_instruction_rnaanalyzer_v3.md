# Installation RNA Analyzer

## Source

/var/www/rnaanalyzer

## Nginx config in

/etc/nginx/sites-available/rnaanalyzer.bioinfo-wuerz.de

### install nginx

copy nginx file to your nginx/site-available folder and adjust the script lines.



## Dependencies to install in Ubuntu

```bash
sudo apt update
sudo apt install fcgiwrap perl
sudo apt install cpanminus
sudo apt install spawn-fcgi # for adding this to distro service
```


## Config fcgiwrap 
open port 9000
test if 9000 port is open and used
```bash
sudo iptables -L | grep 9000
nc -zv 116.202.226.17 9000
sudo ss -tuln | grep 9000
```

edit fcgiwrap config
```bash
sudo nano /etc/systemd/system/fcgiwrap.service
```

enter in the file to listen to port 9000
````
[Unit]
Description=Simple CGI Server
After=nss-user-lookup.target

[Service]
ExecStart=/usr/sbin/fcgiwrap -s tcp:0.0.0.0:9000
#ExecStart=/usr/sbin/fcgiwrap
#ExecStart=/usr/sbin/fcgiwrap -s unix:/var/run/fcgiwrap.socket
User=www-data
Group=www-data
#StandardInput=socket
````


## Configure nginx for fcgiwrap

appended config file 

````
sudo nano /etc/nginx/sites-available/rnaanalyzer.bioinfo-wuerz.de
````



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

in ./rnaanalyzer/bin


### Genscan

Example of extracting and installing Genscan:
#Dont know download path of genescan

```bash
tar -xzf genscan.tar.gz
cd genscan
make
```

### tRNAscan-SE

Or install from source:

```bash
git clone https://github.com/UCSC-LoweLab/tRNAscan-SE.git
cd tRNAscan-SE/
./configure --prefix=/path/rnaanalyzer/bin/tRNAscan-SE/
automake --add-missing
sudo make
cp tRNAscan-SE bin/
cp tRNAscan-SE.conf bin/
```

### Install internal (required for tRNAscan)

```
wget http://eddylab.org/software/infernal/infernal.tar.gz
tar zxf infernal.tar.gz
cd infernal-1.1.5
./configure --prefix=/path/rnaanalyzer/bin/tRNAscan-SE/
cd tRNAscan-SE/bin
ln -s ../infernal-1.1.5/src/cmsearch cmsearch
ln -s ../infernal-1.1.5/src/cmscan cmscan
ln -s ../infernal-1.1.5/src/cmstat cmstat
```

### ViennaRNA Package

Install via package manager:

```bash
sudo apt install viennarna
#apt doesnt know viennarna
```

Or download and install from source:

```bash
wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_6_x/ViennaRNA-2.4.14.tar.gz
tar -xzf ViennaRNA-2.4.14.tar.gz
cd ViennaRNA-2.4.14
./configure
make
sudo make install
```

## ownerhsip of rnaanlyzer folder

sudo chown -R www-data:www-data /var/www/rnaanalyzer
sudo find /var/www/rnanalyzer -type d -exec chmod 755 {} ;
sudo find /var/www/rnanalyzer -type f -exec chmod 644 {} ;


## adjust the path in the rnaanlyzeer files

- find /var/www/rnaanalyzer/cgi-bin -type f \\( -name "\*.cgi" -or -name "\*.pl" \\) -exec sed -i 's|/storage/srv/bioapps/rnaanalyzer|/var/www/rnaanalyzer|g' {} +
- 
- grep -rl '/storage/srv/bioapps/rnaanalyzer/' /var/www/rnaanalyzer/cgi-bin
- grep -rl '/var/www/rnaanalyzer' /var/www/rnaanalyzer/cgi-bin
- 
- grep -rl '/var/www/rnaanalyzer' /var/www/rnaanalyzer/cgi-bin
- 

/var/www/rnaanalyzer/bin/ViennaRNA-1.5/Progs

change to

/var/www/rnaanalyzer/bin/ViennaRNA-2.4.18/src/bin/

- find /var/www/rnaanalyzer/cgi-bin -type f \\( -name "\*.cgi" -or -name "\*.pl" \\) -exec sed -i 's|/var/www/rnaanalyzer/bin/ViennaRNA-1.5/Progs|/var/www/rnaanalyzer/bin/ViennaRNA-2.4.18/src/bin/|g' {} +

/var/www/rnaanalyzer/bin/ViennaRNA-1.5/Utils

change to

/var/www/rnaanalyzer/bin/ViennaRNA-2.4.18/src/Utils

- find /var/www/rnaanalyzer/cgi-bin -type f \\( -name "\*.cgi" -or -name "\*.pl" \\) -exec sed -i 's|/var/www/rnaanalyzer/bin/ViennaRNA-1.5/Utils|/var/www/rnaanalyzer/bin/ViennaRNA-2.4.18/src/Utils|g' {} +
- grep -rl 'var/www/rnaanalyzer/bin/ViennaRNA-2.4.18/src/bin/' /var/www/rnaanalyzer/cgi-bin

correct he path

- find /var/www/rnaanalyzer/cgi-bin -type f \\( -name "\*.cgi" -or -name "\*.pl" \\) -exec sed -i 's|/var/www/rnaanalyzer/bin/ViennaRNA-2.4.18/src/bin/|/var/www/rnaanalyzer/bin/ViennaRNA-2.4.18/src/bin|g' {} +


## adjust path to svg file 

in rnaanalyzer/cgi-bin/webserver_AA.cgi
print "<br><img src='$TEMPDIR/"."$job"."_"."ss".".svg' width='452' height='650' alt='RNA Structure'<br>";
to
print "<br><img src='/tmp/"."$job"."_"."ss".".svg' width='452' height='650' alt='RNA Structure'<br>";
