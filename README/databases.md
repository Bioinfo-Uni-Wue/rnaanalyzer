# get teh relevant databases

# Rfam

```
mkdir /var/www/rnaanalyzer/databases/
mkdir /var/www/rnaanalyzer/databases/rfam
cd ../databases/rfam
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz
cmpress Rfam.cm
```

# miRbase

```
mkdir /var/www/rnaanalyzer/databases/mirbase
wget https://www.mirbase.org/download/hairpin.fa # for miRNA scanning
wget https://www.mirbase.org/download/mature.fa # for miRNA target prediction
```
