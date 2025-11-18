# Use an official Ubuntu base image
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies, including libxml2 libraries for XML::LibXML
RUN apt-get update && apt-get install -y \
    fcgiwrap \
    perl \
    cpanminus \
    spawn-fcgi \
    git \
    wget \
    unzip \
    make \
    automake \
    build-essential \
    libncurses5-dev \
    libxml2 \
    libxml2-dev \
    libssl-dev \
    libc6-i386 \
    openssl \
    python3 \
    python3-pip \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Clone the rnaanalyzer project into /var/www/rnaanalyzer
RUN git clone https://github.com/Department-of-Bioinformatics/rnaanalyzer.git /var/www/rnaanalyzer

# Install prerequisites for BioPerl-related modules (including XML::LibXML)
RUN cpanm --notest \
    XML::LibXML \
    Bio::DB::DBFetch \
    Bio::DB::WebDBSeqI \
    Bio::Root::Version \
    Bio::SearchIO \
    Bio::Root::IO \
    Bio::SeqIO

# Install additional BioPerl dependencies, skipping tests to avoid failure from missing optional modules
RUN cpanm --notest --force \
    Bio::DB::EMBL \
    Bio::DB::GenBank \
    Bio::DB::GenPept \
    Bio::DB::RefSeq \
    Bio::DB::SwissProt \
    Bio::Root::Test \
    Bio::Tools::Run::RemoteBlast || true

    
# Install BioPerl itself, using --notest to skip any problematic tests
RUN cpanm --notest --force Bio::Perl

# Install the CGI module (if needed)
RUN cpanm --notest CGI

# Install required Perl modules
RUN cpan CGI

# Biopython
RUN pip3 install biopython


# --- Install backend tools in a suitable subdirectory ---
# Here we assume you want to install backend programs under /var/www/rnaanalyzer/bin.
RUN mkdir /var/www/rnaanalyzer/tmp

## Installing gene-scan
RUN cd /var/www/rnaanalyzer/bin && \
    wget -O /var/www/rnaanalyzer/bin/genscan.zip "https://www.coreunitrdm.biozentrum.uni-wuerzburg.de/index.php/s/ZspAjbN4nqmpTwd/download" && \
    unzip /var/www/rnaanalyzer/bin/genscan.zip -d /var/www/rnaanalyzer/bin && \
    rm /var/www/rnaanalyzer/bin/genscan.zip && \
    chmod a+x /var/www/rnaanalyzer/bin/genscanlinux && \
    chmod a+r /var/www/rnaanalyzer/bin/genscanlinux/*.smat
RUN chmod a+x /var/www/rnaanalyzer/bin/genscanlinux/genscan

## Installing tRNAscan-SE:


## Installing tRNAscan-SE:
RUN mkdir /tmp/tRNAscan-SE 
RUN git clone https://github.com/UCSC-LoweLab/tRNAscan-SE.git /tmp/tRNAscan-SE && \
    cd /tmp/tRNAscan-SE && \
    autoreconf -fi && \
    ./configure --prefix=/var/www/rnaanalyzer/bin/tRNAscan-SE && \
    make && \
    mkdir -p /var/www/rnaanalyzer/bin/tRNAscan-SE/bin && \
    cp tRNAscan-SE /var/www/rnaanalyzer/bin/tRNAscan-SE/bin/ && \
    cp tRNAscan-SE.conf /var/www/rnaanalyzer/bin/tRNAscan-SE/bin/ && \
    cp -R lib/ /var/www/rnaanalyzer/bin/tRNAscan-SE/

# Set PERL5LIB so Perl can locate tRNAscan-SE's Perl modules.
# This makes Perl search in /var/www/rnaanalyzer/bin/tRNAscan-SE/lib,
# which contains the tRNAscanSE/Configuration.pm module.
ENV PERL5LIB="/var/www/rnaanalyzer/bin/tRNAscan-SE/lib:$PERL5LIB"
    


# Install Infernal 1.1.5 and create symbolic links for cmsearch, cmscan, and cmstat
RUN cd /var/www/rnaanalyzer/bin && \
    wget http://eddylab.org/software/infernal/infernal.tar.gz && \
    tar zxf infernal.tar.gz && \
    cd infernal-1.1.5 && \
    ./configure --prefix=/var/www/rnaanalyzer/bin/tRNAscan-SE/ && \
    make && \
    cd /var/www/rnaanalyzer/bin/tRNAscan-SE/bin && \
    ln -s ../infernal-1.1.5/src/cmsearch cmsearch && \
    ln -s ../infernal-1.1.5/src/cmscan cmscan && \
    ln -s ../infernal-1.1.5/src/cmstat cmstat


## Install ViennaRNA 2.7.0 from source in /var/www/rnaanalyzer/bin
RUN cd /var/www/rnaanalyzer/bin && \
    wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_7_x/ViennaRNA-2.7.0.tar.gz && \
    tar -xzf ViennaRNA-2.7.0.tar.gz && \
    cd ViennaRNA-2.7.0 && \
    ./configure && \
    make && \
    make install


RUN cd /var/www/rnaanalyzer/bin/ViennaRNA-2.7.0/interfaces/Perl && \
    perl Makefile.PL && \
    make && \
    make test && \
    make install
    

## Set up the Rfam database
RUN mkdir -p /var/www/rnaanalyzer/databases/rfam && \
    cd /var/www/rnaanalyzer/databases/rfam && \
    wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz && \
    gunzip Rfam.cm.gz


## Install CPC2 standalone and build its libsvm component
RUN mkdir -p /var/www/rnaanalyzer/bin/cpc2 && \
    cd /var/www/rnaanalyzer/bin/cpc2 && \
    wget https://github.com/gao-lab/CPC2_standalone/archive/refs/tags/v1.0.1.tar.gz && \
    gunzip -c v1.0.1.tar.gz | tar -xvf - && \
    cd CPC2_standalone-1.0.1/libs/libsvm && \
    gzip -dc libsvm-3.18.tar.gz | tar xf - && \
    cd libsvm-3.18 && \
    make clean && make

    
    
# Permissons
## Adjust file permissions for the RNA Analyzer directory
RUN chown -R www-data:www-data /var/www/rnaanalyzer && \
    find /var/www/rnaanalyzer -type d -exec chmod 755 {} \; && \
    find /var/www/rnaanalyzer -type f -exec chmod 644 {} \;


# Code fixes
## Genscan & Vienna absolute path
RUN find /var/www/rnaanalyzer/cgi-bin -type f \( -name "*.cgi" -or -name "*.pl" \) \
    -exec sed -i 's|\.\./bin/genscanlinux|/var/www/rnaanalyzer/bin/genscanlinux|g' {} + && \
RUN find /var/www/rnaanalyzer/cgi-bin -type f \( -name "*.cgi" -or -name "*.pl" \) \
    -exec sed -i 's|\.\./bin/ViennaRNA-2.6.4/src/bin|/var/www/rnaanalyzer/bin/ViennaRNA-2.7.0/src/bin|g' {} +

## Vienna: Update paths in the CGI scripts according to instructions
RUN find /var/www/rnaanalyzer/cgi-bin -type f \( -name "*.cgi" -or -name "*.pl" \) \
    -exec sed -i 's|/path/to/rnaanalyzer|/var/www/rnaanalyzer|g' {} + && \
    find /var/www/rnaanalyzer/cgi-bin -type f \( -name "*.cgi" -or -name "*.pl" \) \
    -exec sed -i 's|/var/www/rnaanalyzer/bin/ViennaRNA-1.5/Progs|/var/www/rnaanalyzer/bin/ViennaRNA-2.4.18/src/bin/|g' {} + && \
    find /var/www/rnaanalyzer/cgi-bin -type f \( -name "*.cgi" -or -name "*.pl" \) \
    -exec sed -i 's|/var/www/rnaanalyzer/bin/ViennaRNA-1.5/Utils|/var/www/rnaanalyzer/bin/ViennaRNA-2.4.18/src/Utils|g' {} +
    
## Vienna: Create the expected ViennaRNA-1.5/Progs directory and symlink the executables from ViennaRNA-2.7.0
RUN mkdir -p /var/www/rnaanalyzer/bin/ViennaRNA-1.5/Progs && \
    ln -s /var/www/rnaanalyzer/bin/ViennaRNA-2.7.0/src/bin/RNAsubopt /var/www/rnaanalyzer/bin/ViennaRNA-1.5/Progs/RNAsubopt && \
    ln -s /var/www/rnaanalyzer/bin/ViennaRNA-2.7.0/src/bin/RNAplot /var/www/rnaanalyzer/bin/ViennaRNA-1.5/Progs/RNAplot
    
### Patch the CGI script to remove the unsupported '-o' option from RNAplot calls
RUN sed -i 's/-o //g' /var/www/rnaanalyzer/cgi-bin/webserver_AA.cgi

### Vienna Version
RUN sed -i 's|\.\./bin/ViennaRNA-2\.6\.4/src/bin|\.\./bin/ViennaRNA-2.7.0/src/bin|g' /var/www/rnaanalyzer/cgi-bin/webserver_AA.cgi


## Replace $TEMPDIR with /tmp in the CGI script
RUN sed -i 's/\$TEMPDIR\//\/tmp\//g' /var/www/rnaanalyzer/cgi-bin/webserver_AA.cgi

## Refien RNAplot
RUN sed -i 's|RNAplot svg|RNAplot --filename-full|g' /var/www/rnaanalyzer/cgi-bin/webserver_AA.cgi

## Vienna: Fix to remove duplicate filenmae-full
RN sed -i 's/--filename-full[[:space:]]\+--filename-full/--filename-full/g' /var/www/rnaanalyzer/cgi-bin/webserver_AA.cgi

## Vienna fix svg html query for tmp directory
sed -i "s|src='\$TEMPDIR/\${job}_ss\.svg'|src='/tmp/\${job}_ss.svg'|g" /var/www/rnaanalyzer/cgi-bin/webserver_AA.cgi


# Expose the port for fcgiwrap
EXPOSE 9000

# Start fcgiwrap on port 9000. 
# (In a production setup, you might use a process manager to also start other processes.)
CMD ["/usr/sbin/fcgiwrap", "-s", "tcp:0.0.0.0:9000"]

RUN chmod 755 /var/www/rnaanalyzer/cgi-bin/webserver_AA.cgi
RUN chown -R www-data:www-data /var/www/rnaanalyzer

