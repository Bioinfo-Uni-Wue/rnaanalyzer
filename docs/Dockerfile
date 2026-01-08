# Use the official Perl image as the base image for Perl dependencies
FROM perl:latest as perl-base

# Install necessary Perl modules
RUN cpanm CGI Bio::Perl

# Use the official Nginx image as the base image for the web server
FROM nginx:latest

# Install fcgiwrap for CGI script handling
RUN apt-get update && apt-get install -y fcgiwrap

# Copy Nginx configuration
COPY nginx.conf /etc/nginx/nginx.conf

# Copy the web application files to the container
COPY rnaanalyzer /var/www/rnaanalyzer
COPY cgi-bin /usr/lib/cgi-bin

# Ensure scripts in cgi-bin are executable
RUN chmod +x /usr/lib/cgi-bin/*

# Expose port 80 to the outside world
EXPOSE 80

# Define the default command to start Nginx
CMD ["nginx", "-g", "daemon off;"]
