FROM rocker/r-base

MAINTAINER Paolo Di Tommaso <paolo.ditommaso@gmail.com>

RUN apt-get update  && apt-get install -y \
		build-essential \
		cmake \
		hdf5-tools \
		libhdf5-dev \
		hdf5-helpers \
		libhdf5-serial-dev \
		curl \
		libcurl4-gnutls-dev \
		libxml2-dev \
		libssl-dev 

		
#
# Install Sleuth 
# https://liorpachter.wordpress.com/2015/08/17/a-sleuth-for-rna-seq/
#
RUN R -e 'source("http://bioconductor.org/biocLite.R"); biocLite("rhdf5"); install.packages("devtools"); devtools::install_github("pachterlab/sleuth")'

#
# Install Kallisto
# 
RUN wget -q https://github.com/pachterlab/kallisto/archive/master.zip && \
    unzip master.zip && \
    mkdir kallisto-master/build && \
    cd kallisto-master/build && \
    cmake .. && \
	make && \
	make install
	
	
ENTRYPOINT ["/bin/sh", "-c"]	
    