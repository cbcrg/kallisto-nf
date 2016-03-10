FROM debian:jessie
MAINTAINER Paolo Di Tommaso <paolo.ditommaso@gmail.com>

RUN apt-get update \ 
	&& apt-get install -y --no-install-recommends \
		ed \
		less \
		locales \
		vim-tiny \
		wget \
		ca-certificates \
	&& rm -rf /var/lib/apt/lists/*

RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8
	
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8


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
# install R 
#
RUN echo "deb http://cran.rstudio.com/bin/linux/debian jessie-cran3/" >>  /etc/apt/sources.list &&\
 apt-key adv --keyserver keys.gnupg.net --recv-key 381BA480 &&\
 apt-get update --fix-missing && \
 apt-get -y install r-base
		
#
# Install Sleuth 
# https://liorpachter.wordpress.com/2015/08/17/a-sleuth-for-rna-seq/
#
RUN R -e 'source("http://bioconductor.org/biocLite.R"); library(BiocInstaller); biocLite(c("XML","biomaRt")); biocLite("rhdf5"); install.packages("devtools", repos="http://cloud.r-project.org/"); devtools::install_github("pachterlab/sleuth")'

#
# Install Kallisto
# 
RUN wget -q https://github.com/pachterlab/kallisto/archive/v0.42.4.zip && \
    unzip v0.42.4.zip && \
    mkdir kallisto-0.42.4/build && \
    cd kallisto-*/build && \
    cmake .. && \
	make && \
	make install
  