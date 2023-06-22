# take this image to create a non-cuda environment
FROM rocker/r-ver:4.3.1

# Install packages dependencies
RUN apt-get update && apt-get install -y \
      bzip2 \
      curl \
      make \
      gcc \
      libncurses-dev \
      libbz2-dev \
      liblzma-dev \
      zlib1g-dev \
      libcurl4-openssl-dev \
      rna-star \
      pigz \
    && rm -rf /var/cache/apt/* /var/lib/apt/lists/*;

# Install SAMTOOLS
RUN curl -kL 'https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2' | tar -C /tmp -jxf - \
    && cd /tmp/samtools-*/ \
    && ./configure \
    && make -j6 \
    && make install \
    && rm -rf /tmp/samtools-*

# Install R packages
#COPY install.R /tmp
#RUN R -q -f "/tmp/install.R" \
#    && rm -rf /tmp/Rtmp*

VOLUME /data
VOLUME /ref
WORKDIR /data
ADD pipeline/ /pipeline

ENTRYPOINT ["make","-f","/pipeline/Makefile"]



