# take this image to create a non-cuda environment
#FROM bioconductor/bioconductor_docker:RELEASE_3_16

FROM rocker/r-ver

RUN TZ=Europe/Madrid && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone


RUN apt-get update && apt-get install -y \
    software-properties-common aptitude make libncurses-dev libbz2-dev liblzma-dev gcc zlib1g-dev wget man-db

#creating user for autamate login into shell
RUN echo "auth requisite pam_deny.so" >> /etc/pam.d/su &&\
    sed 's/required/sufficient/g' -i /etc/pam.d/chsh &&\
    useradd -ms /usr/bin/bash -p "psswd" user

RUN yes | unminimize
    
ENV HOME=/home/user

WORKDIR $HOME/Downloads
RUN wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2
RUN tar xvjf samtools-1.17.tar.bz2 && cd samtools-1.17/ && ./configure && make && make install
ENV PATH="$PATH:~HOME/Downloads/samtools-1.17"


RUN apt-get update && aptitude -y install rna-star
    
RUN aptitude -y install nano

# removing cache files
RUN rm -rf /var/cache/apt/* /var/lib/apt/lists/*;
# removing temporary folders
RUN rm -rf /tmp/downloaded_packages && \
    rm -rf /tmp/RtmpIj1cu8/downloaded_packages && \
    rm -rf /tmp/RtmpmxTsSH/downloaded_packages

WORKDIR $HOME
    
USER user


CMD ["/bin/bash"]


