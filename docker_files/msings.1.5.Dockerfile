#This is a Image for msings
FROM ubuntu:14.04

MAINTAINER Srikar Chamala <srikarchamala@gmail.com>
LABEL version="1.5"

#This was partly picked up from
#https://gitlab.com/sheenamt/singularity_files/blob/master/Singularity-msings-v3.5

RUN apt-get update && \
	apt-get install -y --no-install-recommends \
	make \
	git \
	ca-certificates \
	build-essential \
	software-properties-common \
	wget

RUN apt-get install --yes \
	scons g++ cmake libgsl0-dev libncurses5-dev gfortran

# set default java environment variable
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64

RUN add-apt-repository ppa:openjdk-r/ppa -y && \
    apt-get update && \
    apt-get install -y --no-install-recommends openjdk-8-jdk && \
    rm -rf /var/lib/apt/lists/*

RUN wget https://repo.anaconda.com/miniconda/Miniconda2-4.6.14-Linux-x86_64.sh
RUN sh Miniconda2-4.6.14-Linux-x86_64.sh -b -p /usr/local/anaconda/PYTHON2

#Do not want anaconda bin to be in the path as it may interfere with python 2.* version above
#ENV PATH="/usr/local/anaconda/bin:$PATH"

RUN /usr/local/anaconda/PYTHON2/bin/conda install -c anaconda  -c conda-forge -c bioconda -c r \
	samtools==0.1.18 \
	varscan==2.3.7 \
	numpy==1.11.1 \
	virtualenv==15.1.0 \
	pandas==0.19.2 
	

ENV PATH="/usr/local/anaconda/PYTHON2/bin:$PATH"


RUN cd /opt \
    && git clone  https://bitbucket.org/uwlabmed/msings.git \
    && cd msings && git checkout 0302893 && python setup.py install


ENV PATH="/opt/msings:$PATH"

CMD [“echo”,”Image created”] 
