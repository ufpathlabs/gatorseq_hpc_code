
FROM centos:7.6.1810

MAINTAINER Srikar Chamala <srikarchamala@gmail.com>
LABEL version="1.0"

RUN yum -y install wget nano curl software-properties-common sudo git-core unzip gcc wget bzip2 \
    ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh

RUN bash Miniconda3-4.6.14-Linux-x86_64.sh -b -p /usr/local/anaconda \
         export PATH="/usr/local/anaconda/bin:$PATH" \
         apt-get install -y python3-pip \
         conda config --add channels defaults \
         conda config --add channels bioconda \
         conda config --add channels conda-forge

ENV PATH="/usr/local/anaconda/bin:$PATH"

RUN conda install -c bioconda -c conda-forge \
        sambamba==0.7.0=h89e63da_1 \
        bwa==0.7.17=hed695b0_6

CMD [“bash”]

