
FROM centos:7.6.1810

MAINTAINER Srikar Chamala <srikarchamala@gmail.com>
LABEL version="1.3"

RUN yum -y install wget nano curl software-properties-common sudo git-core unzip gcc wget bzip2 \
    ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 which

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh
RUN sh Miniconda3-4.6.14-Linux-x86_64.sh -b -p /usr/local/anaconda
ENV PATH="/usr/local/anaconda/bin:$PATH"

RUN rm Miniconda3-4.6.14-Linux-x86_64.sh 

RUN conda install -c bioconda -c conda-forge \
        sambamba==0.7.0=h89e63da_1 \
        bwa==0.7.17=hed695b0_6

CMD [“bash”]
