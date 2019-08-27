
FROM centos:7.6.1810

MAINTAINER Srikar Chamala <srikarchamala@gmail.com>
LABEL version="1.0"

RUN yum -y install wget nano curl software-properties-common sudo git-core unzip gcc wget bzip2 \
    ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 which

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh
RUN sh Miniconda3-4.6.14-Linux-x86_64.sh -b -p /usr/local/anaconda
ENV PATH="/usr/local/anaconda/bin:$PATH"

RUN rm Miniconda3-4.6.14-Linux-x86_64.sh 


RUN conda install -c defaults -c bioconda -c conda-forge \
        cnvkit==0.9.6=py_2



CMD [“bash”]
