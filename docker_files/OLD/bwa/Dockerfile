#This is a sample Image 
FROM centos 
MAINTAINER harsh0280@gmail.com 

RUN yum -y install wget nano curl software-properties-common sudo git-core unzip gcc wget bzip2 \
    ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1
#RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN wget https://repo.anaconda.com/miniconda/Miniconda2-4.7.10-Linux-x86_64.sh
RUN bash Miniconda2-4.7.10-Linux-x86_64.sh -b -p /usr/local/anaconda \
         export PATH="/usr/local/anaconda/bin:$PATH" \
         apt-get install -y python3-pip \
         conda config --add channels defaults \
         conda config --add channels bioconda \
         conda config --add channels conda-forge
ENV PATH="/usr/local/anaconda/bin:$PATH"
RUN conda install -c conda-forge -c bioconda \
        bwa==0.7.17 \
        java-jdk==8.0.112 \
        sambamba==0.7.0-1
         
RUN conda install -c conda-forge -c bioconda \
#ENV PATH="/usr/local/anaconda/bin:$PATH"
CMD [“bash”] 
