#DockerFile for GenerateMatrix
FROM centos:7.6.1810

MAINTAINER Srikar Chamala <srikarchamala@gmail.com>
LABEL version="1.3"


RUN yum -y install wget nano curl software-properties-common sudo git-core unzip bc gcc wget bzip2 \
    ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 java-1.8.0-openjdk-devel
ENV JAVA_HOME=/usr/lib/jvm/java-1.8.0-openjdk
ENV PATH=$JAVA_HOME/bin:$PATH

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh
RUN sh Miniconda3-4.6.14-Linux-x86_64.sh -b -p /usr/local/anaconda
ENV PATH="/usr/local/anaconda/bin:$PATH"


RUN conda install -c conda-forge -c bioconda -c r\
         rtg-tools==3.10.1\
         sambamba==0.7.0\
         bamutil==1.0.14\
         datamash==1.4.0\
         bedtools==2.26.0\
         bamtools==2.5.1

#RUN conda install -c r r==3.2.2
RUN rm Miniconda3-4.6.14-Linux-x86_64.sh 

RUN conda install -c conda-forge -c bioconda java-jdk==8.0.112

RUN conda clean --tarballs
RUN yum -y install which
CMD [“echo”,”Image created”]

