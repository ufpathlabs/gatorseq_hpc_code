#Docker file for icall
FROM centos:7.6.1810
MAINTAINER Srikar Chamala <srikarchamala@gmail.com>
LABEL version="1.0"

RUN yum -y install wget nano curl software-properties-common sudo git-core unzip gcc wget bzip2 \
    ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 java-1.8.0-openjdk-deve which

ENV JAVA_HOME=/usr/lib/jvm/java-1.8.0-openjdk
ENV PATH=$JAVA_HOME/bin:$PATH


RUN wget https://repo.anaconda.com/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh
RUN sh Miniconda3-4.6.14-Linux-x86_64.sh -b -p /usr/local/anaconda
ENV PATH="/usr/local/anaconda/bin:$PATH"
RUN conda install -c conda-forge -c bioconda -c defaults -c r \
	r==3.2.2 \
	pandas==0.24.2\
	pysam==0.8.4\
	pyvcf==0.6.7 \
	delly==0.7.6 \
	biopython==1.66 \
	samtools==1.7 \
	bcftools==1.6 \
	openssl=1.0


#RUN conda config --add channels defaults
#RUN conda config --add channels bioconda
#RUN conda config --add channels conda-forge

#RUN conda install -c bioconda -y samtools==1.7 
#RUN conda install -c bioconda -y  bcftools==1.6 
#RUN conda install openssl=1.0


RUN pip install coloredlogs==5.2 && pip install reportlab==3.3.0
RUN cd /usr/local && git clone https://bitbucket.org/srikarchamala/icallsv.git && git clone https://bitbucket.org/srikarchamala/iannotatesv.git
ENV PATH="/usr/local/icallsv/0.0.6/bin:$PATH"
ENV PATH="/usr/local/iannotatesv/1.0.9/bin:$PATH"
ENV PATH="/usr/local/icallsv/0.0.6/lib/python2.7/site-packages/iCallSV-1.0.9-py2.7.egg/iCallSV:$PATH"
ENV PATH="/usr/local/icallsv/0.0.6/iCallSV/R/Rscripts:$PATH"
CMD [“echo”,”Image created”]

