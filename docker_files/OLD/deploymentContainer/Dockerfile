#DockerFile for vardict
FROM centos:7
MAINTAINER harsh0280@gmail.com

RUN yum -y install wget nano curl software-properties-common sudo git-core unzip gcc wget bzip2 \
    ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 java-1.8.0-openjdk-devel
ENV JAVA_HOME=/usr/lib/jvm/java-1.8.0-openjdk
ENV PATH=$JAVA_HOME/bin:$PATH

RUN cd /usr/local && mkdir nextflow && cd nextflow && wget -qO- https://get.nextflow.io | bash
ENV PATH=/usr/local/nextflow:$PATH
RUN yum install -y yum-utils device-mapper-persistent-data lvm2
RUN yum-config-manager --add-repo https://download.docker.com/linux/centos/docker-ce.repo
RUN yum install -y docker-ce docker-ce-cli containerd.io
RUN yum update -y && \
    yum groupinstall -y 'Development Tools' && \
    yum install -y \
    openssl-devel \
    libuuid-devel \
    libseccomp-devel \
    wget \
    squashfs-tools
RUN export VERSION=1.11 OS=linux ARCH=amd64 && \
    wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \
    sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz && \
    rm go$VERSION.$OS-$ARCH.tar.gz && \
    echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
    echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
    source ~/.bashrc && \
    go get -u github.com/golang/dep/cmd/dep && \
    yum update -y && \
    yum install -y epel-release && \
    yum update -y && \
    yum install -y singularity-runtime singularity
#         export PATH="/usr/local/anaconda/bin:$PATH" \
#         apt-get install -y python3-pip \
#         conda config --add channels defaults \
#         conda config --add channels bioconda \
#         conda config --add channels conda-forge
#ENV PATH="/usr/local/anaconda/bin:$PATH"
#RUN conda install -c conda-forge -c bioconda -c r\
#         bwa==0.7.15\
#         vcftools==0.1.14\
#         vt==2015.11.10\
#         r==3.2.2\
#         samtools==1.3.1\
#         htslib==1.3\
#         openssl=1.0
#RUN rm Miniconda3-latest-Linux-x86_64.sh
#RUN git clone -b 'v1.4.6' --recursive https://github.com/AstraZeneca-NGS/VarDictJava.git /home/VarDictJava 
#RUN cd /home/VarDictJava && ./gradlew clean installApp 
#RUN cp -r /home/VarDictJava /usr/local/VarDict 
#RUN conda clean --tarballs
#ENV PATH="/usr/local/VarDict/build/install/VarDict/bin:$PATH"
#ENV PATH="/usr/local/VarDict/VarDict:$PATH"
CMD [“echo”,”Image created”]

