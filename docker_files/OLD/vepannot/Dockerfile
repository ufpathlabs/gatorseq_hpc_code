#Docker file for vepAnnote 
FROM centos 
MAINTAINER harsh0280@gmail.com 

RUN yum -y install wget nano curl software-properties-common sudo git-core unzip gcc wget bzip2 \
    ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 java-1.8.0-openjdk-devel
ENV JAVA_HOME=/usr/lib/jvm/java-1.8.0-openjdk
ENV PATH=$JAVA_HOME/bin:$PATH
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/anaconda \
         export PATH="/usr/local/anaconda/bin:$PATH" \
         apt-get install -y python3-pip \
         conda config --add channels defaults \
         conda config --add channels bioconda \
         conda config --add channels conda-forge
ENV PATH="/usr/local/anaconda/bin:$PATH"
RUN conda install -c conda-forge -c bioconda -c zwdzwd\
         transvar\
         variant-effect-predictor
RUN conda install -c r r==3.2.2
RUN python -m pip install pip==9.0.3 && git clone https://github.com/counsyl/hgvs.git && pip install hgvs/ && pip install pygr==0.8.2 && pip install --upgrade pip
RUN git clone -b 'v0.21' https://github.com/weiyi-bitw/varnorm.git && pip install varnorm/              
RUN conda install -c conda-forge -c bioconda java-jdk==8.0.112
#ENV PATH="/usr/local/anaconda/bin:$PATH"
CMD [“echo”,”Image created”] 
