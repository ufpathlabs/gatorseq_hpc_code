#DockerFile for vardict
FROM centos:7
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
RUN conda install -c conda-forge -c bioconda -c r\
         bwa==0.7.17\
         vcftools==0.1.16\
         vt==0.57721\
         r==3.2.2\
         samtools==1.9\
         htslib==1.9\
         openssl=1.0
RUN rm Miniconda3-latest-Linux-x86_64.sh
RUN git clone -b 'v1.6' --recursive https://github.com/AstraZeneca-NGS/VarDictJava.git /home/VarDictJava 
RUN cd /home/VarDictJava && ./gradlew clean installDist 
RUN cp -r /home/VarDictJava /usr/local/VarDict 
RUN conda clean --tarballs
ENV PATH="/usr/local/VarDict/build/install/VarDict/bin:$PATH"
ENV PATH="/usr/local/VarDict/VarDict:$PATH"
CMD ["bash"]

