#DockerFile for vardict
FROM centos:7.6.1810

MAINTAINER Srikar Chamala <srikarchamala@gmail.com>
LABEL version="1.3"

RUN yum -y install wget nano curl software-properties-common sudo git-core unzip gcc wget bzip2 \
    ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 java-1.8.0-openjdk-devel which
ENV JAVA_HOME=/usr/lib/jvm/java-1.8.0-openjdk
ENV PATH=$JAVA_HOME/bin:$PATH

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh
RUN sh Miniconda3-4.6.14-Linux-x86_64.sh -b -p /usr/local/anaconda
ENV PATH="/usr/local/anaconda/bin:$PATH"

#RUN bash Miniconda3-4.6.14-Linux-x86_64.sh -b -p /usr/local/anaconda \
#         conda config --add channels defaults \
#         conda config --add channels bioconda \
#         conda config --add channels conda-forge


RUN conda install -c conda-forge -c bioconda -c r \
         vcftools==0.1.16=he860b03_3 \
         perl-vcftools-vcf==0.1.16=pl526_2 \
         vt==0.57721=hdf88d34_2 \
         samtools==1.9=h8571acd_11 \
         htslib==1.9=ha228f0b_7 \
         r==3.6.0=r36_0 


         #vardict-java==1.6.0=0
         #openssl=1.0.2s=h7b6447c_0 

RUN rm Miniconda3-4.6.14-Linux-x86_64.sh 

RUN git clone -b 'v1.6' --recursive https://github.com/AstraZeneca-NGS/VarDictJava.git /home/VarDictJava 
RUN cd /home/VarDictJava && ./gradlew clean installDist 
RUN cp -r /home/VarDictJava /usr/local/VarDict 
RUN conda clean --tarballs
ENV PATH="/usr/local/VarDict/build/install/VarDict/bin:$PATH"
ENV PATH="/usr/local/VarDict/VarDict:$PATH"


CMD ["bash"]

