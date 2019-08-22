# Using OpenJDK 8
# This Dockerfile does not require any files that are in the GATK4 repo.
FROM centos 
MAINTAINER harsh0280@gmail.com

#### Basic image utilities
RUN yum update && \
    #apt-get upgrade -y && \
    yum install -y python && \
    yum install -y wget \
    curl \
    bc \
    unzip \
    less \
    bedtools \
    samtools \
    openjdk-8-jdk \
    tabix \
    software-properties-common && \
    yum -y clean  && \
    yum -y autoclean  && \
    yum -y autoremove

#### Specific for google cloud support
RUN export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/yum $CLOUD_SDK_REPO main" | tee -a /etc/yum/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/yum/doc/yum-key.gpg | yum-config-manager --add-repo - && \
    yum update -y && yum install google-cloud-sdk -y && \
    yum -y autoremove && \
    yum -y clean
###########

# Set environment variables.
ENV HOME /root

# Define working directory.
WORKDIR /root

# Define default command.
CMD ["bash"]

ENV JAVA_LIBRARY_PATH /usr/lib/jni
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/

RUN java -version

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9 && \
    add-apt-repository "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" && \
    yum update && \
    yum install -y --force-yes \
        r-base-dev=3.2.5-1xenial \
        r-base-core=3.2.5-1xenial && \
    yum -y clean && \
    yum -y autoremove && \
    yum  -y autoclean

COPY install_R_packages.R install_R_packages.R
RUN Rscript install_R_packages.R
# Deleting unneeded caches
RUN rm -rf /var/lib/yum/lists/*
