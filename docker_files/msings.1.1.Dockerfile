#This is a Image for msings
FROM centos:7.6.1810

MAINTAINER Srikar Chamala <srikarchamala@gmail.com>
LABEL version="1.1"

RUN yum -y install wget nano curl software-properties-common sudo git-core make unzip gcc wget bzip2 \
    ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 which


RUN wget https://www.python.org/ftp/python/2.7.16/Python-2.7.16.tgz
RUN tar xzf Python-2.7.16.tgz
RUN cd Python-2.7.16 && ./configure && make altinstall


RUN yum install -y epel-release 
RUN yum install -y python34 
RUN yum install -y python34-devel python34-setuptools zlib-devel libz-dev ncurses-devel ncurses
RUN easy_install-3.4 pip && pip3 install virtualenv
RUN mkdir msings && cd msings && git clone -b 'v3.6' https://bitbucket.org/uwlabmed/msings.git && cd msings && bash dev/bootstrap.sh && source msings-env/bin/activate

ENV PATH="/msings/msings/msings-env/bin:$PATH"
ENV PATH="/msings/msings:$PATH"

CMD [“echo”,”Image created”] 
