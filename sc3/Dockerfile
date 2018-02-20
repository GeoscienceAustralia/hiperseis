FROM centos:centos7
MAINTAINER Sudipta Basak <basaks@gmail.com>
ENV HOME /root
WORKDIR /root
RUN echo $HOME
RUN yum update -y
RUN yum groupinstall 'Development Tools' -y
RUN yum install -y git

RUN git clone https://github.com/GeoscienceAustralia/passive-seismic/ && \
    cd $HOME/passive-seismic && git submodule update --init --recursive --remote

RUN sed -i 's/sudo //g' $HOME/passive-seismic/sc3/create_perl_python.sh
RUN chmod +x $HOME/passive-seismic/sc3/create_perl_python.sh

# everything happens here
RUN $HOME/passive-seismic/sc3/create_perl_python.sh
