FROM centos:centos7
MAINTAINER Sudipta Basak <basaks@gmail.com>
ENV HOME /root
WORKDIR /root

ADD setup_perl.sh $HOME/setup_perl.sh
RUN sed -i 's/sudo //g' setup_perl.sh
ADD setup_python.sh $HOME/setup_python.sh
RUN sed -i 's/sudo //g' setup_python.sh

RUN chmod +x setup_perl.sh
RUN chmod +x setup_python.sh

RUN mkdir -p $HOME/passive-seismic/sc3/
ADD perl_envs.txt $HOME/passive-seismic/sc3/perl_envs.txt

# everything happens here
RUN ./setup_perl.sh
RUN ./setup_python.sh
RUN rm -rf $HOME/passive-seismic
