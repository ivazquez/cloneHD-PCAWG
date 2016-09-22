############################################################
# Dockerfile to build cloneHD workflow container
# Based on Ubuntu
############################################################

FROM ubuntu

# File Author / Maintainer
MAINTAINER Ignacio Vazquez-Garcia <ivg@sanger.ac.uk>

# Install software
RUN apt-get update && apt-get install -y make gcc build-essential libgsl2 gsl-bin libgsl-dev git

WORKDIR /opt

# Install cloneHD
RUN git clone https://github.com/ivazquez/cloneHD.git && cd cloneHD && git checkout pcawg
RUN cd cloneHD/src && mkdir ../build && make -f Makefile.farm