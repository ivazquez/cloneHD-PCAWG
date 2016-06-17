############################################################
# Dockerfile to build cloneHD workflow container
# Based on Ubuntu
############################################################

FROM ubuntu

# File Author / Maintainer
MAINTAINER Ignacio Vazquez-Garcia <ivg@sanger.ac.uk>

# Setup packages
USER root
RUN apt-get -m update && apt-get install -y make gcc gfortran \
build-essential wget git libgsl2 gsl-bin libgsl-dev libblas-dev \
liblapack-dev python-pip
RUN pip install numpy scipy PyVCF

RUN git clone https://github.com/ivazquez/cloneHD.git && cd cloneHD && git checkout pcawg
RUN cd cloneHD/src && mkdir ../build && make -f Makefile.farm

# RUN git clone git@github.com:ivazquez/cloneHD-tools.git && cd cloneHD-tools && git checkout pcawg
# RUN cd cloneHD-tools && python setup.py install

# Switch back to the ubuntu user so this tool (and the files written) are not owned by root
RUN groupadd -r -g 1000 ubuntu && useradd -r -g ubuntu -u 1000 -m ubuntu
USER ubuntu

# By default /bin/bash is executed
CMD ["/bin/bash"]