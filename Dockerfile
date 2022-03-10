############################################################
# Dockerfile to build cloneHD workflow container
# Based on Ubuntu
############################################################

# docker build -t clonehd-pcawg:latest -t clonehd-pcawg:v0.1 .
# docker tag clonehd-pcawg:v0.1 ivazquez/clonehd-pcawg:v0.1
# docker tag clonehd-pcawg:latest ivazquez/clonehd-pcawg:latest
# docker push ivazquez/clonehd-pcawg:v0.1
# docker push ivazquez/clonehd-pcawg:latest
# docker run ivazquez/clonehd-pcawg /bin/bash -c "/opt/cloneHD_workflow.sh -v"

FROM ubuntu

# File Author / Maintainer
MAINTAINER Ignacio Vazquez-Garcia <ivg@sanger.ac.uk>

# Install software
ENV DEBIAN_FRONTEND=noninteractive

# RUN apt-get update && apt-get install -y make gcc build-essential libgslcblas0 gsl-bin libgsl-dev git
RUN apt-get update && \
    apt-get install -y make gfortran gcc && \
    apt-get install -y build-essential && \
    apt-get install -y libgslcblas0 gsl-bin libgsl-dev && \
    apt-get install -y libboost-all-dev && \
    apt-get install -y libblas-dev liblapack-dev && \
    apt-get install -y git perl python3.6 python3-pip gzip

WORKDIR /opt

# Set python3 as default
RUN ln -s /usr/bin/python3 python

# Install python modules
RUN pip install PyVCF

# Install cloneHD
RUN git clone https://github.com/ivazquez/cloneHD.git && cd cloneHD && git checkout pcawg
RUN cd cloneHD/src && mkdir ../build && make -f Makefile.farm

# Copy scripts to `WORKDIR`
RUN git clone -b docker https://github.com/ivazquez/cloneHD-PCAWG.git
# COPY cloneHD_workflow.sh *.pl *.py ./

# ENTRYPOINT ["/opt/cloneHD_workflow.sh"]
# CMD ["/opt/cloneHD_workflow.sh"]
