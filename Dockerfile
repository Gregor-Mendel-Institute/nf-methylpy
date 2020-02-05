#############################################################
# Dockerfile tools for nf-methylpy
#############################################################
FROM continuumio/anaconda3:latest
MAINTAINER Rahul Pisupati <rahul.pisupati@gmi.oeaw.ac.at>
LABEL authors="rahul.pisupati@gmi.oeaw.ac.at" \
    description="Docker image containing all requirements for the nf-core/methylpy pipeline"
RUN apt-get -y update --fix-missing && apt-get -y install build-essential apt-utils && apt-get -y install zlib1g-dev && apt-get -y install libbz2-dev libcurl4-gnutls-dev libssl-dev
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
# install dependencies for DMRfind
RUN git clone https://github.com/yupenghe/methylpy.git
RUN apt-get -y install libgsl-dev libgslcblas0
RUN ln -s /usr/lib/x86_64-linux-gnu/libgsl.so.23 lib/libgsl.so.0
RUN g++ -O3 -o /opt/conda/envs/methylpy/bin/run_rms_tests.out methylpy/methylpy/rms.cpp -l gsl -l gslcblas
ENV PATH /opt/conda/envs/methylpy/bin:$PATH
