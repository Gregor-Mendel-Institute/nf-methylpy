#############################################################
# Dockerfile tools for methylpy
#############################################################
FROM continuumio/anaconda3
MAINTAINER Rahul Pisupati <rahul.pisupati@gmi.oeaw.ac.at>
LABEL authors="rahul.pisupati@gmi.oeaw.ac.at" \
    description="Docker image containing all requirements for the nf-core/methylpy pipeline"
USER root
RUN apt-get -y -m update && apt-get install -y wget unzip zip libgsl-dev libgsl2 bzip2 build-essential git zlib1g-dev lib32z1-dev
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nfcore-methylseq/bin:$PATH
RUN git clone https://github.com/yupenghe/methylpy.git
RUN g++ -O3 -l gsl -l gslcblas -o /opt/conda/envs/nfcore-methylseq/bin/run_rms_tests.out rms.cpp
