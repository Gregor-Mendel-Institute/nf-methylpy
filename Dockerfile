#############################################################
# Dockerfile tools for methylpy
#############################################################
FROM ubuntu:18.04
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
MAINTAINER Rahul Pisupati <rahul.pisupati@gmi.oeaw.ac.at>
LABEL authors="rahul.pisupati@gmi.oeaw.ac.at" \
    description="Docker image containing all requirements for the nf-core/methylpy pipeline"
# Get basic ubuntu packages needed
ENV PATH /opt/conda/bin:$PATH
RUN apt-get update --fix-missing && apt-get install -y wget bzip2 ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1  \
    git mercurial subversion \
    build-essential libgsl0-dev
RUN wget --quiet https://repo.anaconda.com/archive/Anaconda2-5.3.0-Linux-x86_64.sh -O ~/anaconda.sh && \
    /bin/bash ~/anaconda.sh -b -p /opt/conda && \
    rm ~/anaconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc
# RUN apt-get -y -m update && apt-get install -y libgsl-dev libgsl2  build-essential git zlib1g-dev lib32z1-dev
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nfcore-methylseq/bin:$PATH
RUN git clone https://github.com/yupenghe/methylpy.git
RUN g++ -O3 -l gsl -l gslcblas -o /opt/conda/envs/nfcore-methylseq/bin/run_rms_tests.out methylpy/methylpy/rms.cpp
