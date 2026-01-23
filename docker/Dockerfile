FROM ubuntu:22.04

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update --fix-missing && \
  apt-get install -q -y wget curl bzip2 libbz2-dev git build-essential zlib1g-dev locales vim fontconfig 

# Set the locale
RUN locale-gen en_US.UTF-8  
ENV LANG en_US.UTF-8  
ENV LANGUAGE en_US:en  
ENV LC_ALL en_US.UTF-8     

# Install conda
ARG CONDA=Miniconda3-py39_23.5.2-0-Linux-x86_64.sh
RUN curl -LO https://repo.anaconda.com/miniconda/$CONDA && \
    bash $CONDA -p /miniconda3 -b && \
    rm $CONDA 
ENV PATH /miniconda3/bin:${PATH}

# Install conda dependencies
ADD environment.yaml /
ADD VERSION /

RUN mkdir -p ~/.conda && \
    touch ~/.conda/tos && \
    echo 'terms_of_service_accepted: true' >> ~/.conda/tos && \
    export CONDA_ALWAYS_YES="true" && \
    export ANACONDA_TOS_ACCEPTED="true" && \
    conda config --set always_yes yes --set changeps1 no && \
    conda config --add channels conda-forge && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --get && \
    conda update -q conda && \
    conda install -y conda-libmamba-solver && \
    conda config --set solver libmamba && \
    conda info -a && \
    conda env create -f environment.yaml && \
    conda clean -afy

# Make the "ribo_genome" environment the default.
ENV PATH /miniconda3/envs/ribo_genome/bin:$PATH
RUN echo "source activate ribo_genome" >> ~/.bashrc
