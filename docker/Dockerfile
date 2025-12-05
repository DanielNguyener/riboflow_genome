FROM continuumio/miniconda3:latest

WORKDIR /riboflow

RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    && rm -rf /var/lib/apt/lists/*

RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority flexible && \
    conda config --set remote_connect_timeout_secs 30.0 && \
    conda config --set remote_read_timeout_secs 60.0 && \
    conda config --set remote_max_retries 10

COPY environment.yaml /riboflow/
COPY ribo_bigwig_environment.yaml /riboflow/

COPY modified_rfcommands /riboflow/modified_rfcommands

RUN conda env create -f environment.yaml || \
    (sleep 10 && conda env create -f environment.yaml) || \
    (sleep 20 && conda env create -f environment.yaml) && \
    conda clean -afy

RUN conda env create -f ribo_bigwig_environment.yaml || \
    (sleep 10 && conda env create -f ribo_bigwig_environment.yaml) || \
    (sleep 20 && conda env create -f ribo_bigwig_environment.yaml) && \
    conda clean -afy

RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

COPY . /riboflow/

RUN /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate ribo_genome && \
    pip install -e ./modified_rfcommands && \
    conda deactivate"

ENV PATH="/opt/conda/envs/ribo_genome/bin:/opt/conda/envs/ribo_bigwig/bin:$PATH"
ENV CONDA_DEFAULT_ENV=base
ENV PIP_ROOT_USER_ACTION=ignore

RUN echo 'eval "$(/opt/conda/bin/conda shell.bash hook)"' >> /etc/profile.d/conda.sh

CMD ["/bin/bash"]
