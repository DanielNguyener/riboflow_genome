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
    conda config --set remote_connect_timeout_secs 120.0 && \
    conda config --set remote_read_timeout_secs 300.0 && \
    conda config --set remote_max_retries 30 && \
    conda config --set ssl_verify true && \
    conda config --set always_yes true

COPY environment.yaml /riboflow/
COPY ribo_bigwig_environment.yaml /riboflow/

COPY modified_rfcommands /riboflow/modified_rfcommands

RUN /bin/bash -c '\
    SUCCESS=0; \
    for i in {1..10}; do \
        echo "========================================="; \
        echo "Attempt $i/10 to create conda environment from environment.yaml"; \
        echo "========================================="; \
        conda clean -afy || true; \
        if conda env create -f environment.yaml 2>&1; then \
            echo "Successfully created environment on attempt $i"; \
            SUCCESS=1; \
            break; \
        else \
            echo "Attempt $i failed, cleaning cache and waiting before retry..."; \
            conda clean -afy || true; \
            sleep $((i * 30)); \
        fi; \
    done; \
    if [ "$SUCCESS" -eq 0 ]; then \
        echo "ERROR: Failed to create environment after 10 attempts"; \
        echo "This is likely due to network connectivity issues."; \
        exit 1; \
    fi; \
    conda clean -afy'

RUN /bin/bash -c '\
    SUCCESS=0; \
    for i in {1..10}; do \
        echo "========================================="; \
        echo "Attempt $i/10 to create conda environment from ribo_bigwig_environment.yaml"; \
        echo "========================================="; \
        conda clean -afy || true; \
        if conda env create -f ribo_bigwig_environment.yaml 2>&1; then \
            echo "Successfully created environment on attempt $i"; \
            SUCCESS=1; \
            break; \
        else \
            echo "Attempt $i failed, cleaning cache and waiting before retry..."; \
            conda clean -afy || true; \
            sleep $((i * 30)); \
        fi; \
    done; \
    if [ "$SUCCESS" -eq 0 ]; then \
        echo "ERROR: Failed to create environment after 10 attempts"; \
        echo "This is likely due to network connectivity issues."; \
        exit 1; \
    fi; \
    conda clean -afy'

RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

COPY . /riboflow/

RUN /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate ribo_genome && \
    pip install ./modified_rfcommands && \
    conda deactivate"

ENV PATH="/opt/conda/envs/ribo_genome/bin:/opt/conda/envs/ribo_bigwig/bin:$PATH"
ENV CONDA_DEFAULT_ENV=base
ENV PIP_ROOT_USER_ACTION=ignore

RUN echo 'eval "$(/opt/conda/bin/conda shell.bash hook)"' >> /etc/profile.d/conda.sh

CMD ["/bin/bash"]
