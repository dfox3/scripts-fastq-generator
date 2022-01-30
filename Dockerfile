FROM ubuntu:latest

LABEL maintainer="Dylan Fox"

ARG CONDA_PYTHON_VERSION=3.7
ARG CONDA_DIR=/opt/conda
ARG USERNAME=docker
ARG USERID=1000

# Instal basic utilities
RUN apt-get update && \
    apt-get install -y --no-install-recommends git wget unzip bzip2 sudo build-essential ca-certificates && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install miniconda
ENV PATH $CONDA_DIR/bin:$PATH
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    echo 'export PATH=$CONDA_DIR/bin:$PATH' > /etc/profile.d/conda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p $CONDA_DIR && \
    rm -rf /tmp/*

# Create the user
RUN useradd --create-home -s /bin/bash --no-user-group -u $USERID $USERNAME && \
    chown $USERNAME $CONDA_DIR -R && \
    adduser $USERNAME sudo && \
    echo "$USERNAME ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers

USER $USERNAME
WORKDIR /home/$USERNAME

ADD ./environment.yml .
RUN conda env create -f environment.yml -n fastq
ENV PATH /opt/conda/envs/fastq/bin:$PATH
ENV CONDA_DEFAULT_ENV fastq

# Demonstrate the environment is activated:
RUN echo "Make sure pybedtools is installed:"
RUN python -c "import pyfaidx"

ADD ./ /tool

RUN sudo chmod -R 777 /tool