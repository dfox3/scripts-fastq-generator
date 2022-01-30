#! /bin/bash
docker build -t fastq -f Dockerfile . && \
docker run -it -w /tool fastq