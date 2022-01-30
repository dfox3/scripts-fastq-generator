#! /bin/bash
docker build -t fastq -f Dockerfile . && \
docker run -it --rm --entrypoint /bin/bash -w /tool fastq