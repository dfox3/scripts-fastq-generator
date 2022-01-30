#! /bin/bash
docker build -t fastq -f Dockerfile . && \
docker run -i --rm --entrypoint /bin/bash -w /tool fastq