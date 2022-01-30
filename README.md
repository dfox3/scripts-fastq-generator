# scripts-fastq-generator
A script to generate test fastqs

## How to run
###reqs:

- **python 3.7**
  - pybedtools
  - pyfaidx
  - biopython

###Conda option
Conda install the above requirements, 
and then run script:

``python scripts-fastq-generator/fastq_generator.py -f fasta_ref -b bed_file``

###Docker option
With docker installed, run the shell starter to load the container.

``bash start_dev_docker.sh``

``python scripts-fastq-generator/fastq_generator.py -f fasta_ref -b bed_file``

## Arguments

  -h, --help           show this help message and exit

  -o --output          Output file prefix

  -b --bed             Bed file - range input

  -f --fasta           Fasta reference for bed file - range input

  -i --inject_fasta    Fasta with sequences to inject

  -a --adapter         Adapter sequences to add to the end of fragments

  -n --num_reads       Number of read records printed

  -mn --min_frag_size  Minimum fragment size

  -mx --max_frag_size  Maximum fragment size

  -rl --read_length    Read length

  -e --error_rate      Error rate (decimal). 0 == no errors

  -ot --off_target     Number of off target reads to inject

  -q --q_score         Q score - max value 41, min value 0

  -qm --q_mode         Q mode - exact, thresh-above, or thresh-below


