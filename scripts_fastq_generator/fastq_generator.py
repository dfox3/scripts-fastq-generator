import sys
import pathlib
import argparse
import logging
import random
import pybedtools
import math
import gzip
from Bio import SeqIO
from pyfaidx import Fasta

# -----------------------------------------------------------------------------------------------------------------------
# Arguments
parser = argparse.ArgumentParser(description="Generate test fastqs")
parser.add_argument(
    "-o", metavar="--output", type=str, required=True, help="Output file prefix"
)
parser.add_argument(
    "-b", metavar="--bed", type=str, default="empty", help="Bed file - range input"
)
parser.add_argument(
    "-f",
    metavar="--fasta",
    type=str,
    default="empty",
    help="Fasta reference for bed file - " + "range input",
)
parser.add_argument(
    "-i",
    metavar="--inject_fasta",
    default="empty",
    help="Fasta with sequences to inject",
)
parser.add_argument(
    "-a",
    metavar="--adapter",
    type=str,
    default="",
    help="Adapter sequences to add to the end of " + "fragments",
)
parser.add_argument(
    "-n",
    metavar="--num_reads",
    type=int,
    default=1000,
    help="Number of read records printed",
)
parser.add_argument(
    "-mn",
    metavar="--min_frag_size",
    type=int,
    default=150,
    help="Minimum fragment size",
)
parser.add_argument(
    "-mx",
    metavar="--max_frag_size",
    type=int,
    default=250,
    help="Maximum fragment size",
)
parser.add_argument(
    "-rl", metavar="--read_length", type=int, default=150, help="Read length"
)
parser.add_argument(
    "-e",
    metavar="--error_rate",
    type=float,
    default=0.0,
    help="Error rate (decimal). 0 == no errors",
)
parser.add_argument(
    "-ot",
    metavar="--off_target",
    type=int,
    default=0,
    help="Number of off target reads to inject",
)
parser.add_argument(
    "-q",
    metavar="--q_score",
    type=int,
    default=30,
    help="Q score - max value 41, min value 0",
)
parser.add_argument(
    "-qm",
    metavar="--q_mode",
    type=str,
    default="thresh-above",
    help="Q mode - exact, " + "thresh-above, or thresh-below",
)

# -----------------------------------------------------------------------------------------------------------------------
# Global variables
TO_Q_SCORE = [
    "!",
    '"',
    "#",
    "$",
    "%",
    "&",
    "'",
    "(",
    ")",
    "*",
    "+",
    ",",
    "-",
    ".",
    "/",
    "0",
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "7",
    "8",
    "9",
    ":",
    ";",
    "<",
    "=",
    ">",
    "?",
    "@",
    "A",
    "B",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "J",
]
BASES = ["A", "T", "C", "G"]
DIR = pathlib.Path(__file__).parent.resolve()

# -----------------------------------------------------------------------------------------------------------------------
# Log config
logger = logging.getLogger("fastq_generator")
logging.basicConfig(
    format="%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s",
    datefmt="%Y-%m-%d:%H:%M:%S",
    level=logging.DEBUG,
)
with open(str(DIR) + "/../logs/debug_log.txt", "w") as f:
    pass
fh = logging.FileHandler(str(DIR) + "/../test_out/debug_log.txt")
fh.setLevel(logging.DEBUG)
logger.addHandler(fh)

# -----------------------------------------------------------------------------------------------------------------------
# Main
def main():

    options = parser.parse_args()
    sequences = []
    if options.f != "empty" and options.b == "empty":
        logger.warning(
            "\tNo bed file added along with fasta. If you are looking to add fasta sequences "
            + "independent of bed file, please us -i instead of -f"
        )
    if options.b != "empty":
        if options.f == "empty":
            logger.error(
                "\tPlease specify fasta file reference (-f) for bed file. Exiting script"
            )
            exit(2)
        sequences = sequences + get_sequences_from_bed(
            options.b, options.f, options.mx, num_offtarget_reads=options.ot
        )
    if options.i != "empty":
        sequences = sequences + load_sequences_from_fasta(options.i)
    errored_sequences = [add_error(seq, options.e) for seq in sequences]
    r1_fastq, r2_fastq = generate_random_fastq_set(
        errored_sequences,
        options.rl,
        options.n,
        min_frag_size=options.mn,
        max_frag_size=options.mx,
        score=options.q,
        mode=options.qm,
        adapter_seq=options.a,
    )
    write_gzip_string(r1_fastq, str(options.o) + "_S1_L001_R1_001.fastq.gz")
    write_gzip_string(r2_fastq, str(options.o) + "_S1_L001_R2_001.fastq.gz")


def generate_random_fastq_set(
    sequences,
    read_length,
    num_reads,
    min_frag_size=150,
    max_frag_size=250,
    score=30,
    mode="thresh-above",
    adapter_seq="",
):
    random.shuffle(sequences)
    r1_fastq = ""
    r2_fastq = ""
    try:
        trimmed_sequences = [
            _rand_read_from_seq(
                sequences[r % len(sequences)], min_frag_size, max_frag_size
            )
            for r in range(num_reads)
        ]
        r1_fasta = [
            _generate_read_from_sequence(seq, read_length, adapter_seq)
            for seq in trimmed_sequences
        ]
        r2_fasta = [
            _generate_read_from_sequence(seq, read_length, adapter_seq, r2=True)
            for seq in trimmed_sequences
        ]
        r1_fastq = "".join(
            [
                generate_record(read, str(i), score, mode=mode)
                for i, read, in enumerate(r1_fasta)
            ]
        )
        r2_fastq = "".join(
            [
                generate_record(read, str(i), score, mode=mode)
                for i, read, in enumerate(r2_fasta)
            ]
        )
    except ZeroDivisionError:
        logger.warning(
            "\tFastqs will be empty based on input parameters. Please change parameters to add "
            + "sequences."
        )
    return r1_fastq, r2_fastq


def generate_record(sequence, desc, score, mode="thresh-above"):
    record = (
        str(generate_title(desc))
        + "\n"
        + str(sequence)
        + "\n"
        + "+\n"
        + str(generate_q_string(len(sequence), score, mode=mode))
    )
    return record


def generate_title(desc):
    return "@" + str(desc)


def generate_q_string(sequence_length, score, mode="thresh_above"):
    if score > 41:
        score = 41
    if score < 0:
        score = 0
    q_string = ""
    if mode == "exact":
        q_string = "".join([TO_Q_SCORE[score] for i in range(sequence_length)])
    if mode == "thresh-above":
        q_string = "".join(
            [TO_Q_SCORE[random.randint(score, 41)] for i in range(sequence_length)]
        )
    if mode == "thresh-below":
        q_string = "".join(
            [TO_Q_SCORE[random.randint(0, score)] for i in range(sequence_length)]
        )
    return q_string


def _generate_read_from_sequence(sequence, read_length, adapter_seq, r2=False):
    if r2:
        sequence = reverse_complement(sequence)

    base_legality = [base in BASES for base in adapter_seq]
    if not all(base_legality):
        logger.error(
            "\tSequence contains illegal bases preventing adapter addition: "
            + str(adapter_seq)
        )
        exit(2)
    adapted_sequence = sequence + adapter_seq
    if len(adapted_sequence) > read_length:
        adapted_sequence = adapted_sequence[:read_length]
    return adapted_sequence


def _rand_read_from_seq(sequence, min_frag_size, max_frag_size):
    # ignore randomizing start and stop if read is too small
    if len(sequence) > min_frag_size:
        # randomly choose start
        start = random.randint(0, len(sequence) - min_frag_size)
        sequence = sequence[start:]
        # revalidate seq length before moving on
        if len(sequence) > min_frag_size:
            if len(sequence) < max_frag_size:
                max_frag_size = len(sequence)
            # randomly choose stop
            stop = random.randint(min_frag_size, max_frag_size)
            sequence = sequence[:stop]
    return sequence


def get_sequences_from_bed(bed, fasta, max_frag_size, num_offtarget_reads=0):
    with open("fasta.genome", "w") as f:
        f.write(
            "\n".join([str(rec.name) + "\t" + str(len(rec)) for rec in Fasta(fasta)])
        )
    bed_tool = pybedtools.BedTool(bed)
    slopped_bed_tool = bed_tool.slop(b=max_frag_size, g="fasta.genome")
    fasta_tool = slopped_bed_tool.sequences(fi=fasta)
    sequences = [line.strip() for line in open(fasta_tool.seqfn) if line[0] != ">"]
    if num_offtarget_reads > 0:
        bed_tool = pybedtools.BedTool()
        random_bed_tool = bed_tool.random(
            l=max_frag_size, n=num_offtarget_reads, g="fasta.genome"
        )
        random_fasta_tool = random_bed_tool.sequences(fi=fasta)
        sequences = sequences + [
            line.strip() for line in open(random_fasta_tool.seqfn) if line[0] != ">"
        ]
    return sequences


def add_error(sequence, rate):
    errored_sequence = sequence
    if rate != 0:
        random_max = math.floor((1 / rate))
        if random_max < 1:
            random_max = 1
        errored_sequence = "".join(
            [_rand_error_base(base, random_max) for base in sequence]
        )
    return errored_sequence


def _rand_error_base(base, random_max):
    if random.randint(1, random_max) == 1:
        base = random.choice(BASES)
    return base


def load_sequences_from_fasta(fasta):
    return [str(records.seq) for records in SeqIO.parse(fasta, "fasta")]


def reverse_complement(sequence):
    base_legality = [base in BASES for base in sequence]
    if not all(base_legality):
        logger.error(
            "\tSequence contains illegal bases preventing reverse complement: "
            + str(sequence)
        )
        exit(2)
    comp = str.maketrans("ACGT", "TGCA")
    return sequence.translate(comp)[::-1]


def write_gzip_string(string, out_file):
    with gzip.open(out_file, "wb") as f:
        f.write(string.encode())


# -----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
