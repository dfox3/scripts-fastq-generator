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
fh = logging.FileHandler(str(DIR) + "/../logs/debug_log.txt")
fh.setLevel(logging.DEBUG)
logger.addHandler(fh)

# -----------------------------------------------------------------------------------------------------------------------
# Main
def main():
    """
    Evaluate arguments and print fastqs.

    Parameters
    ----------
    N/A

    Returns
    -------
    N/A

    See Also
    --------
    get_sequences_from_bed
    load_sequences_from_fasta
    add_error
    generate_random_fastq_set
    write_gzip_string

    Examples
    --------
    N/A
    """

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
    """
    Generate random paired fastqs.

    Parameters
    ----------
    sequences: list
        List of sequences to build fastq records from.
    num_reads: int
        Number of reads desired for final paired fastqs.
    min_frag_size: int (optional)
        Minimum simulated insert fragment size.
    max_frag_size: int (optional)
        Maximum simulated insert fragment size.
    score: int (optional)
        Score threshold for phred quality score string building.
    mode: str (optional)
        Mode of thresholding for phred quality score string building.
    adapter_seq: str (optional)
        Adapter sequence to append to 3' end of read.

    Returns
    -------
    list:
        R1 fastq records.
    list:
        R2 fastq records.

    See Also
    --------
    _rand_read_from_seq
    _generate_read_from_sequence
    generate_record

    Examples
    --------
    N/A
    """
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
    """
    Evaluate arguments and print fastqs.

    Parameters
    ----------
    sequence: str
        DNA sequence. A, T, C, G only valid characters.
    desc: str
        Record title.
    score: int
        Phred quality score threshold
    mode: str (optional)
        Phred quality score threshold mode.

    Returns
    -------
    str
        4 lined fastq record

    See Also
    --------
    generate_title
    generate_q_string

    Examples
    --------
    N/A
    """
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
    """
    Make a fastq title from a string.

    Parameters
    ----------
    desc: str
        A string to turn into a title.

    Returns
    -------
    str
        A fastq style title for a fastq record.

    See Also
    --------
    N/A

    Examples
    --------
    N/A
    """
    return "@" + str(desc)


def generate_q_string(sequence_length, score, mode="thresh_above"):
    """
    Generate a phred quality score string.

    Parameters
    ----------
    sequence_length: int
        Length of the sequence used to build quality score.
    score: int
        Phred quality score threshold.
    mode: str (optional)
        Phred quality score thresholding mode.

    Returns
    -------
    str
        A fastq style phred quality score string for a fastq record.

    See Also
    --------
    TO_Q_SCORE

    Examples
    --------
    N/A
    """
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
    """
    Truncate read to a specified length and append a adapter sequence.

    Parameters
    ----------
    sequence: str
        DNA sequence. A, T, C, G are valid characters.
    read_length: int
        Read length.
    adapter_seq: str
        Adapter DNA seq. A, T, C, G are valid characters.
    r2: bool (optional)
        True if sequence is being used to make R2 fastq records.

    Returns
    -------
    str
        A simulated read sequence.

    See Also
    --------
    reverse_complement
    BASES

    Examples
    --------
    N/A
    """
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
    """
    Make a slightly random read based on a range for fragment sizes.

    Parameters
    ----------
    sequence: str
        DNA sequence string.
    min_frag_size: int
        Minimum fragment size.
    max_frag_size: int
        Maximum fragment size.

    Returns
    -------
    str
        A slightly randomized read.

    See Also
    --------
    N/A

    Examples
    --------
    N/A
    """
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
    """
    Use bed tools to get sequences from a bed file and fasta reference.

    Parameters
    ----------
    bed: str
        Bed file name.
    fasta: str
        Fasta reference name
    max_frag_size: int
        Maximum fragment size.
    num_offtarget_reads: int (optional)
        Number of random offtarget sequences to generate.

    Returns
    -------
    list
        A list of DNA sequences.

    See Also
    --------
    pyfaidx
    pybedtools

    Examples
    --------
    N/A
    """
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
    """
    Randomly add an error to sequence bases.

    Parameters
    ----------
    sequence: str
        DNA sequence.
    rate: float
        Decimal error rate.

    Returns
    -------
    str
        A sequence that may contain an error.

    See Also
    --------
    _rand_error_base

    Examples
    --------
    N/A
    """
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
    """
    Randomly add an error to sequence base.

    Parameters
    ----------
    base: str
        Single DNA base.
    random_max: int
        An error rate's denominator.

    Returns
    -------
    str
        A base that may contain an error.

    See Also
    --------
    N/A

    Examples
    --------
    N/A
    """
    if random.randint(1, random_max) == 1:
        base = random.choice(BASES)
    return base


def load_sequences_from_fasta(fasta):
    """
    Load a list of sequences from fasta file.

    Parameters
    ----------
    fasta: str
        Name of fasta file.

    Returns
    -------
    list
        List of sequences from fasta.

    See Also
    --------
    biopython

    Examples
    --------
    N/A
    """
    return [str(records.seq) for records in SeqIO.parse(fasta, "fasta")]


def reverse_complement(sequence):
    """
    Reverse and complement a DNA sequence.

    Parameters
    ----------
    sequence: str
        DNA sequence. A, T, C, G are valid characters.

    Returns
    -------
    str
        Reverse complemented sequence.

    See Also
    --------
    N/A

    Examples
    --------
    N/A
    """
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
    """
    Write a string to a gzip file.

    Parameters
    ----------
    string: str
        A string.
    outfile: str
        Name of outfile.

    Returns
    -------
    N/A

    See Also
    --------
    N/A

    Examples
    --------
    N/A
    """
    with gzip.open(out_file, "wb") as f:
        f.write(string.encode())


# -----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
