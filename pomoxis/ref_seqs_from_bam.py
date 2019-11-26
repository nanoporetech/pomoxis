import argparse
import pysam
import sys
from Bio import SeqIO
from Bio.Seq import Seq


def _gen_seqs(bam):
    """Yield SeqIO records of the reference sequence each query is aligned to."""
    for read in bam:
        ref_seq = read.get_reference_sequence().upper()
        name = '{}_{}'.format(read.reference_name, read.query_name)
        yield SeqIO.SeqRecord(id=name, seq=Seq(ref_seq))


def main():
    parser = argparse.ArgumentParser(
        prog='ref_seqs_from_bam',
        description='Extract reference sequence that queries are aligned to',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('bam', help='input bam file, MD tag must be set (mini_align -m).')

    args = parser.parse_args()

    with pysam.AlignmentFile(args.bam) as bam:
        SeqIO.write(_gen_seqs(bam), sys.stdout, 'fasta')


if __name__ == '__main__':
    main()
