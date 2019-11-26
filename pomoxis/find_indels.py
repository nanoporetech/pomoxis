import argparse
from concurrent.futures import ProcessPoolExecutor
import functools
import itertools
from math import ceil
import sys

import pysam

parser = argparse.ArgumentParser(
        prog='find_indels',
        description='Parse a bamfile and document indels.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('bam', type=str, help='Path to bam file.')
parser.add_argument('-m', '--min_indel_length', type=int, default=0,
                    help='Filter out indels shorter than this length.')
parser.add_argument('-a', '--all_alignments', action='store_true',
                    help='Include secondary and supplementary alignments.')
parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
                    help='Output indels to file instead of stdout.')
parser.add_argument('-b', '--bed', type=argparse.FileType('w'), default=None,
                    help='Additionaly output a .bed file.')
parser.add_argument('-t', '--threads', type=int, default=1,
                    help='Number of threads for parallel processing.')


def _get_indels(pairs, min_len=0):
    start = 0
    def check_is_indel(p):
         return None in p

    for is_indel, group in itertools.groupby(pairs, key=check_is_indel):
        length = sum(1 for x in group)
        if is_indel and length >= min_len:
            yield start, length
        start += length


def get_trimmed_pairs(aln):
    """Trim aligned pairs to the alignment.

    :param aln: `pysam.AlignedSegment` object
    :yields pairs:
    """
    def pos_is_none(x):
        return x[1] is None or x[0] is None

    for qp, rp in itertools.dropwhile(pos_is_none, aln.get_aligned_pairs()):
        if (rp == aln.reference_end or qp == aln.query_alignment_end):
            break
        yield qp, rp


def _get_read_indels(read, min_indel_length=5):

    direction = '-' if read.is_reverse else '+'
    pairs = list(get_trimmed_pairs(read))
    indels = []
    for start_pair_i, n_pair in _get_indels(pairs, min_indel_length):
        qp, rp = pairs[start_pair_i]
        if rp is None:  # insertion
            typ = 'ins'
            qp_start = qp
            qp_end = qp + n_pair
            # insertions are left-aligned, start at previous ref pos
            rp_start = pairs[start_pair_i - 1][1]
            # bed intervals are half-open, so exclude last position
            rp_end = rp_start + 1
        else:
            typ = 'del'
            rp_start = rp
            rp_end = rp + n_pair

            qp_start = pairs[start_pair_i - 1][0]
            qp_end = qp_start + 1

        indel = {
            'type': typ,
            'ref_name': read.reference_name,
            'ref_start': rp_start,
            'ref_end': rp_end,
            'query_name': read.query_name,
            'query_start': qp_start,
            'query_end': qp_end,
            'indel_len': n_pair,
            'direction': direction,
         }
        indels.append(indel)

    return indels


def _process_reads(bam_fp, start_stop, all_alignments=False, min_indel_length=None):
    start, stop = start_stop
    results = []
    with pysam.AlignmentFile(bam_fp, 'rb') as bam:
        for i, read in enumerate(bam):
            if i < start:
                continue
            elif i == stop:
                break
            if read.is_unmapped:
                continue
            if not all_alignments and (read.is_secondary or read.is_supplementary):
                continue

            results.extend(_get_read_indels(read, min_indel_length=min_indel_length))

    return results


def main(arguments=None):
    args = parser.parse_args(arguments)

    if args.threads > 1 and args.bed is not None:
        pool = ProcessPoolExecutor(args.threads)
        mapper = pool.map
    else:
        pool = None
        mapper = map
        args.threads = 1

    headers = ['type', 'indel_len', 'ref_name', 'ref_start', 'ref_end', 'query_name',
               'query_start', 'query_end', 'direction']

    bed_cols = ['ref_name', 'ref_start', 'ref_end']

    args.output.write('\t'.join(headers))
    args.output.write('\n')

    # create a slice of reads to process in each thread to avoid looping through
    # bam n read times and reduce mp overhead
    with pysam.AlignmentFile(args.bam) as bam:
        n_reads = bam.count()
    n_reads_per_proc = ceil(n_reads / args.threads)
    ranges = [(start, min(start + n_reads_per_proc, n_reads))
               for start in range(0, n_reads, n_reads_per_proc)]

    func = functools.partial(_process_reads, args.bam,
                             all_alignments=args.all_alignments,
                             min_indel_length=args.min_indel_length)

    for results in mapper(func, ranges):
        for result in results:
            out_row = (str(result[x]) for x in headers)
            args.output.write('\t'.join(out_row) + '\n')
            if args.bed is not None:
                out_row = (str(result[x]) for x in bed_cols)
                args.bed.write('\t'.join(out_row) + '\n')

    if pool is not None:
        pool.shutdown(wait=True)


if __name__ == '__main__':
    main()
