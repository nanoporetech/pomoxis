"""
Tabulate some simple alignment stats from sam.
"""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import functools
import itertools
from math import ceil
import sys

import pysam

from pomoxis.util import intervaltrees_from_bed, masked_stats_from_aligned_read, stats_from_aligned_read

parser = argparse.ArgumentParser(
        prog='stats_from_bam',
        description="""Parse a bamfile (from a stream) and output summary stats for each read.""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('bam', type=str, help='Path to bam file.')
parser.add_argument('--bed', default=None, help='.bed file of reference regions to include.')
parser.add_argument('-m', '--min_length', type=int, default=None)
parser.add_argument('-a', '--all_alignments', action='store_true',
                    help='Include secondary and supplementary alignments.')
parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
                    help='Output alignment stats to file instead of stdout.')
parser.add_argument('-s', '--summary', type=argparse.FileType('w'), default=sys.stderr,
                    help='Output summary to file instead of stderr.')
parser.add_argument('-t', '--threads', type=int, default=1,
                    help='Number of threads for parallel processing.')



def _process_reads(bam_fp, start_stop, all_alignments=False, min_length=None, bed_file=None):
    start, stop = start_stop
    counts = Counter()
    results = []
    lra_flag = False
    with pysam.AlignmentFile(bam_fp, 'rb') as bam:
        if bed_file is not None:
            trees = intervaltrees_from_bed(bed_file)
        for i, read in enumerate(bam.fetch(until_eof=True)):
            if i < start:
                continue
            elif i == stop:
                break
            if read.is_unmapped:
                counts['unmapped'] += 1
                continue
            if not all_alignments and (read.is_secondary or read.is_supplementary):
                continue

            lra_flag = False
            if bed_file is not None:
                if (read.reference_name not in trees or
                    not trees[read.reference_name].has_overlap(read.reference_start, read.reference_end)):
                    sys.stderr.write('read {} does not overlap with any regions in bedfile\n'.format(read.query_name))
                    counts['masked'] += 1
                    continue
                else:
                    result = masked_stats_from_aligned_read(read, trees[read.reference_name])
                    if result is None:  # no matches within bed regions
                        counts['all_matches_masked'] +=  1
                        continue
            else:
                result, lra_flag = stats_from_aligned_read(read)

            if min_length is not None and result['length'] < min_length:
                counts['short'] += 1
                continue

            results.append(result)

    return counts, results, lra_flag


def main(arguments=None):
    args = parser.parse_args(arguments)

    if args.threads > 1 and args.bed is not None:
        pool = ProcessPoolExecutor(args.threads)
        mapper = pool.map
    else:
        pool = None
        mapper = map
        args.threads = 1

    headers = ['name', 'ref', 'coverage', 'ref_coverage', 'qstart', 'qend',
               'rstart', 'rend', 'aligned_ref_len', 'direction', 'length',
               'read_length', 'match', 'ins', 'del', 'sub', 'iden', 'acc',
               'mapq']

    masked_headers = ['masked']
    if args.bed is not None:
        headers.extend(masked_headers)

    args.output.write('\t'.join(headers))
    args.output.write('\n')

    # create a slice of reads to process in each thread to avoid looping through
    # bam n read times and reduce mp overhead
    ranges = [(0, float('inf'))]
    if args.threads > 1:
        with pysam.AlignmentFile(args.bam) as bam:
            n_reads = bam.count(until_eof=True)
        n_reads_per_proc = ceil(n_reads / args.threads)
        ranges = [(start, min(start + n_reads_per_proc, n_reads))
                   for start in range(0, n_reads, n_reads_per_proc)]

    func = functools.partial(_process_reads, args.bam,
                             all_alignments=args.all_alignments,
                             min_length=args.min_length, bed_file=args.bed)

    detected_lra = False
    counts = Counter()
    for batch_counts, results, lra_flag in mapper(func, ranges):
        counts.update(batch_counts)
        counts['total'] += len(results)
        for result in results:
            out_row = (str(result[x]) for x in headers)
            args.output.write('\t'.join(out_row))
            args.output.write('\n')
        detected_lra = detected_lra or lra_flag
    if pool is not None:
        pool.shutdown(wait=True)

    if detected_lra:
        args.summary.write('Some reads contained an NX tag, assuming these were aligned with LRA.\n')
        args.summary.write('If this is not the case, results may be wrong.\n')

    if counts['total'] == 0:
        args.summary.write('No alignments processed. Check your bam and filtering options.\n')

    args.summary.write('Mapped/Unmapped/Short/Masked/Skipped(all matches masked): {total}/{unmapped}/{short}/{masked}/{all_matches_masked}\n'.format_map(counts))


if __name__ == '__main__':
    main()
