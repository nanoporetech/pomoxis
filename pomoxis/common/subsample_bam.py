import argparse
import logging

from Bio import SeqIO
from intervaltree import IntervalTree, Interval
import numpy as np
import pysam


from pomoxis.common.util import parse_regions, Region


def main():
    logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    parser = argparse.ArgumentParser('subsample bam and fastx to create fastx with uniform depth')
    parser.add_argument('bam',
        help='input bam file.')
    parser.add_argument('fastx',
        help='input .fasta/fastq file of reads.')
    parser.add_argument('depth', nargs='+', type=int,
        help='Target depth.')
    parser.add_argument('-i', '--ifmt', default='fasta', choices=['fasta', 'fastq'],
        help='Input format.')
    parser.add_argument('-o', '--output_prefix', default='sub_sampled',
        help='Output prefix')
    parser.add_argument('-r', '--regions', nargs='+',
        help='Only process given regions.')
    parser.add_argument('-s', '--stride', type=int, default=1000,
        help='Stride in genomic coordinates when searching for new reads. Smaller can lead to more compact pileup.')
    parser.add_argument('-p', '--profile', type=int, default=1000,
        help='Stride in genomic coordinates for depth profile.')
    parser.add_argument('-O', '--orientation', choices=['fwd', 'rev'],
        help='Sample only forward or reverse reads.')

    args = parser.parse_args()

    bam = pysam.AlignmentFile(args.bam)
    ref_lengths = dict(zip(bam.references, bam.lengths))
    fastx_ndx = SeqIO.index(args.fastx, args.ifmt)

    if args.regions is not None:
        regions = parse_regions(args.regions, ref_lengths=ref_lengths)

    else:
        regions = [Region(ref_name=r, start=0, end=ref_lengths[r]) for r in bam.references]

    for region in regions:
        # make interval tree of reads for querying
        tree = IntervalTree()
        for r in bam.fetch(region.ref_name, region.start, region.end):
            if (r.is_reverse and args.direction == 'fwd') or \
               (not r.is_reverse and args.direction == 'rev'):
                continue
            # trim reads to region
            tree.add(Interval(
                max(r.reference_start, region.start), min(r.reference_end, region.end),
                r.query_name))

        logging.info('Starting pileup.')
        coverage = np.zeros(region.end - region.start, dtype=np.uint16)
        reads = set()
        n_reads = 0
        iteration = 0
        while True:
            cursor = 0
            while cursor < ref_lengths[region.ref_name]:
                hits = _nearest_overlapping_point(tree, cursor)
                if hits is None:
                    cursor += args.stride
                else:
                    found = False
                    for read in hits:
                        # one would like to simply remove read from the tree
                        #   but there seems to be a bug in the removal
                        if read.data in reads:
                            continue
                        reads.add(read.data)
                        cursor = read.end
                        found = True
                        coverage[read.begin - region.start:read.end - region.start] += 1
                        break
                    if not found:
                        cursor += args.stride
            iteration += 1
            median_depth = np.median(coverage)
            stdv_depth = np.std(coverage)
            logging.info(u'Iteration {}. reads: {}, depth: {:.0f}X (\u00B1{:.1f}).'.format(
                iteration, len(reads), median_depth, stdv_depth))
            # output when we hit a target
            if median_depth in args.depth:
                logging.info("Hit target depth {}.".format(median_depth))
                _write_pileup(
                    fastx_ndx, args.output_prefix, region, reads,
                    coverage, args.ifmt, args.profile)
                if median_depth == max(args.depth):
                    break
            # exit if nothing happened this iteration
            if n_reads == len(reads):
                logging.warn("No reads added, finishing pileup.")
                break
            else:
                n_reads = len(reads)


def _nearest_overlapping_point(src, point):
    """Find the interval with the closest start point to a given point.

    :param src: IntervalTree instance.
    :param point: query point.

    :returns: Interval instance of interval with closest start.

    """
    items = src.search(point)
    if len(items) == 0:
        return None 
    items = sorted(items, key=lambda x: x.end - x.begin, reverse=True)
    items.sort(key=lambda x: abs(x.begin - point))
    return items


def _write_pileup(source, prefix, region, sequences, coverage, seq_fmt, profile):
    depth = int(np.median(coverage))
    prefix = '{}_{}X'.format(prefix, depth)

    output = '{}_{}.{}'.format(prefix, region.ref_name, seq_fmt)
    logging.info('Writing {} sequences to {}'.format(len(sequences), output))
    with open(output, 'w') as fh:
        seqs = (source[k] for k in sequences)
        SeqIO.write(seqs, fh,  seq_fmt)
    output = '{}_{}.depth'.format(prefix, region.ref_name)
    logging.info('Writing depth profile to {}.'.format(output))
    end = profile * (len(coverage) // profile)
    cov_blocks = coverage[0:end].reshape(-1, profile)
    depth_profile = np.mean(cov_blocks, axis=1, dtype=np.uint32)
    start = region.start + profile // 2
    positions = (start + profile * x for x in range(len(depth_profile)))
    with open(output, 'w') as fh:
        fh.write("position\tdepth\n")
        for pos, depth in zip(positions, depth_profile):
            fh.write("{}\t{}\n".format(pos, depth))


if __name__ == '__main__':
    main()
