import argparse
import logging
import numpy as np
import pandas as pd
import pysam
from Bio import SeqIO
from pomoxis.common.util import parse_regions, Region


def checkpoint(ref, bins, coverage, reads_kept, low_cov_sites, prefix):
    depth_fp = '{}_{}.depth.txt'.format(prefix, ref)
    pd.DataFrame({'pos': bins, 'depth': coverage}).to_csv(depth_fp, sep='\t', index=False)
    reads_fp = '{}_{}.reads.txt'.format(prefix, ref)
    pd.DataFrame({'reads': list(reads_kept)}).to_csv(reads_fp, sep='\t', index=False)
    low_cov_fp = '{}_{}.low_cov.txt'.format(prefix, ref)
    pd.DataFrame({'pos': list(low_cov_sites.keys()),
                  'depth': list(low_cov_sites.values())}).to_csv(low_cov_fp, sep='\t', index=False)


def main():
    logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    parser = argparse.ArgumentParser('subsample bam and fastx to create fastx with uniform depth')
    parser.add_argument('bam', help='input bam file.')
    parser.add_argument('fastx', help='input .fasta/fastq file of reads.')
    parser.add_argument('depth', type=int, help='Target depth.')
    parser.add_argument('-i', '--ifmt', default='fasta', choices=['fasta', 'fastq'], help='Input format.')
    parser.add_argument('-o', '--output_prefix', default='sub_sampled', help='Output prefix')
    parser.add_argument('-d', '--damping', default=0.5, type=float,
                        help='Fraction of reads required to achieve target depth to add at each iteration. ' +
                        'Setting 1 will achieve target depth in one iteration, but will likely result in ' +
                        'excess depth in places, and thus less uniform depth.')
    parser.add_argument('-c', '--chkpnt', default=100, type=int, help='Frequency at which to write depths and reads.')
    parser.add_argument('-r', '--regions', nargs='+', help='Only process given regions.')
    parser.add_argument('-s', '--stride', type=int, default=1, help='Stride in genomic coordinates.')
    parser.add_argument('-D', '--direction', choices=['fwd', 'rev'], help='Sample only forward or reverse reads.')

    args = parser.parse_args()

    bam = pysam.AlignmentFile(args.bam)
    ref_lengths = dict(zip(bam.references, bam.lengths))

    if args.regions is not None:
        regions = parse_regions(args.regions, ref_lengths=ref_lengths)

    else:
        regions = [Region(ref_name=r, start=0, end=ref_lengths[r]) for r in bam.references]

    _read_filter_ = {'fwd': lambda r: not r.is_reverse,
                     'rev': lambda r: r.is_reverse,
                      None: True,
                     }

    for region in regions:
        reads_kept = set()
        bins = np.arange(region.start, region.end, args.stride)
        msg = 'Processing region {}:{}-{}'
        logging.info(msg.format(region.ref_name, bins[0], bins[-1]))
        coverage = np.zeros(len(bins))
        count = 0
        low_cov_sites = {}
        prefix = '{}_{}X'.format(args.output_prefix, args.depth)
        if args.direction is not None:
            prefix = '{}_{}'.format(prefix, args.direction)
        while True:
            bin_i = np.argmin(coverage)
            if coverage[bin_i] >= args.depth:
                break
            count += 1
            if count % args.chkpnt == 0:
                logging.info('Min depth {} at {} (Target depth {})'.format(coverage[bin_i], bins[bin_i], args.depth))
                checkpoint(region.ref_name, bins, coverage, reads_kept, low_cov_sites, prefix)
            pos = bins[bin_i]
            reads = [r for r in bam.fetch(contig=region.ref_name, start=pos, end=pos+1)
                     if _read_filter_[args.direction](r)]
            reads_set = set((r.query_name for r in reads))
            reads_in_common = reads_kept.intersection(reads_set)
            reads_not_used = reads_set.difference(reads_in_common)
            if len(reads_not_used) == 0:
                low_cov_sites[bins[bin_i]] = coverage[bin_i]
                coverage[bin_i] = np.inf
                logging.info('Insufficient depth ({}X) at {}:{}'.format(low_cov_sites[bins[bin_i]], region.ref_name, pos))
                continue
            n_reads_to_add = max(int(round(args.damping * (args.depth - coverage[bin_i]), 0)), 1)
            logging.info('Adding {} reads at {}:{} (depth {})'.format(n_reads_to_add, region.ref_name, pos, coverage[bin_i]))
            n_reads_to_add = min(n_reads_to_add, len(reads_not_used))
            reads_to_add = np.random.choice(list(reads_not_used), n_reads_to_add, replace=False)
            reads_kept.update(reads_to_add)
            # update coverage
            r_objs = [r for r in reads if r.query_name in reads_to_add]
            for r_obj in r_objs:
                start_i = max((r_obj.reference_start - bins[0]) // args.stride, 0)
                end_i = min((r_obj.reference_end - bins[0]) // args.stride, len(bins))
                coverage[start_i: end_i] += 1

        # write final depth and reads
        checkpoint(region.ref_name, bins, coverage, reads_kept, low_cov_sites, prefix)
        fastx_ndx = SeqIO.index(args.fastx, args.ifmt)
        output = '{}_{}.{}'.format(prefix, region.ref_name, args.ifmt)
        logging.info('Writing {} sequences to {}'.format(len(reads_kept), output))
        with open(output, 'w') as fh:
            seqs = (fastx_ndx[k] for k in reads_kept)
            SeqIO.write(seqs, fh,  args.ifmt)


if __name__ == '__main__':
    main()
