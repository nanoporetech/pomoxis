import argparse
import logging

import pysam

from pomoxis.util import filter_args, filter_read, parse_regions, Region


def main():
    logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    parser = argparse.ArgumentParser(
        prog='filter_bam',
        description='Filter a bam',
        parents=[filter_args()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('bam',
        help='input bam file.')
    parser.add_argument('-o', '--output_bam', default='filtered.bam',
        help='Output bam file.')
    parser.add_argument('-r', '--region',
        help='Only process given region.')

    args = parser.parse_args()
    logger = logging.getLogger('filter')

    with pysam.AlignmentFile(args.bam, check_sq=False) as bam:

        if args.region is not None:
            ref_lengths = dict(zip(bam.references, bam.lengths))
            region = parse_regions([args.region], ref_lengths=ref_lengths)[0]
            reads = bam.fetch(region.ref_name, region.start, region.end)
        else:
            reads = bam
        out_bam = pysam.AlignmentFile(args.output_bam, "wb", header=bam.header)
        for read in reads:
            if filter_read(read, args, logger):
                continue
            out_bam.write(read)
    out_bam.close()
    pysam.index(args.output_bam)
