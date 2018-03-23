import argparse
import logging
import numpy as np
import os
import pandas as pd
import pysam
from pomoxis.common.util import parse_regions, Region


def main():
    logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    parser = argparse.ArgumentParser('Calculate read coverage depth from a bam.')
    parser.add_argument('bam', help='.fasta/fastq file.')
    parser.add_argument('-r', '--regions', nargs='+', help='Only process given regions.')
    parser.add_argument('-p', '--prefix', help='Prefix for output, defaults to basename of bam.')
    parser.add_argument('-s', '--stride', type=int, default=1, help='Stride in genomic coordinate.')

    args = parser.parse_args()

    bam = pysam.AlignmentFile(args.bam)
    ref_lengths = dict(zip(bam.references, bam.lengths))

    if args.regions is not None:
        regions = parse_regions(args.regions, ref_lengths=ref_lengths)

    else:
        regions = [Region(ref_name=r, start=0, end=ref_lengths[r]) for r in bam.references]

    for region in regions:
        ref_len = ref_lengths[region.ref_name]
        bins = np.arange(region.start, region.end, args.stride)
        msg = 'Processing reference {}:{}-{}'
        logging.info(msg.format(region.ref_name, bins[0], bins[-1]))
        coverage = np.zeros(len(bins))
        for r_obj in bam.fetch(contig=region.ref_name, start=region.start, end=region.end):
            start_i = max((r_obj.reference_start - bins[0]) // args.stride, 0)
            end_i = min((r_obj.reference_end - bins[0]) // args.stride, len(bins))
            coverage[start_i: end_i] += 1

        # write final depth
        prefix = args.prefix
        if prefix is None:
            prefix = os.path.splitext(os.path.basename(args.bam))[0]

        depth_fp = '{}_{}_{}_{}.depth.txt'.format(prefix, region.ref_name, region.start, region.end)
        pd.DataFrame({'pos': bins,
                      'depth': coverage}).to_csv(depth_fp, sep='\t', index=False)


if __name__ == '__main__':
    main()
