import argparse
import logging
import numpy as np
import os
import pandas as pd
import pysam
from pomoxis.util import parse_regions, Region


def coverage_of_region(region, bam_fp, stride):
    """Get coverage data for a region"""

    bins = np.arange(region.start, region.end, stride)
    msg = 'Processing region {}:{}-{}'
    logging.info(msg.format(region.ref_name, bins[0], bins[-1]))
    coverage_by_is_rev = {True: np.zeros(len(bins)), False: np.zeros(len(bins))}
    for r_obj in pysam.AlignmentFile(bam_fp).fetch(contig=region.ref_name, start=region.start, end=region.end):
        start_i = max((r_obj.reference_start - bins[0]) // stride, 0)
        end_i = min((r_obj.reference_end - bins[0]) // stride, len(bins))
        coverage_by_is_rev[r_obj.is_reverse][start_i: end_i] += 1
    total_depth = coverage_by_is_rev[True] + coverage_by_is_rev[False]

    df = pd.DataFrame({'pos': bins,
                       'depth': total_depth,
                       'depth_fwd': coverage_by_is_rev[False],
                       'depth_rev': coverage_by_is_rev[True],
                       }
    )
    return df


def coverage_summary_of_region(*args):
    df = coverage_of_region(*args)
    return df.describe()


def main():
    logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    parser = argparse.ArgumentParser(
        prog='coverage_from_bam',
        description='Calculate read coverage depth from a bam.',
        epilog='By default a file is written per reference sequence, this can be changed with the `--one_file` '
        'option. If overlapping regions are specified, `--one_file` should not be used.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('bam', help='.fasta/fastq file.')
    parser.add_argument('-r', '--regions', nargs='+', help='Only process given regions.')
    grp = parser.add_mutually_exclusive_group()
    grp.add_argument('-p', '--prefix', help='Prefix for output, defaults to basename of bam.')
    grp.add_argument('-o', '--one_file', help='Single output file with "region" column.')
    parser.add_argument('-s', '--stride', type=int, default=1000, help='Stride in genomic coordinate.')
    parser.add_argument('--summary_only', action='store_true',
                        help='Output only the depth_summary.txt file')

    args = parser.parse_args()

    bam = pysam.AlignmentFile(args.bam)
    ref_lengths = dict(zip(bam.references, bam.lengths))

    if args.regions is not None:
        regions = parse_regions(args.regions, ref_lengths=ref_lengths)
    else:
        regions = [Region(ref_name=r, start=0, end=ref_lengths[r]) for r in bam.references]

    summary = {}
    first = True
    if args.one_file is not None:
        args.one_file = open(args.one_file, 'w')
    for region in regions:
        prefix = args.prefix
        if prefix is None:
            prefix = os.path.splitext(os.path.basename(args.bam))[0]

        region_str = '{}_{}_{}'.format(region.ref_name, region.start, region.end)
        df = coverage_of_region(region, args.bam, args.stride)
        summary[region_str] = df['depth'].describe()

        if not args.summary_only:
            if args.one_file is not None:
                df['rname'] = region.ref_name
                df.to_csv(args.one_file, sep='\t', index=False, header=first)
            else:
                depth_fp = '{}_{}.depth.txt'.format(prefix, region_str)
                df.to_csv(depth_fp, sep='\t', index=False)
        first = False
    if args.one_file is not None:
        args.one_file.close()


    summary_fp = '{}_depth_summary.txt'.format(prefix)
    summary_df = pd.DataFrame(summary).T.reset_index().rename(columns={'index': 'region'})
    summary_df.to_csv(summary_fp, index=False, sep='\t')


if __name__ == '__main__':
    main()
