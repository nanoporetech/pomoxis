"""
Tabulate some simple alignment stats from sam.
"""

import argparse
import sys
import pysam
import numpy as np
from collections import Counter, defaultdict


parser = argparse.ArgumentParser(
        description="""Parse a bamfile (from a stream) and output summary stats for each read.""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# TODO (aw): Ideally we would use argparse.FileType to allow parsing from stream or file, but this will
#     only be supported in an upcoming release of pysam (https://github.com/pysam-developers/pysam/pull/137)
parser.add_argument('--bam', type=argparse.FileType('r'), default=sys.stdin, nargs='+')
parser.add_argument('--min_length', type=int, default=None)


def stats_from_aligned_read(read):
    """Create summary information for an aligned read

    :param read: :class:`pysam.AlignedSegment` object
    """

    tags = dict(read.tags)
    try:
        tags.get('NM')
    except:
        raise IOError("Read is missing required 'NM' tag. Try running 'samtools fillmd -S - ref.fa'.")

    name = read.qname
    if read.flag == 4:
        return None
    counts = defaultdict(int)
    for (i, j) in read.cigar:
        counts[i] += j
    match = counts[0]
    ins = counts[1]
    delt = counts[2]
    # NM is edit distance: NM = INS + DEL + SUB
    sub = tags['NM'] - ins - delt
    length = match + ins + delt
    iden = 100*float(match - sub)/match
    acc = 100 - 100*float(tags['NM'])/length

    read_length = read.infer_read_length()
    coverage = 100*float(read.query_alignment_length) / read_length
    direction = '-' if read.is_reverse else '+'

    results = {
        "name":name,
        "coverage":coverage,
        "direction":direction,
        "length":length,
        "read_length":read_length,
        "ins":ins,
        "del":delt,
        "sub":sub,
        "iden":iden,
        "acc":acc
    }

    return results


def main(arguments=None):
    args = parser.parse_args(arguments)

    qual_counter = Counter()
    count = 0
    unmapped = 0
    total_id = 0.0
    total_acc = 0.0
    headers = ['name', 'ref', 'coverage', 'ref_coverage', 'qstart', 'qend', 'rstart', 'rend', 'direction', 'length', 'read_length', 'ins', 'del', 'sub', 'iden', 'acc', 'mean_q']
    sys.stdout.write('\t'.join(headers))
    sys.stdout.write('\n')

    if not isinstance(args.bam, list):
        args.bam = ('-',)

    for samfile in (pysam.AlignmentFile(x) for x in args.bam):
        for read in samfile:
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            results = stats_from_aligned_read(read)
            if results is None:
                unmapped += 1
                continue
            if args.min_length is not None and results['length'] < args.min_length:
                continue

            results['qstart'] = read.query_alignment_start
            results['qend'] = read.query_alignment_end
            results['rstart'] = read.reference_start
            results['rend'] = read.reference_end

            results["ref"] = samfile.references[read.reference_id]
            results["ref_coverage"] = 100*float(read.reference_length) / samfile.lengths[read.reference_id]
            try:
                results["mean_q"] = np.mean(read.query_alignment_qualities)
                qual_counter.update(read.query_alignment_qualities)
            except:
                results["mean_q"] = 0

            count += 1
            out_row = (str(results[x]) for x in headers)
            sys.stdout.write('\t'.join(out_row))
            sys.stdout.write('\n')
            total_id += results["iden"]
            total_acc += results["acc"]
        samfile.close()

    sys.stderr.write('Mean identity: {}\n'.format(total_id / count))
    sys.stderr.write('Mean accuracy: {}\n'.format(total_acc / count))
    sys.stderr.write('Mapped/Unmapped: {}/{}\n'.format(count, unmapped))
    sys.stderr.write('Qualities: {}\n'.format(qual_counter))


if __name__ == '__main__':
    main()
