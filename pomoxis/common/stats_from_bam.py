"""
Tabulate some simple alignment stats from sam.
"""

import argparse
import sys
import pysam


parser = argparse.ArgumentParser(
        description="""Parse a bamfile (from a stream) and output summary stats for each read.""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# TODO (aw): Ideally we would use argparse.FileType to allow parsing from stream or file, but this will
#     only be supported in an upcoming release of pysam (https://github.com/pysam-developers/pysam/pull/137)
parser.add_argument('--bam', type=argparse.FileType('r'), default=sys.stdin, nargs='+')
parser.add_argument('-m', '--min_length', type=int, default=None)
parser.add_argument('-a', '--all_alignments', action='store_true',
                    help='Include secondary and supplementary alignments.')
parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
                    help='Output alignment stats to file instead of stdout.')
parser.add_argument('-s', '--summary', type=argparse.FileType('w'), default=sys.stderr,
                    help='Output summary to file instead of stderr.')


def stats_from_aligned_read(read, references, lengths):
    """Create summary information for an aligned read

    :param read: :class:`pysam.AlignedSegment` object
    """

    tags = dict(read.tags)
    try:
        tags.get('NM')
    except:
        raise IOError("Read is missing required 'NM' tag. Try running 'samtools fillmd -S - ref.fa'.")

    name = read.qname
    counts, _ = read.get_cigar_stats()
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
        "name": name,
        "coverage": coverage,
        "direction": direction,
        "length": length,
        "read_length": read_length,
        "ins": ins,
        "del": delt,
        "sub": sub,
        "iden": iden,
        "acc": acc,
        "qstart": read.query_alignment_start,
        "qend": read.query_alignment_end,
        "rstart": read.reference_start,
        "rend": read.reference_end,
        "ref": references[read.reference_id],
        "ref_coverage": 100*float(read.reference_length) / lengths[read.reference_id],
    }

    return results


def main(arguments=None):
    args = parser.parse_args(arguments)

    count = 0
    n_unmapped = 0
    n_short = 0
    headers = ['name', 'ref', 'coverage', 'ref_coverage', 'qstart', 'qend', 'rstart', 'rend', 'direction', 'length', 'read_length', 'ins', 'del', 'sub', 'iden', 'acc']
    args.output.write('\t'.join(headers))
    args.output.write('\n')

    if not isinstance(args.bam, list):
        args.bam = ('-',)

    for samfile in (pysam.AlignmentFile(x) for x in args.bam):
        for read in samfile:
            if read.is_unmapped:
                n_unmapped += 1
                continue
            if not args.all_alignments and (read.is_secondary or read.is_supplementary):
                continue
            results = stats_from_aligned_read(read, samfile.references, samfile.lengths)
            if args.min_length is not None and results['length'] < args.min_length:
                n_short += 1
                continue

            count += 1
            out_row = (str(results[x]) for x in headers)
            args.output.write('\t'.join(out_row))
            args.output.write('\n')
        samfile.close()

    if count == 0:
        raise ValueError('No alignments processed. Check your bam and filtering options.')

    args.summary.write('Mapped/Unmapped/Short: {}/{}/{}\n'.format(count, n_unmapped, n_short))


if __name__ == '__main__':
    main()
