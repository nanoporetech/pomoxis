"""
Tabulate some simple alignment stats from sam.
"""

import argparse
import functools
import itertools
import sys

import intervaltree
import pybedtools
import pysam

from pomoxis.common.util import intervaltree_from_bed

parser = argparse.ArgumentParser(
        description="""Parse a bamfile (from a stream) and output summary stats for each read.""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# TODO (aw): Ideally we would use argparse.FileType to allow parsing from stream or file, but this will
#     only be supported in an upcoming release of pysam (https://github.com/pysam-developers/pysam/pull/137)
parser.add_argument('--bam', type=argparse.FileType('r'), default=sys.stdin, nargs='+')
parser.add_argument('--bed', default=None, help='.bed file of reference regions to include.')
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
    match = counts[0] # alignment match (can be a sequence match or mismatch)
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
        "match": match,
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
        "aligned_ref_len": read.reference_length,
        "ref_coverage": 100*float(read.reference_length) / lengths[read.reference_id],
    }

    return results


def masked_stats_from_aligned_read(read, references, lengths, bed_file):
    """Create summary information for an aligned read over regions in bed file.

    :param read: :class:`pysam.AlignedSegment` object
    """

    tags = dict(read.tags)
    try:
        tags.get('MD')
    except:
        raise IOError("Read is missing required 'MD' tag. Try running 'samtools callmd - ref.fa'.")

    tree = intervaltree_from_bed(bed_file, read.reference_name)
    if not tree.overlaps(read.reference_start, read.reference_end):
        sys.stderr.write('read {} does not overlap with any regions in bedfile\n'.format(read.query_name))
        return None

    correct, delt, ins, sub, aligned_ref_len, masked = 0, 0, 0, 0, 0, 0
    pairs = read.get_aligned_pairs(with_seq=True)
    qseq = read.query_sequence
    pos_is_none = lambda x: (x[1] is None or x[0] is None)
    pos = None
    insertions = []
    # TODO: refactor to use get_trimmed_pairs (as in catalogue_errors)?
    for qp, rp, rb in itertools.dropwhile(pos_is_none, pairs):
        if rp == read.reference_end or (qp == read.query_alignment_end):
            break
        pos = rp if rp is not None else pos
        if not tree.overlaps(pos) or (rp is None and not tree.overlaps(pos + 1)):
            # if rp is None, we are in an insertion, check if pos + 1 overlaps
            # (ref position of ins is arbitrary)
            # print('Skipping ref {}:{}'.format(read.reference_name, pos))
            masked += 1
            continue
        else:
            if rp is not None:  # we don't have an insertion
                aligned_ref_len += 1
                if qp is None:  # deletion
                    delt += 1
                elif qseq[qp] == rb:  # correct
                    correct += 1
                elif qseq[qp] != rb:  # sub
                    sub += 1
            else:  # we have an insertion
                ins += 1

    name = read.qname
    match = correct + sub
    length = match + ins + delt
    iden = 100 * float(match - sub) / match
    acc = 100 - 100 * float(sub + ins + delt) / length

    read_length = read.infer_read_length()
    masked_query_alignment_length = correct + sub + ins
    coverage = 100*float(masked_query_alignment_length) / read_length
    direction = '-' if read.is_reverse else '+'

    results = {
        "name": name,
        "coverage": coverage,
        "direction": direction,
        "length": length,
        "read_length": read_length,
        "match": match,
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
        "aligned_ref_len": aligned_ref_len,
        "ref_coverage": 100 * float(aligned_ref_len) / lengths[read.reference_id],
        "masked": masked,
    }
    return results


def main(arguments=None):
    args = parser.parse_args(arguments)

    count = 0
    n_unmapped = 0
    n_short = 0
    n_masked = 0
    headers = ['name', 'ref', 'coverage', 'ref_coverage', 'qstart', 'qend',
               'rstart', 'rend', 'aligned_ref_len', 'direction', 'length',
               'read_length', 'match', 'ins', 'del', 'sub', 'iden', 'acc']

    masked_headers = ['masked']
    if args.bed is not None:
        headers.extend(masked_headers)
        func = functools.partial(masked_stats_from_aligned_read, bed_file=args.bed)
    else:
        func = stats_from_aligned_read

    args.output.write('\t'.join(headers))
    args.output.write('\n')

    if not isinstance(args.bam, list):
        args.bam = ('-',)

    # TODO: multiprocessing over alignments (as in catalogue_errors)
    for samfile in (pysam.AlignmentFile(x) for x in args.bam):
        for read in samfile:
            if read.is_unmapped:
                n_unmapped += 1
                continue
            if not args.all_alignments and (read.is_secondary or read.is_supplementary):
                continue
            results = func(read, samfile.references, samfile.lengths)
            if results is None:
                n_masked += 1
                continue
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

    args.summary.write('Mapped/Unmapped/Short/Masked: {}/{}/{}\n'.format(count, n_unmapped, n_short, n_masked))

if __name__ == '__main__':
    main()
