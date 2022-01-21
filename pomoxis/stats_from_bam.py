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

from pomoxis.util import intervaltrees_from_bed

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


def stats_from_aligned_read(read, references, lengths):
    """Create summary information for an aligned read

    :param read: :class:`pysam.AlignedSegment` object
    """
    try:
        NM = read.get_tag('NM')
    except:
        raise IOError("Read is missing required 'NM' tag. Try running 'samtools fillmd -S - ref.fa'.")

    name = read.qname
    start_offset = 0
    if read.is_secondary or read.is_supplementary:
        first_cig = read.cigartuples[0]
        if first_cig[0] == 5:  # should always be true for minimap2
            start_offset = first_cig[1]
    counts, _ = read.get_cigar_stats()
    match = counts[0] + counts[7] + counts[8]  # total of M, =, and X
    ins = counts[1]
    delt = counts[2]

    lra_flag = False
    if read.has_tag('NX'):
        # likely from lra
        # NM is number of matches, see https://github.com/ChaissonLab/LRA/issues/32
        sub = counts[8]
        lra_flag = True
    else:
        # likely from minimap2
        # NM is edit distance: NM = INS + DEL + SUB
        sub = NM - ins - delt

    length = match + ins + delt
    iden = 100 * float(match - sub) / match
    acc = 100 * float(match - sub) / length

    read_length = read.infer_read_length()
    coverage = 100 * float(read.query_alignment_length) / read_length
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
        "qstart": read.query_alignment_start + start_offset,
        "qend": read.query_alignment_end + start_offset,
        "rstart": read.reference_start,
        "rend": read.reference_end,
        "ref": references[read.reference_id],
        "aligned_ref_len": read.reference_length,
        "ref_coverage": 100*float(read.reference_length) / lengths[read.reference_id],
        "mapq": read.mapping_quality,
    }

    return results, lra_flag


def masked_stats_from_aligned_read(read, references, lengths, tree):
    """Create summary information for an aligned read over regions in bed file.

    :param read: :class:`pysam.AlignedSegment` object
    """
    try:
        MD = read.get_tag('MD')
    except:
        raise IOError("Read is missing required 'MD' tag. Try running 'samtools callmd - ref.fa'.")

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

    if match == 0:
        # no matches within bed regions - all bed ref positions were deleted.
        # skip this alignment.
        return None

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
        "mapq": read.mapping_quality,
        "masked": masked,
    }
    return results


def _process_reads(bam_fp, start_stop, all_alignments=False, min_length=None, bed_file=None):
    start, stop = start_stop
    counts = Counter()
    results = []
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
                if not trees[read.reference_name].overlaps(read.reference_start, read.reference_end):
                    sys.stderr.write('read {} does not overlap with any regions in bedfile\n'.format(read.query_name))
                    counts['masked'] += 1
                    continue
                else:
                    result = masked_stats_from_aligned_read(read, bam.references, bam.lengths, trees[read.reference_name])
                    if result is None:  # no matches within bed regions
                        counts['all_matches_masked'] +=  1
                        continue
            else:
                result, lra_flag = stats_from_aligned_read(read, bam.references, bam.lengths)

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
