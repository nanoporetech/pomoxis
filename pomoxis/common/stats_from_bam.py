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

from pomoxis.common.util import intervaltrees_from_bed

parser = argparse.ArgumentParser(
    description="""Parse a bamfile (from a stream) and output summary stats for each read.""",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

from collections import defaultdict
insert_len_frequency = defaultdict(int)
delete_len_frequency = defaultdict(int)

parser.add_argument('bam', type=str, help='Path to bam file.')
parser.add_argument('--bed', default=None, help='.bed file of reference regions to include.')
parser.add_argument('-m', '--min_length', type=int, default=None)
parser.add_argument('-l', '--longest_indel_length', type=int, default=0,
                    help='Ignore all indels less than the set value. Default: 0 which counts all indels.')
parser.add_argument('-a', '--all_alignments', action='store_true',
                    help='Include secondary and supplementary alignments.')
parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
                    help='Output alignment stats to file instead of stdout.')
parser.add_argument('-s', '--summary', type=argparse.FileType('w'), default=sys.stderr,
                    help='Output summary to file instead of stderr.')
parser.add_argument('-t', '--threads', type=int, default=1,
                    help='Number of threads for parallel processing.')


def get_total_substitutions(read):
    """
    returns the total number of substitutions in a read
    :param read:
    :return:
    """
    tags = dict(read.tags)
    try:
        tags.get('NM')
    except:
        raise IOError("Read is missing required 'NM' tag. Try running 'samtools fillmd -S - ref.fa'.")

    counts, _ = read.get_cigar_stats()
    ins = counts[1]
    delt = counts[2]
    # NM is edit distance: NM = INS + DEL + SUB
    sub = tags['NM'] - ins - delt
    return sub


def count_from_cigartuples(cigartuples, longest_indel):
    """Count number of insertions and deletions from the cigartuples. If any insertion or deletion is longer than
    longest indel then it is ignored. If longest_indel is 0 then everything is counted.

    CIGAR OP LIST:
    M	BAM_CMATCH	0
    I	BAM_CINS	1
    D	BAM_CDEL	2
    N	BAM_CREF_SKIP	3
    S	BAM_CSOFT_CLIP	4
    H	BAM_CHARD_CLIP	5
    P	BAM_CPAD	6
    =	BAM_CEQUAL	7
    X	BAM_CDIFF	8
    B	BAM_CBACK	9

    :param cigartuples: List of cigar tuples (cigar_op, cigar_len)
    :param longest_indel: Integer number defining the highest indel to count in.
    :return: A list containing match, ins, delt counts

    """
    match, ins, delt = 0, 0, 0

    for cigar_op, cigar_len in cigartuples:
        # match
        if cigar_op == 0:
            match += cigar_len
        # insert
        if cigar_op == 1:
            insert_len_frequency[cigar_len] += 1

            if longest_indel == 0 or cigar_len <= longest_indel:
                ins += cigar_len
            # else:
            #     print('Skipping IN {}:{}'.format(cigar_op, cigar_len))
        # delt
        if cigar_op == 2:
            delete_len_frequency[cigar_len] += 1
            if longest_indel == 0 or cigar_len <= longest_indel:
                delt += cigar_len
            # else:
            #     print('Skipping DEL {}:{}'.format(cigar_op, cigar_len))

    return match, ins, delt


def stats_from_aligned_read(read, references, lengths, longest_indel_length):
    """Create summary information for an aligned read

    :param read: :class:`pysam.AlignedSegment` object
    """
    name = read.qname

    # count the number of match, insert and delete. A match can be a match or a mismatch.
    match, ins, delt = count_from_cigartuples(read.cigartuples, longest_indel_length)
    # get the number of substitute bases
    sub = get_total_substitutions(read)

    length = match + ins + delt
    total_mismatch_in_del = sub + ins + delt
    iden = 100 * float(match - sub)/match
    acc = 100 - 100 * float(total_mismatch_in_del)/length

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


def masked_stats_from_aligned_read(read, references, lengths, tree):
    """Create summary information for an aligned read over regions in bed file.

    :param read: :class:`pysam.AlignedSegment` object
    """

    tags = dict(read.tags)
    try:
        tags.get('MD')
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

def masked_stats_from_aligned_read_excluding_longindels(read, references, lengths, tree, longest_indel):
    """Create summary information for an aligned read over regions in bed file.

    :param read: :class:`pysam.AlignedSegment` object
    """

    tags = dict(read.tags)
    try:
        tags.get('MD')
    except:
        raise IOError("Read is missing required 'MD' tag. Try running 'samtools callmd - ref.fa'.")

    correct, delt, ins, sub, aligned_ref_len, masked = 0, 0, 0, 0, 0, 0
    pairs = read.get_aligned_pairs(with_seq=True)
    qseq = read.query_sequence
    pos_is_none = lambda x: (x[1] is None or x[0] is None)
    pos = None
    insertions = []

    # use a flag that indicates which operation we are in now to help count the length
    # 0 = MATCH, 1 = INSERT, 2 = DELETE, -1 = START

    running_operation_flag = -1
    running_operation_length = 0
    # TODO: refactor to use get_trimmed_pairs (as in catalogue_errors)?
    for qp, rp, rb in itertools.dropwhile(pos_is_none, pairs):
        if rp == read.reference_end or (qp == read.query_alignment_end):
            break
        pos = rp if rp is not None else pos
        if not tree.overlaps(pos) or (rp is None and not tree.overlaps(pos + 1)):
            # if rp is None, we are in an insertion, check if pos + 1 overlaps
            # (ref position of ins is arbitrary)
            print('Skipping ref {}:{}'.format(read.reference_name, pos))
            masked += 1

            # check if the previous operation is complete due to this
            if running_operation_flag != -1:
                if running_operation_flag == 1: # it was an insert
                    if longest_indel == 0 or running_operation_length <= longest_indel:
                        ins += running_operation_length
                elif running_operation_flag == 2: # it was a delete
                    if longest_indel == 0 or running_operation_length <= longest_indel:
                        delt += running_operation_length
                        aligned_ref_len += running_operation_length
                # we count the sub and mismatch separately so no need to handle that here

                # unset the running operation flag now
                running_operation_flag = -1
                running_operation_length = 0

            continue
        else:
            if rp is not None:  # we don't have an insertion
                if qp is None:  # deletion
                    if running_operation_flag == -1:
                        # the flag is unset, so set the flag to delete and start counting
                        running_operation_flag = 2
                        running_operation_length = 1

                    elif running_operation_flag == 1:
                        # the running operation is a insert but this is del, so count the insert and change the flag
                        if longest_indel == 0 or running_operation_length <= longest_indel:
                            ins += running_operation_length
                        # set the flag
                        running_operation_flag = 2
                        running_operation_length = 1

                    elif running_operation_flag == 2:
                        # the running operation is a delete so increase the length
                        running_operation_length += 1

                elif qseq[qp] == rb:  # correct
                    correct += 1
                    aligned_ref_len += 1

                    # in match and sub we don't have to maintain the flag as we count them directly
                    # just check if we need to take previous operation into account
                    if running_operation_flag == 1:
                        # it was an insert
                        if longest_indel == 0 or running_operation_length <= longest_indel:
                            ins += running_operation_length
                    elif running_operation_flag == 2:
                        # it was a delete
                        if longest_indel == 0 or running_operation_length <= longest_indel:
                            delt += running_operation_length
                            aligned_ref_len += running_operation_length

                    # unset the running operation flag now
                    running_operation_flag = -1
                    running_operation_length = 0

                elif qseq[qp] != rb:  # sub
                    sub += 1
                    aligned_ref_len += 1
                    # in match and sub we don't have to maintain the flag as we count them directly
                    # just check if we need to take previous operation into account
                    if running_operation_flag == 1:
                        # it was an insert
                        if longest_indel == 0 or running_operation_length <= longest_indel:
                            ins += running_operation_length
                    elif running_operation_flag == 2:
                        # it was a delete
                        if longest_indel == 0 or running_operation_length <= longest_indel:
                            delt += running_operation_length
                            aligned_ref_len += running_operation_length

                    # unset the running operation flag now
                    running_operation_flag = -1
                    running_operation_length = 0

            else:  # we have an insertion
                if running_operation_flag == -1:
                    # the flag is unset, so set the flag to delete and start counting
                    running_operation_flag = 1
                    running_operation_length = 1

                elif running_operation_flag == 1:
                    # the running operation is a insert so increase the running operation length
                    running_operation_length += 1

                elif running_operation_flag == 2:
                    # the running operation is a delete but this is an insert, so count the delete and change flag
                    if longest_indel == 0 or running_operation_length <= longest_indel:
                        delt += running_operation_length
                        aligned_ref_len += running_operation_length

                    running_operation_flag = 1
                    running_operation_length = 1

    # for the last operation that we recorded
    if running_operation_flag != -1:
        if running_operation_flag == 1:
            # it was an insert
            if longest_indel == 0 or running_operation_length <= longest_indel:
                ins += running_operation_length
        elif running_operation_flag == 2:
            # it was a delete
            if longest_indel == 0 or running_operation_length <= longest_indel:
                delt += running_operation_length
                aligned_ref_len += running_operation_length

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


def masked_stats_from_aligned_read_excluding_longindels_2(read, references, lengths, tree, longest_indel):
    """Create summary information for an aligned read over regions in bed file.

    :param read: :class:`pysam.AlignedSegment` object
    """

    tags = dict(read.tags)
    try:
        tags.get('MD')
    except:
        raise IOError("Read is missing required 'MD' tag. Try running 'samtools callmd - ref.fa'.")

    correct, delt, ins, sub, aligned_ref_len, masked = 0, 0, 0, 0, 0, 0
    cigartuples = read.cigartuples
    reference_position = read.reference_start
    read_index = 0
    reference_index = 0
    reference_sequence = read.get_reference_sequence()
    read_sequence = read.query_sequence
    """
    M	BAM_CMATCH	0
    I	BAM_CINS	1
    D	BAM_CDEL	2
    N	BAM_CREF_SKIP	3
    S	BAM_CSOFT_CLIP	4
    H	BAM_CHARD_CLIP	5
    P	BAM_CPAD	6
    =	BAM_CEQUAL	7
    X	BAM_CDIFF	8
    B	BAM_CBACK	9
    """
    for cigar_op, cigar_len in cigartuples:
        # match or mismatch
        if cigar_op == 7 or cigar_op == 8 or cigar_op == 0:
            for i in range(0, cigar_len):
                read_base = read_sequence[read_index]
                ref_base = reference_sequence[reference_index]

                # match
                if read_base == ref_base:
                    # check if position overlaps the bed region
                    if tree.overlaps(reference_position):
                        correct += 1
                        aligned_ref_len += 1
                    else:
                        masked += 1
                else:
                    # base mismatch
                    # check if position overlaps the bed region
                    if tree.overlaps(reference_position):
                        sub += 1
                        aligned_ref_len += 1
                    else:
                        masked += 1

                read_index += 1
                reference_position += 1
                reference_index += 1
        # insert
        elif cigar_op == 1:
            # insert
            # first check if the insert is longer than allowed longes indel length
            if longest_indel == 0 or cigar_len <= longest_indel:
                # then check if the anchor and the next reference position overlaps the range
                if tree.overlaps(reference_position-1) and tree.overlaps(reference_position):
                    ins += cigar_len
                else:
                    masked += cigar_len
            else:
                # longer than longest insert
                masked += cigar_len

            # for insert only move the read index
            read_index += cigar_len
        # delete
        elif cigar_op == 2:
            # delete
            # first check if the deletion is longer than longest allowed indel length
            if longest_indel == 0 or cigar_len <= longest_indel:
                for i in range(0, cigar_len):
                    # see if the position overlaps the bed region
                    if tree.overlaps(reference_position):
                        delt += 1
                        aligned_ref_len += 1
                    else:
                        masked += 1

                    reference_position += 1
                    reference_index += 1
            else:
                masked += cigar_len
                reference_position += cigar_len
                reference_index += cigar_len

        # ref_skip or PAD
        elif cigar_op == 3 or cigar_op == 6:
            reference_position += cigar_len
            reference_index += cigar_len
        # soft_clip
        elif cigar_op == 4:
            read_index += cigar_len
        # hard_clip
        elif cigar_op == 5:
            continue

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


def _process_reads(bam_fp, longest_indel_length, start_stop, all_alignments=False, min_length=None, bed_file=None):
    start, stop = start_stop
    counts = Counter()
    results = []
    with pysam.AlignmentFile(bam_fp, 'rb') as bam:
        if bed_file is not None:
            trees = intervaltrees_from_bed(bed_file)
        for i, read in enumerate(bam):
            if i < start:
                continue
            elif i == stop:
                break
            if read.is_unmapped:
                counts['unmapped'] += 1
                continue
            if not all_alignments and (read.is_secondary or read.is_supplementary):
                continue

            if bed_file is not None:
                if not trees[read.reference_name].overlaps(read.reference_start, read.reference_end):
                    sys.stderr.write('read {} does not overlap with any regions in bedfile\n'.format(read.query_name))
                    counts['masked'] += 1
                    continue
                else:
                    result = masked_stats_from_aligned_read(read, bam.references, bam.lengths, trees[read.reference_name])
                    result_1 = masked_stats_from_aligned_read_excluding_longindels(read, bam.references, bam.lengths,
                                                                                   trees[read.reference_name], longest_indel_length)
                    result_2 = masked_stats_from_aligned_read_excluding_longindels_2(read, bam.references, bam.lengths,
                                                                                     trees[read.reference_name], longest_indel_length)
                    assert result == result_1 == result_2
            else:
                result = stats_from_aligned_read(read, bam.references, bam.lengths, longest_indel_length)

            if min_length is not None and result['length'] < min_length:
                counts['short'] += 1
                continue

            results.append(result)

    return counts, results


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
               'read_length', 'match', 'ins', 'del', 'sub', 'iden', 'acc']

    masked_headers = ['masked']
    if args.bed is not None:
        headers.extend(masked_headers)

    args.output.write('\t'.join(headers))
    args.output.write('\n')

    # create a slice of reads to process in each thread to avoid looping through
    # bam n read times and reduce mp overhead
    with pysam.AlignmentFile(args.bam) as bam:
        n_reads = bam.count()
    n_reads_per_proc = ceil(n_reads / args.threads)
    ranges = [(start, min(start + n_reads_per_proc, n_reads))
              for start in range(0, n_reads, n_reads_per_proc)]

    func = functools.partial(_process_reads, args.bam,
                             args.longest_indel_length,
                             all_alignments=args.all_alignments,
                             min_length=args.min_length, bed_file=args.bed)

    counts = Counter()
    for batch_counts, results in mapper(func, ranges):
        counts.update(batch_counts)
        counts['total'] += len(results)
        for result in results:
            out_row = (str(result[x]) for x in headers)
            args.output.write('\t'.join(out_row))
            args.output.write('\n')
    if pool is not None:
        pool.shutdown(wait=True)

    if counts['total'] == 0:
        raise ValueError('No alignments processed. Check your bam and filtering options.')

    args.summary.write('Mapped/Unmapped/Short/Masked: {total}/{unmapped}/{short}/{masked}\n'.format_map(counts))

    # for in_len, count in sorted(insert_len_frequency.items()):
    #     print("INSERT LENGTH: ", in_len, " COUNT: ", count)
    # print("############################")
    # for del_len, count in sorted(delete_len_frequency.items()):
    #     print("DELETE LENGTH: ", del_len, " COUNT: ", count)


if __name__ == '__main__':
    main()
