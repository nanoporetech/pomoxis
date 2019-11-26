import argparse
import pysam
from collections import Counter
from pomoxis.summary_from_stats import qscore
from pomoxis.stats_from_bam import stats_from_aligned_read


def get_errors(aln, ref_seq=None):
    seq = aln.query_sequence
    errors = {}
    insertions = ''
    pairs = aln.get_aligned_pairs(with_seq=True)
    if pairs[0][0] is None or pairs[0][1] is None:
        raise ValueError('It does not look like bam is trimmed to a common alignment window')
    for qp, rp, rb in pairs[::-1]:  # process pairs in reverse to easily accumulate insertions
        if qp is None:  # deletion
            errors[rp] = (rb, '-')
        elif rp is None:  # insertion
            insertions += seq[qp]
        # if we reach here, qp is not None and rp is not None
        elif len(insertions) > 0:
            # this also includes cases where the ref and query don't agree
            # e.g. ref   A-TGC
            #      query GTTGC  would emit, A -> GT
            errors[rp] = (rb.upper(), seq[qp] + insertions[::-1])
            insertions = ''
        elif seq[qp] != rb:  # mismatch
            errors[rp] = (rb.upper(), seq[qp])
        if ref_seq is not None and rp is not None:
            assert ref_seq[rp] == rb.upper()

    if ref_seq is not None:
        # check we can scuff up reference using errors and get back our query
        scuffed = scuff_ref(ref_seq, errors)
        if not scuffed == aln.query_sequence:
            raise ValueError('Scuffing up reference with errors did not recreate query')
    return errors


def count_errors(errors):
    counts = Counter()
    for rp, (ref, query) in errors.items():
        if query == '-':
            counts['del'] += 1
        elif ref != query[0]:
            counts['sub'] += 1
        if len(query) > 1:
            counts['ins'] += len(query) - 1
    return counts


def scuff_ref(ref_seq, errors):
    orig_seq = list(ref_seq)
    for rp, (rb, alt) in errors.items():
        assert orig_seq[rp] == rb
        if alt == '-':
            orig_seq[rp] = ''
        else:
            orig_seq[rp] = alt
    return ''.join(orig_seq)


def get_qscores(counts, ref_len):
    return {
        'Q(acc)': qscore(sum(counts.values()) / ref_len),
        'Q(iden)': qscore(counts['sub'] / (ref_len - counts['del'])),
        'Q(ins)': qscore(counts['ins'] / ref_len),
        'Q(del)': qscore(counts['del'] / ref_len),
    }


def main():
    parser = argparse.ArgumentParser(
        prog='common_errors_from_bam',
        description='Get errors common to multiple assemblies aligned to ref.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('bam', help='input bam file containing assemblies trimmed to a common alignment window')
    parser.add_argument('ref_fasta', help='reference fasta file of the reference over that alignment window')
    parser.add_argument('-o', '--output_prefix', default='common_errors',
                        help='Prefix for outputs.')

    args = parser.parse_args()

    bam = pysam.AlignmentFile(args.bam)
    if len(bam.references) > 1:
        raise ValueError('Bam should have just one reference')
    ref_lengths = dict(zip(bam.references, bam.lengths))

    ref_seq = pysam.FastaFile(args.ref_fasta).fetch(bam.references[0])

    # reads should already be trimmed to a common aligment start and end point
    reads = [r for r in bam]
    ref_end, ref_start = reads[0].reference_end, reads[0].reference_start
    ref_len = ref_end - ref_start
    if not (all([r.reference_end == ref_end for r in reads]) and
            all([r.reference_start == ref_start for r in reads])):
        raise ValueError('Alignments have not been trimmed to a common overlap window, try trim_alignments')

    # get errors in each read
    data = {}
    qscores = []
    for aln in reads:
        errors = get_errors(aln, ref_seq)
        counts = count_errors(errors)
        # check we got the same error counts as stats_from_aligned_read
        stats = stats_from_aligned_read(aln, list(ref_lengths.keys()), list(ref_lengths.values()))
        for k in counts.keys():
            if stats[k] != counts[k]:
                msg = "Error counts {} don't match those from the CIGAR str {}."
                raise ValueError(msg.format(counts, {k: stats[k] for k in counts.keys()}))
        qscores.append((aln.query_name, get_qscores(counts, ref_len)))
        data[aln.query_name] = errors

    # get intersection of errors
    names = list(data.keys())
    common_errors = set(data[names[0]].keys())  # set of reference positions
    for name in names[1:]:
        common_errors = common_errors.intersection(set(data[name].keys()))
    remaining_errors = {}
    # loop through common errors, checking ref is the same and retaining the
    # error with the shortest edit distance
    for rp in common_errors:
        ref = data[names[0]][rp][0]
        assert all([d[rp][0] == ref for d in data.values()])  # refs should be same
        alts = [d[rp][1] for d in data.values()]
        if len(set([len(alt) for alt in alts])) > 1:
            # we should take the best one
            alts = sorted(alts, key=lambda x: len(x))
            shortest = alts[0]
            others = alts[1:]
            if shortest == '-' and any([alt[0] == ref for alt in others]):
                # the alt with the insertion contained the ref,
                # and the alt with the deletion did not contain the insertion
                # so two wrongs make a right!
                continue
            else:
                remaining_errors[rp] = (ref, shortest)
        else:  # pick most common, arbitrary for equal numbers
            remaining_errors[rp] = (ref, Counter(alts).most_common()[0][0])

    # write fasta of ref scuffed with just common errors
    ref_scuffed = scuff_ref(ref_seq, remaining_errors)
    with open('{}.fasta'.format(args.output_prefix), 'w') as fh:
        fh.write('>{}\n'.format(args.output_prefix))
        fh.write(ref_scuffed)

    remaining_counts = count_errors(remaining_errors)
    qscores.append(('common_errors', get_qscores(remaining_counts, ref_len)))

    # print qscores of individual reads and overlapping errors
    cols = ['Q(acc)', 'Q(iden)', 'Q(del)', 'Q(ins)']
    with open('{}.txt'.format(args.output_prefix), 'w') as fh:
        fh.write('\t'.join(['name'] + cols) + '\n')
        for name, d in qscores:
            fh.write('\t'.join([name] + [str(d[c]) for c in cols]) + '\n')


if __name__ == '__main__':
    main()
