import argparse
import collections
import os
import pickle
import re

import matplotlib
matplotlib.use('Agg')

import  matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
import pysam

from pomoxis import bio


class AverageScore:
    """Keep track of a simple average of inputs."""

    def __init__(self, cumulative_score=0, count=0):
        self.cumulative_score = cumulative_score
        self.count = count

    def __add__(self, score):
        if isinstance(score, AverageScore):
            self.cumulative_score += score.cumulative_score
            self.count += score.count
        else:
            self.cumulative_score += score
            self.count += 1
        return self

    @property
    def average_score(self):
        """Return the average of the accumulated samples."""
        try:
            return self.cumulative_score / self.count
        except ZeroDivisionError:
            return 0

    def __repr__(self):
        return str(self.average_score)


def get_fraction_correct(counts):
    """Calculate fraction correct vs length from counts."""
    corr_res = []
    for base, c in counts.items():
        for rlen, sub_counts in c.items():
            corr_res.append({
                'qbase': base,
                'rlen': rlen,
                'n': sum(sub_counts.values()),
                'n_corr': sub_counts[rlen],
                'fr_corr': sub_counts[rlen] / max(sum(sub_counts.values()), 1)
                })
    df = pd.DataFrame(corr_res)
    df_fr = df.pivot(index='rlen', columns='qbase', values='fr_corr').reset_index()
    df_n = df.pivot(index='rlen', columns='qbase', values='n').reset_index()
    df = pd.merge(df_fr, df_n, on='rlen', suffixes=['', '_n'])
    for b, c in [('A', 'T'), ('G', 'C')]:
        # get total count and weighted-fraction correct for AT and GC
        sum_n_col = b + c + '_n'
        df[sum_n_col] = df[b + '_n'] + df[c + '_n']
        df[b + c] = ((df[b + '_n'] * df[b]) + (df[c + '_n'] * df[c])) / df[sum_n_col]
    return df


def plot_fraction_correct(df, fname, cols=('A', 'T', 'G', 'C')):
    """Plot fraction correct vs length from counts."""
    fig, ax = plt.subplots()
    for col in cols:
        ax.plot(df['rlen'], df[col], 'o-', label=col)
    ax.set_ylabel('Fraction Correct')
    ax.set_xlabel('Homopolymer length')
    ax.legend(frameon=False)
    fig.savefig(fname)


def get_relative_lengths(counts):
    rows = []
    for q_base, c in counts.items():
        for ref_len, qc in c.items():
            data = {'ref_len': ref_len, 'q_base': q_base}
            for query_len, count in qc.items():
                rows.append(data.copy())
                rows[-1]['query_len'] = query_len
                rows[-1]['count'] = count
                rows[-1]['rel_len'] = query_len - ref_len

    return pd.DataFrame(rows)


def combine_complementary_relative_lengths(data):
    data = data.copy()
    data['is_AT'] = data['q_base'].str.contains(r'[AT]')

    at_gc = data.groupby(['is_AT', 'ref_len', 'rel_len'], as_index=False).agg({"count": "sum"})
    at_gc['q_base'] = ''
    at_gc.loc[np.where(at_gc['is_AT'])[0], 'q_base'] = 'AT'
    at_gc.loc[np.where(np.logical_not(at_gc['is_AT']))[0], 'q_base'] = 'GC'
    del at_gc['is_AT']
    return at_gc


def plot_relative_lengths(data, fname):
    # some panels will be empty for longer HPs so keep track of what lengths and bases are in each row/column
    rows = sorted(data['ref_len'].unique())
    cols = sorted(data['q_base'].unique())
    fig, axes = plt.subplots(
        ncols=len(cols), nrows=len(rows), sharex=True,
        figsize=(4 * len(cols), 2 * len(rows)))
    for rl, rl_df in data.groupby(['ref_len']):
        i = rows.index(rl)
        for qb, qb_df in rl_df.groupby('q_base'):
            j = cols.index(qb)
            ax = axes[i][j]
            ax.bar(qb_df['rel_len'], qb_df['count'])
            ax.set_title('{}{}'.format(qb, rl))
            ax.xaxis.set_tick_params(which='both', labelbottom=True)
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    for ax in axes[-1,:]:
        ax.set_xlabel('Query length relative to reference')
    for ax in axes[:, 0]:
        ax.set_ylabel('Counts')
    fig.tight_layout()
    fig.savefig(fname)


def count_homopolymers(args):
    os.mkdir(args.output_dir)

    score = AverageScore()
    headers = [
        'query_name', 'ref_name', 'ref_start', 'rev_comp', 'query_start',
        'ref_base', 'query_base', 'ref_len', 'query_len']

    prefix = os.path.join(args.output_dir, 'hp')
    with open(prefix + '_catalogue.txt', 'w') as txt_fh:
        txt_fh.write('\t'.join(headers) + "\n")
        score, counts = process_bam(args.bam, txt_fh, args.homo_len, score)
    print("Found {} homopolymers in {}. Average score: {}".format(
        score.count, args.bam, score))
    save_counts(counts, prefix + '_counts.pkl')
    analyse_counts(counts, prefix)


def save_counts(counts, fp):
    # convert counts from nested defaultdict with factory function to plain dict
    counts = {k: dict(v) for k,v in counts.items()}
    with open(fp, 'wb') as pkl_fh:
        pickle.dump(counts, pkl_fh)


def load_counts(pkl):
    with open(pkl, 'rb') as fh:
        counts = pickle.load(fh)
    return counts


def merge_counts(pkl_files):
    # load and merge existing counts
    counts = collections.defaultdict(counts_factory)
    for pkl in pkl_files:
        c = load_counts(pkl)
        for base, base_counts in c.items():
            for rlen, rlen_counts in base_counts.items():
                counts[base][rlen].update(rlen_counts)
    return counts


def analyse_homopolymers(args):
    os.mkdir(args.output_dir)
    prefix = os.path.join(args.output_dir, 'hp')
    counts = merge_counts(args.pkl)
    save_counts(counts, prefix + '_counts.pkl')
    analyse_counts(counts, prefix)


def analyse_counts(counts, prefix):
    fr_correct = get_fraction_correct(counts)
    fr_correct.to_csv(prefix + '_correct_vs_len.txt', sep='\t', index=False)
    plot_fraction_correct(fr_correct, prefix + '_correct_vs_len.png')
    plot_fraction_correct(fr_correct, prefix + '_correct_vs_len_comp_pairs.png', cols=('AT', 'GC'))

    data = get_relative_lengths(counts)
    data.to_csv(prefix + '_rel_len_counts.txt', sep='\t', index=False)
    plot_relative_lengths(data, prefix + '_rel_len_counts.png')

    data_comp = combine_complementary_relative_lengths(data)
    data_comp.to_csv(prefix + '_rel_len_counts_comp_pairs.txt', sep='\t', index=False)
    plot_relative_lengths(data_comp, prefix + '_rel_len_counts_comp_pairs.png')


def get_next_aligned_base(align, seq_start):
    next_base = seq_start
    try:
        while align[next_base][0] is None:
            next_base += 1
    except IndexError:
        # We don't have any more alignments in the read!
        next_base = max(align, key=lambda x: x[0])
    return align[next_base][0]


def find_homopolymers(seq, min_length, alphabet='ACGT'):
    """
    :param seq: the sequence to search
    :param min_length: minimum hp length to find
    :param alphabet: alphabet of bases to look for

    :returns:
    """
    # backreference to force separation between different letter patterns
    # e.g. 'GTAAAAACCCCCTG' -> 'AAAAA' 'CCCCC' not 'AAAAACCCCC'
    homopolymer_expr = r"([%s])\1{%s,}" % (alphabet, min_length - 1)
    return re.finditer(homopolymer_expr, seq)


def get_longest_homopolymer(seq, base):
    try:
        length, offset = max(
            [(len(m.group()), m.start())
                for m in find_homopolymers(seq, 1, alphabet=base)],
            key=lambda x: x[0])
    except ValueError:
        # No matches found
        length, offset = (0, 0)
    return length, offset


def get_query_coords(align_pairs, seq_start, hp_len, base):
    # Find the bases around the seq_start which align to the reference and match the base
    alignment_region = align_pairs[seq_start - 1: seq_start + hp_len + 1]
    alignment_matches = [x[0] for x in alignment_region if x[0] is not None and x[2] == base]
    qstart = min(alignment_matches)
    qend = max(alignment_matches)
    return qstart, qend


def get_query_region(query_sequence, query_start, query_end, base):
    extend_fwd = 1
    extend_backward = -1
    while (query_end + extend_fwd + 1) < len(query_sequence) and \
            query_sequence[query_end + extend_fwd] == base:
        extend_fwd += 1
    while query_sequence[query_start + extend_backward] == base and \
            (query_start + extend_backward) > 0:
        extend_backward -= 1
    return query_sequence[query_start + extend_backward + 1:query_end + extend_fwd]


def counts_factory():
    return collections.defaultdict(collections.Counter)


def process_bam(input_file, out_fh, homo_len, score):
    # counts[query_base][ref_len][query_len]
    counts = collections.defaultdict(counts_factory)

    with pysam.AlignmentFile(input_file, 'r') as bam:
        for seq in bam:
            if seq.is_secondary or seq.is_supplementary or seq.is_unmapped:
                continue

            refseq = seq.get_reference_sequence().upper()
            align_pairs = seq.get_aligned_pairs(with_seq=True)
            # lookup is a mapping of the refseq indices to the query sequence indicies
            lookup = np.array(
                [i for i, v in enumerate(align_pairs) if v[1] is not None])

            # find all the homopolymers in the reference and then attempt to match
            # each one in the query sequence
            for homo_ref in find_homopolymers(refseq, min_length=homo_len):
                base = homo_ref.group()[0]
                ref_hp_len = len(homo_ref.group())
                ref_start = homo_ref.start()
                align_start = lookup[ref_start]

                try:
                    query_start, query_end = get_query_coords(
                        align_pairs, align_start, ref_hp_len, base)
                    query_region = get_query_region(
                        seq.query_sequence, query_start, query_end, base)

                    # We find the longest hp in the region to avoid insertions
                    # breaking our results
                    # e.g. ACAAAAA -> (5mer, offset=2) not (1mer offset=0)
                    # nor 7mer!
                    call_len, offset = get_longest_homopolymer(query_region, base)
                    query_start += offset
                except ValueError:
                    # If we don't find any matching base then we return the
                    # next query position that aligns
                    query_start = get_next_aligned_base(align_pairs, align_start)
                    call_len = 0
                orginal_base = base
                # Output the results
                if seq.is_reverse:
                    query_start = seq.query_length - query_start - call_len
                    base = bio.comp[base]

                score += min((call_len - ref_hp_len) ** 2, 1000)
                out_str = '\t'.join([str(i) for i in (
                    seq.query_name,
                    seq.reference_name,
                    seq.reference_start + ref_start,
                    '-' if seq.is_reverse else '+',
                    query_start,
                    orginal_base,
                    base,
                    ref_hp_len,
                    call_len)])
                out_fh.write(out_str + '\n')
                # update counts of calls
                counts[base][ref_hp_len][call_len] += 1
    return score, counts


def main():
    """Entry point for homopolymer accuracy counting."""
    parser = argparse.ArgumentParser(
        prog='homopolymer',
        description='Analyse homopolymer query and reference lengths.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(
        title='subcommands', description='valid commands',
        help='additional help', dest='command')
    subparsers.required = True

    cparser = subparsers.add_parser('count',
        help='Count homopolymers starting from a bam. ',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    cparser.set_defaults(func=count_homopolymers)
    cparser.add_argument('bam', help='Input bam file.')
    cparser.add_argument('-o', '--output_dir', default='homopolymers',
        help="Output directory (will be created).")
    cparser.add_argument('-l', '--homo_len', default=3, type=int,
        help='Minimum homopolymer length, default 3')

    aparser = subparsers.add_parser('analyse',
        help='Analyse existing counts, optionally merging multiple counters.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    aparser.set_defaults(func=analyse_homopolymers)
    aparser.add_argument('pkl', nargs='+', help='Input .pkl file(s).')
    aparser.add_argument('-o', '--output_dir', default='homopolymers',
        help="Output directory (will be created).")

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
