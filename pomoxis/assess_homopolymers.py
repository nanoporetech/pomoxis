import argparse
import collections
import concurrent.futures
from functools import partial
from math import ceil
import os
import pickle
import re
import shutil

import matplotlib
matplotlib.use('Agg', force=True)

import  matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
import pysam

from pomoxis import bio
import pomoxis.util
from pomoxis.summary_from_stats import qscore

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
        if base != 'length':
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
    df_n[np.isnan(df_n)] = 0
    df = pd.merge(df_fr, df_n, on='rlen', suffixes=['', '_n'])
    for b in ['A', 'T', 'G', 'C']:
        if b not in df:
            df.loc[:, b] = np.NaN
            df.loc[:, b + '_n'] = 0
    for b, c in [('A', 'T'), ('G', 'C')]:
        # get total count and weighted-fraction correct for AT and GC
        sum_n_col = b + c + '_n'
        df[sum_n_col] = df[b + '_n'].add(df[c + '_n'], fill_value=0)
        df[b + c] = sum([df[a + '_n'].mul(df[a], fill_value=0) for a in [b, c]]) / df[sum_n_col]
    df['ATGC_n'] = df.loc[:, ('A_n', 'T_n', 'G_n', 'C_n')].sum(axis=1)
    df['ATGC'] = sum([df[b + '_n'].mul(df[b], fill_value=0) for b in ('AT', 'GC')]) / df['ATGC_n']
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
    plt.close(fig)


def get_relative_lengths(counts):
    rows = []
    for q_base, c in counts.items():
        if q_base != 'length':
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
        figsize=(4 * len(cols), 2 * len(rows)), squeeze=False)
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
    plt.close(fig)


def get_heatmap(data, data_comp):
    """Transform relative length data into a 2d array of called length vs true length."""
    df = data.pivot(index='query_len', columns=['q_base', 'ref_len'],
        values='count').sort_index(axis=0).sort_index(axis=1)
    column_reindex = pd.MultiIndex.from_product((['A','T','G','C'],
            range(df.columns.levels[1].min(), df.columns.levels[1].max() + 1)))
    df = df.reindex(range(df.index.min(), df.index.max() + 1))
    df[np.isnan(df)] = 0
    df = df.reindex(column_reindex, axis=1)

    df_comp = data_comp.copy()
    df_comp['query_len'] = df_comp['ref_len'] + df_comp['rel_len']
    df_comp = df_comp.pivot(index='query_len', columns=['q_base', 'ref_len'],
        values='count').sort_index(axis=0).sort_index(axis=1)
    column_reindex = pd.MultiIndex.from_product((df_comp.columns.levels[0],
            range(df_comp.columns.levels[1].min(), df_comp.columns.levels[1].max() + 1)),
            names=['q_base', 'ref_len'])
    df_comp = df_comp.reindex(range(df_comp.index.min(), df_comp.index.max() + 1))
    df_comp[np.isnan(df_comp)] = 0
    df_comp = df_comp.reindex(column_reindex, axis=1)

    df_all = df_comp.sum(level='ref_len', axis=1)
    df_all.columns = pd.MultiIndex.from_product((['ATGC',], df_all.columns), 
        names=['q_base', 'ref_len'])

    df = pd.concat((df, df_comp, df_all), axis=1)
    df[np.isnan(df)] = 0
    return df


def plot_heatmap(df, fname, cols=('A', 'T', 'G', 'C')):
    """Plot heatmaps of called length vs true length."""
    fig, axes = plt.subplots(1, len(cols), figsize=(6 * len(cols), 4),
            squeeze=False)
    for col, ax in zip(cols, axes.flatten()):
        if col in df.columns.levels[0]:
            extent = [df[col].columns.min() - 0.5, df[col].columns.max() + 0.5,
                      df[col].index.min() - 0.5, df[col].index.max() + 0.5]
            im = ax.imshow(df[col] / df[col].sum(), origin='lower',
                           extent=extent, vmin=0, vmax=1)
        ax.set_ylabel('Called Length')
        ax.set_xlabel('True Length')
        plt.colorbar(im, ax=ax, label='probability')
        ax.title.set_text(col)
    fig.tight_layout()
    fig.savefig(fname)
    plt.close(fig)


def get_errors_by_length(data, aligned_len):
    """Calculate HP errors and Q scores separately for each length."""
    e = data.copy()
    e['total_err'] = e['count'] * np.abs(e['rel_len'])
    e = e[['ref_len', 'count', 'total_err']].groupby('ref_len').sum().reset_index()
    e['total_err_rate'] = e['total_err'] / aligned_len
    e['total_err_to'] = np.cumsum(e['total_err'])
    e['total_err_from'] = np.cumsum(e['total_err'].iloc[::-1])
    e['total_err_rate_to'] = e['total_err_to'] / aligned_len
    e['total_err_rate_from'] = e['total_err_from'] / aligned_len
    e['QHP'] = qscore(e['total_err_rate'])
    e['QHP_to'] = qscore(e['total_err_rate_to'])
    e['QHP_from'] = qscore(e['total_err_rate_from'])
    e['QHPmax'] = qscore(1. / aligned_len)
    return e


def plot_errors_by_length(e, fname):
    """Plot Q score associated with HP errors as a function on HP length."""
    fig, ax = plt.subplots()
    e.plot('ref_len', ['QHP', 'QHP_to', 'QHP_from', 'QHPmax'],
            ax=ax, style=['-o', '-o', '-o', ':k'],
            label=['QHP(n)', 'QHP(<=n)', 'QHP(>=n)', 'acc. lim.'],
            ylabel="Q score",
            xlabel="Homopolymer length")
    fig.savefig(fname)
    plt.close(fig)


def count_homopolymers(args):
    os.mkdir(args.output_dir)

    # create a slice of reads to process in each thread to avoid looping through
    # bam n read times and reduce mp overhead
    ranges = [(0, 0, float('inf'))]
    if args.threads > 1:
        with pysam.AlignmentFile(args.bam) as bam:
            n_reads = bam.count(until_eof=True)
        n_reads_per_proc = ceil(n_reads / args.threads)
        ranges = [(thread, thread * n_reads_per_proc, min((thread + 1) * n_reads_per_proc, n_reads))
                   for thread in range(args.threads)]

    headers = ['query_name', 'ref_name', 'ref_start', 'rev_comp', 'query_start',
               'ref_base', 'query_base', 'ref_len', 'query_len']

    prefix = os.path.join(args.output_dir, 'hp')

    total_score = AverageScore()
    total_counts = collections.defaultdict(counts_factory)
    total_aligned_length = 0

    if args.bed is not None:
        filter_trees = pomoxis.util.intervaltrees_from_bed(args.bed)
    else:
        filter_trees = None

    f = partial(process_bam, args.bam, prefix, args.homo_len, filter_trees=filter_trees)

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as ex:
        for returned in ex.map(f, ranges):
            if returned is None:
                continue
            else:
                score, counts, aligned_length = returned
                for base in counts:
                    for rlen in counts[base]:
                        total_counts[base][rlen].update(counts[base][rlen])
                total_score += score
                total_aligned_length += aligned_length

    print("Found {} homopolymers in {}. Average score: {}".format(
        total_score.count, args.bam, total_score))
    total_counts['length'] = total_aligned_length

    with open(prefix + '_catalogue.txt', 'w') as txt_fh:
        txt_fh.write('\t'.join(headers) + "\n")
        for f in [prefix + '_catalogue_{}.txt'.format(n) for n in range(args.threads)]:
            try:
                with open(f, 'r') as infile:
                    shutil.copyfileobj(infile, txt_fh)
            except:
                print("Error merging file {}".format(f))
            else:
                print("Merged file {}, deleting".format(f))
                os.remove(f)

    save_counts(total_counts, prefix + '_counts.pkl')
    analyse_counts(total_counts, prefix, total_aligned_length)


def save_counts(counts, fp):
    # convert counts from nested defaultdict with factory function to plain dict
    counts = {k: v if type(v) == int else dict(v) for k,v in counts.items()}
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
            if base == 'length':
                if 'length' not in counts:
                    counts['length'] = base_counts
                else:
                    counts['length'] += base_counts
            else:
                for rlen, rlen_counts in base_counts.items():
                    counts[base][rlen].update(rlen_counts)
    return counts


def analyse_homopolymers(args):
    os.mkdir(args.output_dir)
    prefix = os.path.join(args.output_dir, 'hp')
    counts = merge_counts(args.pkl)
    save_counts(counts, prefix + '_counts.pkl')
    if 'length' in counts:
        aligned_length = counts['length']
    else:
        aligned_length = None
    analyse_counts(counts, prefix, aligned_length)


def analyse_counts(counts, prefix, length=None):
    fr_correct = get_fraction_correct(counts)
    fr_correct.to_csv(prefix + '_correct_vs_len.txt', sep='\t', index=False)
    plot_fraction_correct(fr_correct, prefix + '_correct_vs_len.png')
    plot_fraction_correct(fr_correct, prefix + '_correct_vs_len_comp_pairs.png', cols=('AT', 'GC'))
    plot_fraction_correct(fr_correct, prefix + '_correct_vs_len_all.png', cols=('ATGC',))

    data = get_relative_lengths(counts)
    data.to_csv(prefix + '_rel_len_counts.txt', sep='\t', index=False)
    plot_relative_lengths(data, prefix + '_rel_len_counts.png')

    data_comp = combine_complementary_relative_lengths(data)
    data_comp.to_csv(prefix + '_rel_len_counts_comp_pairs.txt', sep='\t', index=False)
    plot_relative_lengths(data_comp, prefix + '_rel_len_counts_comp_pairs.png')

    heatmap = get_heatmap(data, data_comp)
    heatmap.loc[:, ['A', 'T', 'G', 'C']].to_csv(prefix + '_heatmap.txt', sep='\t')
    heatmap.loc[:, ['AT', 'GC']].to_csv(prefix + '_heatmap_comp_pairs.txt', sep='\t')
    heatmap.loc[:, 'ATGC'].to_csv(prefix + '_heatmap_all.txt', sep='\t')
    plot_heatmap(heatmap, prefix + '_call_distribution.png')
    plot_heatmap(heatmap, prefix + '_call_distribution_comp_pairs.png', cols=('AT', 'GC'))
    plot_heatmap(heatmap, prefix + '_call_distribution_all.png', cols=('ATGC',))

    if length is not None:
        err_by_length = get_errors_by_length(data, length)
        err_by_length.to_csv(prefix + '_qhp_vs_len.txt', sep='\t', index=False)
        plot_errors_by_length(err_by_length, prefix + '_qhp_vs_len.png')


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


def process_bam(bam, prefix, homo_len, read_range, filter_trees=None):
    """Count HPs in a chunk of reads.

    :param bam: str, path to bam file
    :param prefix: str, prefix for output paths
    :param homo_len: int, minimum ref length of HPs to count
    :param read_range: (int, int, int), chunk number, read index start and end
    :param filter_trees: intervalTree, regions to analyse
    :returns: score, counts, aligned_ref_len
        score, AverageScore object
        counts, defaultdict of defaultdicts of Counters, counts by contig, ref length,
            query length
        aligned_ref_len, int
    """
    counts = collections.defaultdict(counts_factory)
    score = AverageScore()
    aligned_ref_len = 0

    catalogue_path = prefix + '_catalogue_{}.txt'.format(read_range[0])
    out_fh = open(catalogue_path, 'w')

    with pysam.AlignmentFile(bam, 'r') as bam_obj:
        aligned_ref_len = 0
        for i, seq in enumerate(bam_obj.fetch(until_eof=True)):
            if i < read_range[1]:
                continue
            elif i == read_range[2]:
                break

            if seq.is_secondary or seq.is_supplementary or seq.is_unmapped:
                continue
            seq_reference_start = seq.reference_start  # avoid recomputing for each HP
            seq_reference_end = seq.reference_end  # avoid recomputing for each HP

            if (filter_trees is not None and not
                filter_trees[seq.reference_name].overlaps(seq_reference_start, seq_reference_end)):
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
                    # If we are at the end of the alignment, skip this hp
                    if align_start == 0 or align_start + ref_hp_len == len(align_pairs):
                        continue
                    # If we don't find any matching base then we return the
                    # next query position that aligns
                    query_start = get_next_aligned_base(align_pairs, align_start)
                    call_len = 0
                orginal_base = base
                # Output the results
                if seq.is_reverse:
                    query_start = seq.query_length - query_start - call_len
                    base = bio.comp[base]

                ref_hp_pos = seq_reference_start + ref_start

                # skip hps at the end of a segment since we don't know their true length
                if ref_hp_pos == seq_reference_start or ref_hp_pos + ref_hp_len == seq_reference_end:
                    continue

                if filter_trees is not None:
                    # skip hps that at the end of an interval or spanning multiple intervals
                    # since we don't know their true length
                    interval_left_end = filter_trees[seq.reference_name][ref_hp_pos - 1]
                    interval_right_end = filter_trees[seq.reference_name][ref_hp_pos + ref_hp_len]
                    if (len(interval_left_end) == 0 or len(interval_right_end) == 0 or
                        interval_left_end != interval_right_end):
                        continue

                aligned_ref_len += ref_hp_len

                score += min((call_len - ref_hp_len) ** 2, 1000)
                out_str = '\t'.join([str(i) for i in (
                    seq.query_name,
                    seq.reference_name,
                    ref_hp_pos,
                    '-' if seq.is_reverse else '+',
                    query_start,
                    orginal_base,
                    base,
                    ref_hp_len,
                    call_len)])
                out_fh.write(out_str + '\n')
                # update counts of calls
                counts[base][ref_hp_len][call_len] += 1

    out_fh.close()

    return score, counts, aligned_ref_len


def main():
    """Entry point for homopolymer accuracy counting."""
    parser = argparse.ArgumentParser(
        prog='assess_homopolymers',
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
    cparser.add_argument('-t', '--threads', type=int, default=1,
        help='Number of threads for parallel execution.')
    cparser.add_argument('-l', '--homo_len', default=3, type=int,
        help='Minimum homopolymer length, default 3')
    cparser.add_argument('-b', '--bed',
        help='Bed file to limit search.')

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
