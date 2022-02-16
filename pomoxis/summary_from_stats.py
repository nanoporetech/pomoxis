"""
Tabulate some simple alignment stats from sam.
"""

import argparse
import sys
import numpy as np
import pandas as pd
import warnings
from collections import OrderedDict
import re

parser = argparse.ArgumentParser(
        description="""Summarise output of `stats_from_bam`.""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stderr,
                    help='Output summary to file instead of stderr.')
parser.add_argument('-p', '--percentiles', type=int, nargs='+', default=(10,50,90),
                    help='Percentiles for summary.')
parser.add_argument('-pr', '--per_reference', action='store_true',
                    help='Also output a summary for each reference.')
parser.add_argument('-a', '--accumulate', type=int, default=1,
                    help='Accumulate chunks in the summary.')

def qscore(d):
    with warnings.catch_warnings():
        # some values might be zero
        warnings.simplefilter("ignore")
        q = -10 * np.log10(d)
    return q

def summarise_stats(d, percentiles=(10, 50, 90), accumulate=1):

    def sum_every_n(x):
        return np.add.reduceat(np.array(x), np.arange(0, len(x), accumulate))

    s = []
    def summarise_stat(name, f):
        # per alignment errors from which to calculate metrics, bounds etc.
        r = OrderedDict()
        r['name'] = name
        r['mean'] = f[0](d)
        vals = f[1](d)
        for p, pval in zip(percentiles, np.percentile(vals, percentiles)):
            r['q{}'.format(p)] = pval
        return r

    # create functions to get overall length-weighted averages
    # and per-alignment errors from which to calculate percentiles
    f = OrderedDict()

    if accumulate > 1:
        per_align_func = sum_every_n
    else:
        per_align_func = lambda x: x

    f['err_ont'] = (lambda d: np.sum(d['sub'] + d['ins'] + d['del']) / np.sum(d['length']), # weighted
                    lambda d: per_align_func(d['sub'] + d['ins'] + d['del']) / per_align_func(d['length'])) # per n alignments
    f['err_bal'] = (lambda d: np.sum(d['sub'] + d['ins'] + d['del']) / np.sum(d['aligned_ref_len']),
                    lambda d: per_align_func(d['sub'] + d['ins'] + d['del']) / per_align_func(d['aligned_ref_len']))
    f['iden'] = (lambda d: np.sum(d['sub']) / (np.sum(d['rend'] - d['rstart'] - d['del'])),
                 lambda d: per_align_func(d['sub']) / per_align_func(d['rend'] - d['rstart'] - d['del']))
    f['del'] = (lambda d: np.sum(d['del']) / (np.sum(d['rend'] - d['rstart'])),
               lambda d: per_align_func(d['del']) / per_align_func(d['rend'] - d['rstart']))
    f['ins'] = (lambda d: np.sum(d['ins']) / (np.sum(d['rend'] - d['rstart'])),
               lambda d: per_align_func(d['ins']) / per_align_func(d['rend'] - d['rstart']))

    for name, f in f.items():
        s.append(summarise_stat(name, f))

    return s

def main(arguments=None):
    args = parser.parse_args(arguments)
    n=args.accumulate

    if n<1:
        raise ValueError('accumulate should be at least 1.')

    stats = pd.read_csv(args.input, sep='\t')
    l=len(stats)

    if l == 0:
        raise ValueError('No alignments stats found.')
        
    # Warn if there are not many chunks to produce meaningful accumulated summary
    if n>1 and l<10*n:
        sys.stdout.write("WARNING! Accumulating "+str(l)+" chunks in groups of "+str(n)+" might not produce meaningful summary stats.\n")

    # Sort by query chromosome and chunk
    if n>1:
        qchr=[re.sub("_chunk.+","",s) for s in stats['name']]
        chunk=[re.sub(".+_chunk","",s) for s in stats['name']]
        # If chunk numbers are present convert to integers for sorting
        if all([r.isdigit() for r in chunk]):
            chunk=[int(r) for r in chunk]
        stats=stats.assign(qchr=qchr,chunk=chunk)
        stats.sort_values(by=["qchr","chunk"],inplace=True)

    def output_summary(stats, prefix=''):
        to_str_opts = {'buf': args.output, 'col_space': 8, 'index':False, 'justify': 'center'}
        errors = pd.DataFrame(summarise_stats(stats, percentiles=args.percentiles, accumulate=n))

        pcnt_formatter = {}
        for col in errors.columns:
            if np.issubdtype(errors[col].dtype, np.number):
                pcnt_formatter[col] = '{:5.3%}'.format
            else:
                pcnt_formatter[col] = '{}'.format
        args.output.write('\n\n# {} Percentage Errors\n'.format(prefix))
        errors.to_string(formatters=pcnt_formatter, **to_str_opts)

        args.output.write('\n\n# {} Q Scores\n'.format(prefix))
        for col in errors.columns:
            if np.issubdtype(errors[col].dtype, np.number):
                errors[col] = qscore(errors[col])
        errors.to_string(float_format='%5.2f', **to_str_opts)

    output_summary(stats)

    args.output.write('\n\n# Ref Coverage\n')
    for ref, grp in stats.groupby('ref'):
        args.output.write('{} {:5.3f}\n'.format(ref, grp['ref_coverage'].sum()))

    if args.per_reference:
        for ref, grp in stats.groupby('ref'):
            output_summary(grp, prefix=ref)

    args.output.write('\n')


if __name__ == '__main__':
    main()
