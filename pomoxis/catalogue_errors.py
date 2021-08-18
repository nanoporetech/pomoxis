import argparse
from collections import defaultdict, Counter, namedtuple
import concurrent.futures
from functools import partial
import itertools
import logging
import shutil
from math import ceil
from operator import attrgetter
import os
import pickle
import re
import unittest
import warnings

import matplotlib; matplotlib.use('Agg', force=True)  # enforce non-interactive backend
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import pysam

from pomoxis.util import get_trimmed_pairs, intervaltrees_from_bed

AlignSeg = namedtuple('AlignSeg', ('rname', 'qname', 'pairs', 'rlen'))
Error = namedtuple('Error', ('rp', 'rname', 'qp', 'qname', 'ref', 'match', 'read', 'counts', 'klass', 'aggr_klass'))
Context = namedtuple('Context', ('p_i', 'qb', 'rb'))

_match_ = ' '
_indel_ = ':'
_sub_ = '|'

_sep_ = '\t'

_indel_sizes_=(5,10,100)
_error_groups_ = [
                  ('HP sub swp', lambda x: 'HP sub swp' in x),
                  ('HP sub ext', lambda x: 'HP sub ext' in x),
                  ('HP sub trc', lambda x: 'HP sub trc' in x),
                  ('sub HP split', lambda x: 'sub HP split' in x),
                  ('ins HP split', lambda x: 'ins HP split' in x),
                  ('sub HP join', lambda x: 'sub HP join' in x),
                  ('del HP join', lambda x: 'del HP join' in x),
                  ('multi ins >', lambda x: bool(re.search('multi .*ins >', x))),
                  ('multi ins <', lambda x: bool(re.search('multi .*ins <', x))),
                  ('multi del >', lambda x: bool(re.search('multi .*del >', x))),
                  ('multi del <', lambda x: bool(re.search('multi .*del <', x))),
                  ('HP del', lambda x: 'HP del' in x),
                  ('HP ins', lambda x: 'HP ins' in x),
                  ('HP sub', lambda x: 'HP sub' in x),
                  ('fwd repeat ins', lambda x: 'fwd repeat ins' in x),
                  ('rev repeat ins', lambda x: 'rev repeat ins' in x),
                  ('fwd repeat del', lambda x: 'fwd repeat del' in x),
                  ('rev repeat del', lambda x: 'rev repeat del' in x),
                  ('ins', lambda x: 'ins' in x),
                  ('del', lambda x: 'del' in x),
                  ('sub', lambda x: 'sub' in x),
]

def get_errors(aln, tree=None):
    """Find positions of errors in an aligment.

    :param aln: iterable of `AlignPos` objects.
    :param bed_file: path to .bed file of regions to include in analysis.
    :param tree: `intervaltree.IntervalTree` object of regions to analyse.
    :returns: ( [(ri, qi, 'error_type', last_ri, last_qi)], aligned_ref_len)
        ri, qi: ref and query positions
        error_type: 'D', 'I' or 'S'
        last_ri, last_qi: ref and query positions of the last match
        aligned_ref_len: total aligned reference length (taking account of masking tree)
    """

    err = []
    last_qi = None
    last_ri = None
    pos = None
    n_masked = 0
    aligned_ref_len = 0
    match = 0
    for (qi, qb, ri, rb) in aln:
        if tree is not None:
            pos = ri if ri is not None else pos
            if not tree.overlaps(pos) or (ri is None and not tree.overlaps(pos + 1)):
                # if ri is None, we are in an insertion, check if pos + 1 overlaps
                # (ref position of ins is arbitrary)
                # print('Skipping ref {}:{}'.format(read.reference_name, pos))
                n_masked += 1
                continue
        if qi is None:  # deletion
            last_ri = ri  # ri will not be None
            err.append((ri, qi, 'D', (last_ri, last_qi)))
            aligned_ref_len += 1
        elif ri is None:
            last_qi = qi  # qi will not be None
            err.append((ri, qi, 'I', (last_ri, last_qi)))
        else:
            last_qi = qi
            last_ri = ri
            if qb != rb:
                err.append((ri, qi, 'S', (last_ri, last_qi)))
            aligned_ref_len += 1
            match += 1

    if match == 0:
        # no matches within bed regions - all bed ref positions were deleted.
        # skip this alignment.
        return None
    
    return err, aligned_ref_len, n_masked


def rle(it):
    """Calculate a run length encoding (rle), of an input vector.

    :param it: iterable.
    :returns: structured array with fields `start`, `length`, and `value`.
    """
    val_dtype = np.array(it[0]).dtype
    dtype = [('length', int), ('start', int), ('value', val_dtype)]

    def _gen():
        start = 0
        for key, group in itertools.groupby(it):
            length = sum(1 for x in group)
            yield length, start, key
            start += length
    return np.fromiter(_gen(), dtype=dtype)


def get_run(i, runs):
    """Find run to which the i'th element belongs.

    :param i: int, element index wihin input to `rle`.
    :returns: int, element index within runs to which i belongs.
    """
    ends = runs['start'] + runs['length']
    start_i = min(np.searchsorted(runs['start'], i), len(runs) - 1)
    end_i = min(np.searchsorted(ends, i), len(runs) - 1)
    start = runs['start'][start_i]
    end = ends[end_i] - 1
    run_i = start_i if np.argmin(np.abs([start - i, end - i])) == 0 else end_i

    return run_i


def _get_context_bounds(p, aln, search_by_q, offset):
    """Find start and en d of context.

    In the simplest case this will be p-offset:p+offset, but we adjust
    for boundaries and to ensure we don't start/end on an error or within a HP.

    :param p: int, position (ref position, or query position)
    :param aln: iterable of `AlignPos` objects.
    :param search_by_q: bool, whether to search by query position (typically done if qi is None).
    :returns: (int start index, int end index, int index of p within aln[start:end])
    """
    if search_by_q:
        pos = [x.qpos for x in aln]
    else:
        pos = [x.rpos for x in aln]


    offset_tmp = offset
    start_p = max(p - offset, pos[0])
    end_p = min(p + offset, pos[-1])
    s = pos.index(start_p)
    e = pos.index(end_p)
    # extend context until it does not start/end with an error or in a HP.
    while True:
        s_is_match = aln[s].qbase == aln[s].rbase
        sq_not_hp = aln[s+1].qbase != aln[s].qbase
        sr_not_hp = aln[s+1].rbase != aln[s].rbase
        if (s_is_match and sq_not_hp and sr_not_hp) or s == 0:
            break
        else:
            s -= 1
    while True:
        e_is_match = aln[e-1].qbase == aln[e-1].rbase
        eq_not_hp = aln[e-1].qbase != aln[e-2].qbase
        er_not_hp = aln[e-1].rbase != aln[e-2].rbase
        if (e_is_match and eq_not_hp and er_not_hp) or e == len(aln):
            break
        else:
            e += 1

    p_i = pos.index(p)

    return s, e, p_i - s


def are_adjacent(inds):
    """"Check if all int indices in the interable are consecutive.

    :param inds: iterable of ints
    :returns: bool
    """
    if len(inds) == 1:
        adjacent = True
    else:
        adjacent = False
        i = inds[0]
        for n in inds[1:]:
            if n - i != 1:
                break
            i = n
        else:
            adjacent = True
    return adjacent


def is_in_hp(seq, p_i):
    n_indels_fwd = len([x for x in itertools.takewhile(lambda x: x=='-', seq[p_i+1:])])
    n_indels_rev = len([x for x in itertools.takewhile(lambda x: x=='-', seq[:p_i][::-1])])
    if seq[p_i] == '-':
        p_is_hp = seq[p_i - n_indels_rev - 1].upper() == seq[p_i + n_indels_fwd + 1].upper()
    else:
        p_is_hp = (seq[p_i].upper() == seq[p_i - 1 - n_indels_rev].upper() or
                    seq[p_i].upper() == seq[p_i + 1 + n_indels_fwd].upper())
    return p_is_hp


def classify_hp_sub(p_i, adjacent, errors, match_line, rb_runs,
                    qb_runs, qp_is_hp, rp_is_hp):
    hp_kls = None
    # sub errors
    # 'swap', e.g. r TTCCC -> q TTTCC or  r CCTTCCC -> q CTTTTCC
    # 'split', e.g. r TTTTTT -> q TTCCTT or r TTTTTT -> q TCTTCT
    # 'join', e.g. r TTCCTT -> TTTTTT
    # 'trunc', e.g. r TTTTTT -> q TTTTTC, r TTTTTT -> q CTTTTC,
    #               r TT -> q TC, r TT -> q GC
    # 'ext', e.g. r TTTTTC -> q TTTTTT
    rp_run_ind = get_run(p_i, rb_runs)
    qp_run_ind = get_run(p_i, qb_runs)
    rp_run = rb_runs[rp_run_ind]
    qp_run = qb_runs[qp_run_ind]
    sub_rb_inds = [get_run(i, rb_runs) for i in errors['sub']]
    sub_qb_inds = [get_run(i, qb_runs) for i in errors['sub']]
    # find query runs ref position run (if if ref is HP, get q over extent)
    q_run_inds_in_r_run = set([get_run(i, qb_runs) for i in range(rp_run['start'], rp_run['start'] + rp_run['length'])])
    q_runs_in_r_run = qb_runs[list(q_run_inds_in_r_run)]
    r_run_inds_in_q_run = set([get_run(i, rb_runs) for i in range(qp_run['start'], qp_run['start'] + qp_run['length'])])
    r_runs_in_q_run = rb_runs[list(r_run_inds_in_q_run)]
    cols = ['value', 'length']
    if qp_run[cols] == rp_run[cols]:
        # this is a false alarm, though some HP in the context has changed,
        # it was not as this position
        return hp_kls
    # if more than 1 query run is in ref run, we have a split / trunc
    # swap is a special case of truncation in which another HP has gained a base
    if len(q_run_inds_in_r_run) > 1:

        # all subs should belong to same RHP, otherwise we have a mess.
        if all([i == rp_run_ind for i in sub_rb_inds]):
            hp_start_ind = rp_run['start']
            hp_end_ind = hp_start_ind + rp_run['length'] - 1
            match_runs = rle(match_line)
            match_sub_inds = [get_run(i, match_runs) for i in errors['sub']]
            # if subs are not flanking we have a split
            if min(errors['sub']) > hp_start_ind and max(errors['sub']) < hp_end_ind:
                hp_kls = 'sub HP split ({}{})'.format(rp_run['value'], rp_run['length'])

            # if (possibly multiple adjacent) subs are all flanking the HP
            elif set(match_sub_inds).issubset([get_run(hp_start_ind, match_runs),
                                                get_run(hp_end_ind, match_runs)]):
                # if query is a HP, this is a swap
                if qp_run['length'] > 1:
                    hp_kls = 'HP sub swp ({}{}->{}{},{})'.format(rp_run['value'], rp_run['length'],
                                                                qp_run['value'], qp_run['length'] - len(errors['sub']), len(errors['sub']))
                else:
                    hp_kls = 'flk HP sub trc ({}{}->{})'.format(rp_run['value'], rp_run['length'],
                                                        rp_run['length'] - len(errors['sub']))
            else:
                hp_kls = 'complex HP sub trc ({}{})'.format(rp_run['value'], rp_run['length'])
        else:
            hp_kls = 'messy HP sub trc'

    # if more than 1 ref run is in query run, we have a join or extension
    # swap is a special case of extension in which another HP has lost a base
    elif len(r_run_inds_in_q_run) > 1:
        # all subs should belong to same QHP, else we have a mess
        if all([i == qp_run_ind for i in sub_qb_inds]):
            hp_start_ind = qp_run['start']
            hp_end_ind = hp_start_ind + qp_run['length'] - 1
            match_runs = rle(match_line)
            match_sub_inds = [get_run(i, match_runs) for i in errors['sub']]
            # if subs are not flanking we have a join
            if min(errors['sub']) > hp_start_ind and max(errors['sub']) < hp_end_ind:
                hp_kls = 'sub HP join ({}{})'.format(qp_run['value'], qp_run['length'])

            # if (possibly multiple adjacent) subs are all flanking the HP
            elif set(match_sub_inds).issubset([get_run(hp_start_ind, match_runs),
                                                get_run(hp_end_ind, match_runs)]):
                # if ref is a HP, this is a swap
                if rp_run['length'] > 1:
                    hp_kls = 'HP sub swp ({}{}->{}{},{})'.format(rp_run['value'], rp_run['length'],
                                                                qp_run['value'], qp_run['length'] - len(errors['sub']), len(errors['sub']))
                else:
                    hp_kls = 'flk HP sub ext ({}{}->{})'.format(qp_run['value'], qp_run['length']- len(errors['sub']), qp_run['length'])
            else:
                hp_kls = 'complex HP sub ext ({}{})'.format(qp_run['value'], qp_run['length'])
        else:
            hp_kls = 'messy HP sub ext'

    # if we just have 1 run of each, we might have a complete sub of the HP
    # e.g. TT->CC
    elif len(r_run_inds_in_q_run) == 1 and len(q_run_inds_in_r_run) == 1:
        if rp_run['value'] != qp_run['value']:
            hp_kls = 'HP sub swp ({}{}->{}{},{})'.format(rp_run['value'], rp_run['length'],
                                                          qp_run['value'], 0, qp_run['length'])

    return hp_kls


def classify_hp_indel(p_i, key, errors, runs1, seq2):
    """Look for a specific kind of HP indel that splits or joins two HPs

    :param p_i: int, index of error
    :param key: key of error type within errors (should be 'ins' or 'del')
    :param errors: dict of error positions
    :runs1: np.ndarray, rle encoding of sequence1 (query if deletions
           join two HPs, or ref if insertions split a HP
    :seq2: iterable of str of sequence2 (ref if deletions join two HPs,
           or query if insertions split a HP.
    :returns: str classification or None
    """
    # look for deletions in query HPs which indicate two ref HPs are joined.
    # or look for insertions in ref HPs which split them
    assert key in ['ins', 'del']
    _type_ = {'ins': 'split', 'del': 'join'}

    hp_kls = None
    p_run_ind = get_run(p_i, runs1)

    # if ref/query is not a HP, return
    if runs1[p_run_ind - 1]['value'] != runs1[p_run_ind + 1]['value']:
        logging.debug('Not a HP {} {}'.format(runs1[p_run_ind - 1]['value'],
                                              runs1[p_run_ind + 1]['value']))
        return hp_kls

    # if we have multiple non-adjacent indels, we could have
    # HP which should be split into >2 runs or HP joined from >2 HPs
    hp_base = runs1[p_run_ind - 1]['value']
    helper = lambda x: runs1[x]['value'] in ['-', hp_base]
    run_inds_in_hp = [i for i in itertools.takewhile(helper, range(p_run_ind - 1, 0, -1))]
    run_inds_in_hp += [i for i in itertools.takewhile(helper, range(p_run_ind, len(runs1)))]
    runs1_in_hp = runs1[run_inds_in_hp]
    runs1_in_hp = runs1_in_hp[np.where(runs1_in_hp['value'] == hp_base)]
    hp_len = np.sum(runs1_in_hp['length'])
    hp_start_ind = runs1_in_hp[0]['start']
    hp_end_ind = runs1_in_hp[-1]['start'] + runs1_in_hp[-1]['length'] - 1

    # check that the HPs in the query and ref actually differ (some insertions
    # might preserve the HP length e.g. CG--GC -> CGGCGC, in which case the HP
    # is not a split/join
    runs2_in_hp = rle(seq2[hp_start_ind: hp_end_ind + 1])
    runs2_in_hp = runs2_in_hp[np.where(runs2_in_hp['length'] > 1)]
    if (len(runs2_in_hp) > 0 and
        np.any(np.logical_and(runs2_in_hp['length'] == hp_len,
                              runs2_in_hp['value'] == hp_base))):

        return hp_kls

    # we have a split/join if any indels are within hp
    if any([i > hp_start_ind and i < hp_end_ind for i in errors[key]]):
        # if any indels are outside HP run, we have a mess
        if min(errors[key]) < hp_start_ind or max(errors[key]) > hp_end_ind:
            subtype = 'messy'
        elif are_adjacent(errors[key]):
            if len(errors[key]) == 1:
                subtype = 'simple'
            else:
                subtype = 'multi'
        else:
            subtype = 'complex'
        hp_kls = '{} {} HP {} ({}{})'.format(subtype, key, _type_[key], hp_base, hp_len)

    return hp_kls


def get_match_line_and_err_index(context):
    match_line = ''
    errors = defaultdict(list)
    k = ''
    for i, (q, r) in enumerate(zip(context.qb, context.rb)):
        if q == r:
            match_line += _match_
        elif q == '-' or r == '-':  # indel
            match_line += _indel_
            k = 'del' if q == '-' else 'ins'
        else:  # sub
            match_line += _sub_
            k = 'sub'
        if match_line[-1] != _match_:
            errors[k].append(i)
        if i == context.p_i:  # this is the central error we are classifying
            p_k = k
    return match_line, errors, p_k


def preprocess_error(p, aln, search_by_q, offset=10):
    """
    :param p: int, position (ref position, or query position)
    :param aln: iterable of `AlignPos` objects.
    :param search_by_q: bool, whether to search by query position (typically done if qi is None).
    :returns: `Context` object
    """
    s, e, p_i = _get_context_bounds(p, aln, search_by_q, offset)
    sl = aln[s:e]
    qi, qb, ri, rb = zip(*sl)

    return Context(p_i, qb, rb)


def simple_klass(adjacent, n, err_type, indel_sizes):
    if not adjacent:
        descr = 'complex'
    else:
        descr = 'simple' if n == 1 else 'multi'
    if 'ins' in err_type or 'del' in err_type:
        size = _get_size(n, indel_sizes)
        klass = "{} {} {}".format(descr, err_type, size)
    else:
        klass = "{} {}".format(descr, err_type)

    return klass


def classify_error(context, indel_sizes=None):
    """Classify error within an alignment.

    :param context: `Context` object
    :indel_sizes: iterable of int, for binning indel sizes.
        indels >= to indel_sizes[0] will not be considered as HP splitting/joining indels
    :returns: (str reference_context,
               str match_line,
               str query_context,
               dict counts of sub/ins/del within context)
    """
    if indel_sizes is None:
        indel_sizes = _indel_sizes_

    p_i, qb, rb = context
    match_line, errors, p_k = get_match_line_and_err_index(context)

    # find homopolymers within context
    qb_runs = rle([b.upper() for b in qb])
    rb_runs = rle([b.upper() for b in rb])
    qb_hps = qb_runs[np.where(qb_runs['length'] > 1)]
    rb_hps = rb_runs[np.where(rb_runs['length'] > 1)]
    cols = ['value', 'length']
    hp_changed = len(qb_hps) != len(rb_hps) or np.any(qb_hps[cols] != rb_hps[cols])
    # hp might also be changed if we have an ins in the middle of a ref hp
    if not hp_changed:
        hp_changed = any([(rb_runs[i]['value'] == '-' and
                           rb_runs[i - 1]['value'] == rb_runs[i + 1]['value'])
                          for i in range(1, len(rb_runs) - 1)])
    # hp might also be changed if we have an del in the middle of a query hp
    if not hp_changed:
        hp_changed = any([(qb_runs[i]['value'] == '-' and
                           qb_runs[i - 1]['value'] == qb_runs[i + 1]['value'])
                          for i in range(1, len(qb_runs) - 1)])

    hp_kls = None

    # Try to class the error, make it up as we go along!
    dels = len(errors['del'])
    ins = len(errors['ins'])
    indels = dels + ins
    subs = len(errors['sub'])
    # is this position part of a hp in query or ref?
    qp_is_hp = is_in_hp(qb, p_i)
    rp_is_hp = is_in_hp(rb, p_i)
    hp_subtypes = {(False, False): '',
                   (True, False): 'RHP',
                   (False, True): 'QHP',
                   (True, True): 'HP',
    }

    hp_subtype = hp_subtypes[(rp_is_hp, qp_is_hp)]
    err_type = '{} {}'.format(hp_subtype, p_k) if hp_subtype != '' else p_k
    kls = 'unknown ({})'.format(err_type)

    if (subs > 0 and indels > 0) or (ins > 0 and dels >0):
        kls = 'complex mess ({})'.format(err_type)

    elif indels == 0:  # we have a sub
        adjacent = are_adjacent(errors['sub'])
        kls = simple_klass(adjacent, subs, err_type, indel_sizes)
        if hp_changed:
            hp_kls = classify_hp_sub(p_i, adjacent, errors, match_line, rb_runs,
                                     qb_runs, qp_is_hp, rp_is_hp)
            kls = hp_kls if hp_kls is not None else kls

    elif dels == 0:
        adjacent = are_adjacent(errors['ins'])
        if hp_changed and ins < indel_sizes[0]:
            hp_kls = classify_hp_indel(p_i, 'ins', errors, rb_runs, qb)
            kls = hp_kls if hp_kls is not None else kls
        if hp_kls is None:
            kls = simple_klass(adjacent, ins, err_type, indel_sizes)
            if adjacent:
                i = errors['ins'][0]
                b = qb[i]
                q_hp_len = len([x for x in itertools.takewhile(lambda x: x==b, qb[i:])])
                r_hp_len = len([x for x in itertools.takewhile(lambda x: x==b, rb[i+ins:])])
                inserted = qb[i:i + ins]
                # get ref_after insertions to check if we have a repeat
                repeat_start = i + ins
                repeat_end = repeat_start + ins
                if repeat_end <= len(rb):
                    ref_after = rb[repeat_start: repeat_end]
                else:  # we don't have enough context after the insertion
                    ref_after = []

                if (q_hp_len > 1 or r_hp_len > 1) and all([qb[j] == b and rb[j] == '-' for j in errors['ins']]):
                    kls = "HP ins ({}{}->{})".format(b, r_hp_len, q_hp_len)
                elif inserted == ref_after:
                    kls = "fwd repeat ins len {}".format(ins)
                elif inserted == ref_after[::-1]:
                    kls = "rev repeat ins len {}".format(ins)

    elif ins == 0:
        adjacent = are_adjacent(errors['del'])
        if hp_changed and dels < indel_sizes[0]:
            hp_kls = classify_hp_indel(p_i, 'del', errors, qb_runs, rb)
            kls = hp_kls if hp_kls is not None else kls
        if hp_kls is None:
            kls = simple_klass(adjacent, dels, err_type, indel_sizes)
            if adjacent:
                i = errors['del'][0]
                b = rb[i]
                # we use aln_rb instead of rb, so we can find the length
                # of HPs which extent beyond the context.
                hp_len = len([x for x in itertools.takewhile(lambda x: x==b, rb[i:])])
                deleted = rb[i:i + dels]
                # get ref_after deletions to check if we have a repeat
                repeat_start = i + dels
                repeat_end = repeat_start + dels
                if repeat_end <= len(rb):
                    ref_after = rb[repeat_start: repeat_end]
                else:  # we don't have enough context after the insertion
                    ref_after = []
                if hp_len > 1:
                    if all([rb[j] == b for j in errors['del']]):
                        kls = "HP del ({}{}->{})".format(b, hp_len, hp_len - dels)
                elif deleted == ref_after:
                    kls = "fwd repeat del len {}".format(dels)
                elif deleted == ref_after[::-1]:
                    kls = "rev repeat del len {}".format(dels)

    return ''.join(rb), match_line, ''.join(qb), errors, kls


def _get_size(n, sizes):
    l_i = np.searchsorted(sizes, n)
    if l_i == len(sizes):
        size = "> {}".format(sizes[-1])
    else:
        size = "<= {}".format(sizes[l_i])
    return size


def _process_read(bam, outdir, read_range, bed_file=None):
    """Load an alignment from bam and return result of `_process_seg`.

    :param bam: str, bam file.
    :param outdir: str, output path
    :param read_range: (int, int, int), process number, range of alignments to process.
    :param bed_file: path to .bed file of regions to include in analysis.
    :returns: result of `_process_seg`.
    """
    trees = None
    if bed_file is not None:
        trees = intervaltrees_from_bed(bed_file)

    total_ref_length = defaultdict(int)
    total_n_ref_sites_masked = defaultdict(int)
    error_count = defaultdict(Counter)

    db_path = os.path.join(outdir, 'error_catalogue_db_{}.txt'.format(read_range[0]))
    txt_path = os.path.join(outdir, 'error_catalogue_{}.txt'.format(read_range[0]))
    db_fh = open(db_path, 'w')
    txt_fh = open(txt_path, 'w')
    headers = get_headers()

    with pysam.AlignmentFile(bam, 'rb') as bam_obj:
        for i, rec in enumerate(bam_obj.fetch(until_eof=True)):
            if i < read_range[1]:
                continue
            elif i == read_range[2]:
                break

            if rec.is_unmapped or rec.is_supplementary or rec.is_secondary:
                return

            if bed_file is not None:
                tree = trees[rec.reference_name]

                if not tree.overlaps(rec.reference_start, rec.reference_end):
                    #sys.stderr.write('read {} does not overlap with any regions in bedfile\n'.format(rec.query_name))
                    continue
            else:
                tree = None

            seg = AlignSeg(rname=rec.reference_name, qname=rec.query_name,
                           pairs=list(get_trimmed_pairs(rec)), rlen=rec.reference_length
                          )
            logging.debug('Loaded query {}'.format(seg.qname))

            seg_result = _process_seg(seg, db_fh, txt_fh, headers, tree)
            if seg_result is None:
                continue

            ref_name, ref_length, counts, n_masked = seg_result
            
            error_count[ref_name].update(counts)
            total_ref_length[ref_name] += ref_length
            total_n_ref_sites_masked[ref_name] += n_masked

    db_fh.close()
    txt_fh.close()

    return error_count, total_ref_length, total_n_ref_sites_masked


def _process_seg(seg, db_fh, txt_fh, headers, tree=None):
    """Classify and count errors within an `AlignSeg` object.

    :param seg: `AlignSeg` object.
    :param db_fh: file object, catalogue db file
    :param txt_fh: file object, catalogue file
    :param headers: list of tuples, error properties and their label ids
    :param tree: `intervaltree.IntervalTree` object of regions to analyse.
    :returns: (seg.rname, aligned_ref_len, error_count, n_masked)
        error_count: `Counter` of error classes
        n_masked: number of reference positions excluded by tree.
    """
    error_count = Counter()
    err_result = get_errors(seg.pairs, tree)
    if err_result is None:  # no matches within bed regions
        logging.debug('Skipping {} since all bed regions were deleted'.format(
                      seg.qname))
        return None

    pos_and_errors, aligned_ref_len, n_masked = err_result
    for ri, qi, error, approx_pos in pos_and_errors:
        ref, match, read, counts, klass = classify_error(preprocess_error(
            ri if ri is not None else qi, seg.pairs, search_by_q=(ri is None)
        ))

        rp = ri
        qp = qi
        if rp is None:
            rp = "~{}".format(approx_pos[0])
        if qp is None:
            qp = "~{}".format(approx_pos[1])

        e = Error(rp=rp, rname=seg.rname, qp=qp, qname=seg.qname, ref=ref,
                  match=match, read=read, counts=counts, klass=klass,
                  aggr_klass=get_aggr_klass(klass))
        error_count[klass] += 1

        db_fh.write(_sep_.join((str(h[1](e)) for h in headers)) + '\n')
        txt_fh.write("Ref Pos: {}, {} Pos {}, {}, {}\n".format(e.rp, e.qname, e.qp, e.klass, e.aggr_klass))
        txt_fh.write(e.ref + "\n")
        txt_fh.write(e.match + "\n")
        txt_fh.write(e.read + "\n")
        txt_fh.write(".\n")

    if tree is None:
        assert seg.rlen == aligned_ref_len

    logging.debug('Done processing {} aligned to {}'.format(seg.qname, seg.rname))
    return seg.rname, aligned_ref_len, error_count, n_masked


def qscore(d):
    """Calculate a qscore"""
    with warnings.catch_warnings():
        # some values might be zero
        warnings.simplefilter("ignore")
        q = -10 * np.log10(d)
    return q


def analyze_counts(counts, total_ref_length):
    df = pd.DataFrame({'klass': list(counts.keys()),
                       'count': list(counts.values())})
    df['err_rate'] = df['count'] / total_ref_length
    df['fraction_of_total'] = df['count'] / df['count'].sum()
    df.sort_values('count', ascending=True, inplace=True)
    df['remaining_count'] = df['count'].cumsum()
    df['remaining_err_rate'] = df['err_rate'].cumsum()
    df['err_rate_q'] = qscore(df['err_rate'])
    df['remaining_err_rate_q'] = qscore(df['remaining_err_rate'])
    df.sort_values('remaining_count', ascending=False, inplace=True)
    return df


def plot_summary(df, outdir, prefix, ref_len):
    """Create a plot showing Q-scores as largest remaining
    error klass is removed"""

    fig, ax = plt.subplots()
    fig.subplots_adjust(left=0.3)
    y_pos = np.arange(len(df) + 1)
    no_error_score = -10 * np.log10(1/ref_len)
    ax.barh(y_pos, df['remaining_err_rate_q'].append(pd.Series(no_error_score)), align='center', color='green', ecolor='black')
    ax.set_xlabel('Q(Accuracy)')
    ax.set_ylabel('Error Class')
    ax.set_ylim((y_pos[0]-0.5, y_pos[-1]+0.5))
    ax.set_yticks(y_pos)
    ax.set_yticklabels(['total error'] + list(df['klass']))
    ax.invert_yaxis()  # labels read top-to-bottom
    xstart, xend = ax.get_xlim()
    ystart, yend = ax.get_ylim()
    ax.text(xend - 2.25, ystart - 0.25, '+')
    ax.set_title('Q-score after removing error class')
    fp = os.path.join(outdir, '{}_remaining_errors.png'.format(prefix))
    fig.savefig(fp)
    plt.close()


def get_aggr_counts(total_counts):
    # get errors per type
    aggregate_counts = Counter()
    max_indel = max(_indel_sizes_)
    for key, val in total_counts.items():
        aggregate_counts[get_aggr_klass(key)] += val
    return aggregate_counts


def get_aggr_klass(klass):
    max_indel = max(_indel_sizes_)
    for error_type, is_type in _error_groups_:
        if is_type(klass):
            if '>' in error_type or '<' in error_type:
                aggr_klass = '{} {}'.format(error_type, max_indel)
            else:
                aggr_klass = error_type
            break
    return aggr_klass


def analyse_errors(args):
    # load and merge existing counts
    error_count = defaultdict(Counter)
    total_ref_length = Counter()
    total_n_ref_sites_masked = Counter()
    for pkl in args.pkl:
        with open(pkl, 'rb') as fh:
            counts = pickle.load(fh)
        
        # pkl format:
        # {'ref_lengths': total_ref_length,
        #  'n_ref_sites_masked': total_n_ref_sites_masked,
        #  'counts': {'by_ref': error_count,
        #             'by_ref_aggr': aggr_by_ref,
        #             'total': total_counts,
        #             'total_aggr': aggregate_counts,
        #            }
        # }
        for ref, d_counts in counts['counts']['by_ref'].items():
            error_count[ref].update(d_counts)
        total_ref_length.update(counts['ref_lengths'])
        total_n_ref_sites_masked.update(counts['n_ref_sites_masked'])

    aggr_and_output(args, error_count, total_ref_length, total_n_ref_sites_masked)


def get_headers():
    # make an approximate position into int
    helper = lambda x: int(x.replace('~','')) if isinstance(x, str) else x

    return [('ref_name', attrgetter('rname')),
            ('ref_pos', helper(attrgetter('rp'))),
            ('ref_context', attrgetter('ref')),
            ('query_name', attrgetter('qname')),
            ('query_pos', helper(attrgetter('qp'))),
            ('query_context', attrgetter('read')),
            ('class', attrgetter('klass')),
            ('aggr_class', attrgetter('aggr_klass')),
            ('n_ins', lambda e: len(e.counts['ins'])),
            ('n_del', lambda e: len(e.counts['del'])),
            ('n_sub', lambda e: len(e.counts['sub'])),
            ('context_len', lambda e: len(e.ref)),
    ]


def count_errors(args):
    # create a slice of reads to process in each thread to avoid looping through
    # bam n read times and reduce mp overhead
    ranges = [(0, 0, float('inf'))]
    if args.threads > 1:
        with pysam.AlignmentFile(args.bam) as bam:
            n_reads = bam.count(until_eof=True)
        n_reads_per_proc = ceil(n_reads / args.threads)
        ranges = [(thread, thread * n_reads_per_proc, min((thread + 1) * n_reads_per_proc, n_reads))
                   for thread in range(args.threads)]

    total_ref_length = defaultdict(int)
    total_n_ref_sites_masked = defaultdict(int)
    error_count = defaultdict(Counter)

    # record draft start position of each long multi indel
    multi_errs = {}

    f = partial(_process_read, args.bam, args.outdir, bed_file=args.bed)

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as ex:
        for returned in ex.map(f, ranges):
            if returned is None:
                continue
            else:
                counts, ref_length, n_masked = returned

            for ref_name in counts:
                error_count[ref_name].update(counts[ref_name])
                total_ref_length[ref_name] += ref_length[ref_name]
                total_n_ref_sites_masked[ref_name] += n_masked[ref_name]

    db_fh = open(os.path.join(args.outdir, 'error_catalogue_db.txt'), 'w')
    db_fh.write(_sep_.join((h[0] for h in get_headers())) + '\n')
    merge_catalogues(db_fh, [os.path.join(args.outdir, 'error_catalogue_db_{}.txt'.format(n)) for n in range(args.threads)])
    db_fh.close()

    txt_fh = open(os.path.join(args.outdir, 'error_catalogue.txt'), 'w')
    merge_catalogues(txt_fh, [os.path.join(args.outdir, 'error_catalogue_{}.txt'.format(n)) for n in range(args.threads)])
    txt_fh.close()

    aggr_and_output(args, error_count, total_ref_length, total_n_ref_sites_masked)


def merge_catalogues(outfile, filenames):
    for f in filenames:
        try:
            with open(f, 'r') as infile:
                shutil.copyfileobj(infile, outfile)
        except:
            logging.info("Error merging file {}".format(f))
        else:
            logging.info("Merged file {}, deleting".format(f))
            os.remove(f)


def aggr_and_output(args, error_count, total_ref_length, total_n_ref_sites_masked):
    total_counts = Counter()
    aggr_by_ref = {}
    for ref_name, counts in error_count.items():
        df = analyze_counts(counts, total_ref_length[ref_name])
        fp = os.path.join(args.outdir, '{}_error_summary.txt'.format(ref_name))
        df.to_csv(fp, sep=_sep_, index=False)

        aggr_by_ref[ref_name] = get_aggr_counts(counts)
        df = analyze_counts(aggr_by_ref[ref_name], total_ref_length[ref_name])
        fp = os.path.join(args.outdir, '{}_aggr_error_summary.txt'.format(ref_name))
        df.to_csv(fp, sep=_sep_, index=False)
        plot_summary(df, args.outdir, '{}_aggr'.format(ref_name),
                     ref_len=total_ref_length[ref_name])

        total_counts.update(counts)

    df = analyze_counts(total_counts, sum(total_ref_length.values()))
    fp = os.path.join(args.outdir, '{}_error_summary.txt'.format('total'))
    df.to_csv(fp, sep=_sep_, index=False)

    aggregate_counts = get_aggr_counts(total_counts)
    df = analyze_counts(aggregate_counts, sum(total_ref_length.values()))
    fp = os.path.join(args.outdir, '{}_aggr_error_summary.txt'.format('total'))
    df.to_csv(fp, sep=_sep_, index=False)
    plot_summary(df, args.outdir, '{}_aggr'.format('total'),
                 ref_len=sum(total_ref_length.values()))

    # save counts to yaml for any further analysis
    to_save = {'ref_lengths': total_ref_length,
               'n_ref_sites_masked': total_n_ref_sites_masked,
               'counts': {'by_ref': error_count,
                          'by_ref_aggr': aggr_by_ref,
                          'total': total_counts,
                          'total_aggr': aggregate_counts,
                         }
               }
    with open(os.path.join(args.outdir, 'counts.pkl'), 'wb') as fh:
         pickle.dump(to_save, fh)


def main():
    logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    parser = argparse.ArgumentParser(
        prog='catalogue_errors',
        description='Create a catalogue of all query errors in a bam.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(
        title='subcommands', description='valid commands',
        help='additional help', dest='command')
    subparsers.required = True
    
    cparser = subparsers.add_parser('count',
        help='Count query errors in a bam. ',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    cparser.set_defaults(func=count_errors)
    cparser.add_argument('bam', help='Input alignments (aligned to ref).')
    cparser.add_argument('--bed', default=None, help='.bed file of reference regions to include.')
    cparser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads for parallel execution.')
    cparser.add_argument('-o', '--outdir', default='error_catalogue', help='Output directory.')

    aparser = subparsers.add_parser('analyse',
        help='Analyse existing counts, optionally merging multiple counters.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    aparser.set_defaults(func=analyse_errors)
    aparser.add_argument('pkl', nargs='+', help='Input .pkl file(s).')
    aparser.add_argument('-o', '--outdir', default='error_catalogue',
        help="Output directory (will be created).")

    args = parser.parse_args()
    os.mkdir(args.outdir)
    args.func(args)
    logging.info('All done, check {} for output.'.format(args.outdir))


if __name__ == '__main__':
    main()


class ClassifyErrorTest(unittest.TestCase):

    def setUp(self):
        pass


    def test_hp_del_6_5(self):
        rb = 'ACAACAGCAGAAAAAACAGGA'
        qb = 'ACAACAGCAG-AAAAACAGGA'
        p_i = 10
        expected = 'HP del (A6->5)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_hp_del_2_0(self):
        rb = 'CACTTTCGGCTTGAGGATCA'
        qb = 'CACTTTCGGC--GAGGATCA'
        p_i = 10
        expected = 'HP del (T2->0)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_multi_ins(self):
        rb = 'ATGTAATGCC---AAGCTTA'
        qb = 'ATGTAATGCCAGAAAGCTT'
        p_i = 10
        expected = 'multi ins <= 5'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_multi_del(self):
        rb = 'ATGTAATGCCAGAAAGCTT'
        qb = 'ATGTAATGCC---AAGCTTA'
        p_i = 10
        expected = 'multi del <= 5'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_hp_ins(self):
        rb = 'CACCTGGTGC-AAAAGAGAG'
        qb = 'CACCTGGTGCAAAAAGAGAG'
        p_i = 10
        expected = 'HP ins (A4->5)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_HP_sub_split(self):
        rb = 'TAATCTGGCCcCTGCAATGC'
        qb = 'TAATCTGGCCTCTGCAATGC'
        p_i = 10
        expected = 'sub HP split (C4)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_HP_sub_swap_trc(self):
        rb = 'ACTGCGTACCtTTTGTATAAT'
        qb = 'ACTGCGTACCCTTTGTATAAT'
        p_i = 10
        expected = 'HP sub swp (T4->C2,1)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_HP_sub_swap_trc2(self):
        rb = 'ACTGCGTACCttTTGTATAAT'
        qb = 'ACTGCGTACCCCTTGTATAAT'
        p_i = 10
        expected = 'HP sub swp (T4->C2,2)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_HP_sub_swap_ext(self):
        rb = 'ACTGCGTACcTTTTGTATAAT'
        qb = 'ACTGCGTACTTTTTGTATAAT'
        p_i = 9
        expected = 'HP sub swp (T4->C2,1)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_HP_sub_swap_ext(self):
        rb = 'ACTGCGTAccTTTTGTATAAT'
        qb = 'ACTGCGTATTTTTTGTATAAT'
        p_i = 9
        expected = 'HP sub swp (C2->T4,2)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_HP_swap_sub_trc(self):
        rb = 'CGGGCCTTCCCCtTGCCATTCA'
        qb = 'CGGGCCTTCCCCCTGCCATTCA'
        p_i = 12
        expected = 'HP sub swp (T2->C4,1)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_HP_flk_sub_trc(self):
        rb = 'CGAGAAAATCgGGATCGTTG'
        qb = 'CGAGAAAATCAGGATCGTTG'
        p_i = 10
        expected = 'flk HP sub trc (G3->2)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_HP_flk_sub_ext(self):
        rb = 'ATGCAACAAGcTTACGCTGC'
        qb = 'ATGCAACAAGTTTACGCTGC'
        p_i = 10
        expected = 'flk HP sub ext (T2->3)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_HP_complex_sub(self):
        rb = 'CGAGAAAAAAAGGATCGTTG'
        qb = 'CGAGTATAAAAGGATCGTTG'
        p_i = 4
        expected = 'complex HP sub trc (A7)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_HP_messy_sub(self):
        rb = 'CTAaACtGCcGTG'
        qb = 'CTACACCGCTGTG'
        p_i = 3
        expected = 'messy HP sub trc'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_HP_sub_join(self):
        rb = 'AGGAACGAATcTCTGAAGCG'
        qb = 'AGGAACGAATTTCTGAAGCG'
        p_i = 10
        expected = 'sub HP join (T3)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_HP_not_complete_sub(self):
        rb = 'AGGAACGAATTTCTGAAGCG'
        qb = 'AGGAACGAACCCCTGAAGCG'
        p_i = 10
        expected = 'HP sub swp (T3->C1,3)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_HP_complete_sub(self):
        rb = 'AGGAACGAATTTGTGAAGCG'
        qb = 'AGGAACGAACCCGTGAAGCG'
        p_i = 10
        expected = 'HP sub swp (T3->C0,3)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_HP_messy_sub_ext(self):
        rb = 'CGGGTCTTTTctTTTTCATC'
        qb = 'CGGGTCTTTTTCTTTTCATC'
        p_i = 10
        expected = 'messy HP sub ext'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_HP_del_join(self):
        rb = 'CGGGTCTTTTCTTTTCATC'
        qb = 'CGGGTCTTTT-TTTTCATC'
        p_i = 10
        expected = 'simple del HP join (T8)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_HP_del_join2(self):
        rb = 'CGGGTCTTTTCCTTTTCATC'
        qb = 'CGGGTCTTTT--TTTTCATC'
        p_i = 10
        expected = 'multi del HP join (T8)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_HP_del_join3(self):
        rb = 'CGGGTCTTTTCATTTCATC'
        qb = 'CGGGTCTTTT--TTTCATC'
        p_i = 10
        expected = 'multi del HP join (T7)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_HP_del_join4(self):
        rb = 'CGGGTCTTTTCATTTCATC'
        qb = 'CGG-TCTTTT--TTTCATC'
        p_i = 10
        expected = 'messy del HP join (T7)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_HP_ins_split(self):
        rb = 'CGGGTCTT-TTTCATC'
        qb = 'CGGGTCTTGTTTCATC'
        p_i = 8
        expected = 'simple ins HP split (T5)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_HP_ins_split2(self):
        rb = 'CGGGTCTT--TTTCATC'
        qb = 'CGGGTCTTGGTTTCATC'
        p_i = 8
        expected = 'multi ins HP split (T5)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_HP_ins_split3(self):
        rb = 'CGGGTCTT--TTTCATC'
        qb = 'CGGGTCTTGATTTCATC'
        p_i = 8
        expected = 'multi ins HP split (T5)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_HP_ins_split4(self):
        rb = 'CGGGTCTT--TT-TCATC'
        qb = 'CGGGTCTTGATTGTCATC'
        p_i = 8
        expected = 'complex ins HP split (T5)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_HP_ins_split5(self):
        rb = 'CGGG-TCTT--TT-TCATC'
        qb = 'CGGGGTCTTGATTGTCATC'
        p_i = 9
        expected = 'messy ins HP split (T5)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_fwd_repeat_ins(self):
        rb = 'ACCTATAACG--GCGCGCTG'
        qb = 'ACCTATAACGGCGCGCGCTG'
        p_i = 10
        expected = 'fwd repeat ins len 2'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_fwd_repeat_ins(self):
        rb = 'ACCTATAACG--GCGCGCTG'
        qb = 'ACCTATAACGCGGCGCGCTG'
        p_i = 10
        expected = 'rev repeat ins len 2'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_HP_ins_split6(self):
        rb = 'GCCGATTTTT-TCTCCCGTA'
        qb = 'GCCGATTTTTCTCTCCCGTA'
        p_i = 10
        expected = 'simple ins HP split (T6)'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_complex_HP_del(self):
        rb = 'AGGGGGGGGACTTGAACCCCCACGTC'
        qb = 'A-GGGGGGGACTTGAA-CCCCACGTC'
        p_i = 16
        expected = 'complex RHP del <= 5'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_long_multi_del(self):
        rb = 'ACCCACACACCACACCCACACACCACACCCACACCACACCCACACCACACCCACACACCACACCCACACCACACCCACACACCACACCCACACACCACACCCACACCACACCCACACCACACCCACACACCACACCACACCCACACACCCACACACCACACACTCTCTTACATCTACCTCTACTCTCGCTGTCACACCTTACCCGGCTTTCTGACCGAAATTAAAAAAAATGAAAATGAAATCCTCTTCTTTAGCCCTACAACACTTTTACATAGCCCTAAATAGCCCTAAATAGCCCTCATGTACGTCTCCTCCAAGCCCTATTGACTCTTACCCGGAGTTTCAGCTAAAGGCTATACTTACT'
        qb = 'ACCCACACACCA---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T'
        p_i = 12
        expected = 'multi del > 100'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)

    def test_long_multi_ins(self):
        rb = 'ACCCACACACCA---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T'
        qb = 'ACCCACACACCACACCCACACACCACACCCACACCACACCCACACCACACCCACACACCACACCCACACCACACCCACACACCACACCCACACACCACACCCACACCACACCCACACCACACCCACACACCACACCACACCCACACACCCACACACCACACACTCTCTTACATCTACCTCTACTCTCGCTGTCACACCTTACCCGGCTTTCTGACCGAAATTAAAAAAAATGAAAATGAAATCCTCTTCTTTAGCCCTACAACACTTTTACATAGCCCTAAATAGCCCTAAATAGCCCTCATGTACGTCTCCTCCAAGCCCTATTGACTCTTACCCGGAGTTTCAGCTAAAGGCTATACTTACT'
        p_i = 12
        expected = 'multi ins > 100'
        found = classify_error(Context(p_i=p_i, qb=list(qb), rb=list(rb)))[-1]
        self.assertEqual(found, expected)
