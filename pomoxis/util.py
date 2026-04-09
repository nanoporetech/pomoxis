import argparse
from collections import namedtuple, defaultdict
import itertools
import logging
import os

from intervaltree import IntervalTree
import numpy as np
import pandas as pd
import pysam

Region = namedtuple('Region', 'ref_name start end')
AlignPos = namedtuple('AlignPos', ('qpos', 'qbase', 'rpos', 'rbase'))


def parse_regions(regions, ref_lengths=None):
    """Parse region strings into `Region` objects.

    :param regions: iterable of str
    :param ref_lengths: {str ref_names: int ref_lengths}, if provided Region.end
        will default to the reference length instead of None.

    >>> parse_regions(['Ecoli'])[0]
    Region(ref_name='Ecoli', start=0, end=None)
    >>> parse_regions(['Ecoli:1000-2000'])[0]
    Region(ref_name='Ecoli', start=1000, end=2000)
    >>> parse_regions(['Ecoli:-1000'])[0]
    Region(ref_name='Ecoli', start=0, end=1000)
    >>> parse_regions(['Ecoli:500-'])[0]
    Region(ref_name='Ecoli', start=500, end=None)
    >>> parse_regions(['Ecoli'], ref_lengths={'Ecoli':4800000})[0]
    Region(ref_name='Ecoli', start=0, end=4800000)
    >>> parse_regions(['NC_000921.1:10000-20000'])[0]
    Region(ref_name='NC_000921.1', start=10000, end=20000)
    """
    def _get_end(end, ref_name):
        return end if end is not None else (
            ref_lengths[ref_name] if ref_lengths is not None else None
        )

    decoded = []
    for region in regions:
        if ':' not in region:
            ref_name, start, end = region, 0, None
        else:
            start, end = 0, None
            ref_name, bounds = region.split(':')
            if bounds[0] == '-':
                start = 0
                end = int(bounds[1:])
            elif bounds[-1] == '-':
                start = int(bounds[:-1])
                end = None
            else:
                start, end = [int(b) for b in bounds.split('-')]
        decoded.append(Region(ref_name, start, _get_end(end, ref_name)))
    return tuple(decoded)

def get_pairs(aln):
    """Return generator of pairs.

    :param aln: `pysam.AlignedSegment` object.
    :returns: generator of `AlignPos` objects.
    """
    seq = aln.query_sequence
    pairs = (AlignPos(qpos=qp,
                      qbase=seq[qp] if qp is not None else '-',
                      rpos=rp,
                      rbase=rb if rp is not None else '-'
                      )
             for qp, rp, rb in aln.get_aligned_pairs(with_seq=True)
             )
    return pairs


def get_trimmed_pairs(aln):
    """Trim aligned pairs to the alignment.

    :param aln: `pysam.AlignedSegment` object
    :yields pairs:
    """

    pairs = get_pairs(aln)
    for pair in itertools.dropwhile(lambda x: x.rpos is None or x.qpos is None, pairs):
        if (pair.rpos == aln.reference_end or pair.qpos == aln.query_alignment_end):
            break
        yield pair


def yield_from_bed(bedfile):
    with open(bedfile) as fh:
        for line in fh:
            if line.lstrip()[0] == '#':
                continue
            split_line = line.split()
            if split_line[0] in {'browser', 'track'} or len(split_line) < 3:
                continue
            chrom = split_line[0]
            start = int(split_line[1])
            stop = int(split_line[2])
            yield chrom, start, stop


def intervaltrees_from_bed(path_to_bed):
    """Created dict of intervaltrees from a .bed file, indexed by chrom.

    :param path_to_bed: str, path to .bed file.
    :returns: { str chrom: IntervalTreeWrapper obj }.
    """
    trees = defaultdict(IntervalTree)
    for chrom, start, stop in yield_from_bed(path_to_bed):
        trees[chrom].addi(start, stop)
    return dict(trees)


def tag_bam():
    """Command line tool to add tags to a bam."""
    parser = argparse.ArgumentParser(
        prog='tag_bam',
        description='Add a tag to all alignments in a bam.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input', help='Input bam file.')
    parser.add_argument('output', help='Output output file.')
    parser.add_argument('tag_name', help='Tag name.')
    parser.add_argument('tag_value', type=int, help='Tag value.')
    args = parser.parse_args()
    with pysam.AlignmentFile(args.input) as bam_in:
        with pysam.AlignmentFile(args.output, 'wb', header=bam_in.header) as bam_out:
            for r in bam_in:
                r.set_tag(args.tag_name, args.tag_value)
                bam_out.write(r)


def reverse_bed():
    """Convert bed-file coordinates to coordinates on the reverse strand."""
    parser = argparse.ArgumentParser(
        prog='reverse_bed',
        description='Convert bed-file coordinates to coordinates on the reverse strand.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('bed_in', help='Input bed file.')
    parser.add_argument('ref_fasta', help='Input reference fasta file.')
    parser.add_argument('bed_out', help='Output bed file.')
    args = parser.parse_args()

    fasta = pysam.FastaFile(args.ref_fasta)
    lengths = dict(zip(fasta.references, fasta.lengths))
    d = pd.read_csv(args.bed_in, sep='\t', names=['chrom', 'start', 'stop'])

    d['chrom_length'] = d['chrom'].map(lambda x: lengths[x])
    d['rc_stop'] = d['chrom_length'] - d['start']
    d['rc_start'] = d['chrom_length'] - d['stop']
    d['chrom_rc'] = d['chrom'] + '_rc'
    d[['chrom_rc', 'rc_start', 'rc_stop']].to_csv(args.bed_out, index=False, header=False, sep='\t')


class PrimSupAction(argparse.Action):
    """Parse primary / supplementary option."""

    def __call__(self, parser, namespace, values, option_string=None):
        if self.dest == 'primary_only':
            setattr(namespace, self.dest, True)
            setattr(namespace, 'keep_supplementary', False)
        elif self.dest == 'keep_supplementary':
            setattr(namespace, self.dest, True)
            setattr(namespace, 'primary_only', False)


def stats_from_aligned_read(read):
    """Create summary information for an aligned read

    :param read: :class:`pysam.AlignedSegment` object
    """
    try:
        NM = read.get_tag('NM')
    except:
        raise IOError("Read is missing required 'NM' tag. Try running 'samtools fillmd -S - ref.fa'.")

    name = read.query_name
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

    if read.query_qualities is None:
        mean_quality = None
    else:
        mean_quality = round(mean_q_from_qualities(read.query_qualities), ndigits=2)

    results = {
        "name": name,
        "coverage": coverage,
        "direction": direction,
        "length": length,
        "read_length": read_length,
        "mean_quality": mean_quality,
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
        "ref": read.reference_name,
        "aligned_ref_len": read.reference_length,
        "ref_coverage": 100*float(read.reference_length) / read.header.lengths[read.reference_id],
        "mapq": read.mapping_quality,
        "flag": read.flag,
    }

    return results, lra_flag


def masked_stats_from_aligned_read(read, tree):
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
        if not tree.has_overlap(pos, pos + 1) or (rp is None and not tree.has_overlap(pos + 1, pos + 2)):
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

    name = read.query_name
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

    if read.query_qualities is None:
        mean_quality = None
    else:
        mean_quality = round(mean_q_from_qualities(read.query_qualities), ndigits=2)


    results = {
        "name": name,
        "coverage": coverage,
        "direction": direction,
        "length": length,
        "read_length": read_length,
        "mean_quality": mean_quality,
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
        "ref": read.reference_name,
        "aligned_ref_len": aligned_ref_len,
        "ref_coverage": 100 * float(aligned_ref_len) / read.header.lengths[read.reference_id],
        "mapq": read.mapping_quality,
        "flag": read.flag,
        "masked": masked,
    }
    return results


def filter_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)
    group = parser.add_argument_group('Read filtering options')
    group.add_argument('-O', '--orientation', choices=['fwd', 'rev'],
        help='Sample only forward or reverse reads.')
    group.add_argument('-q', '--quality', type=float,
        help='Filter reads by mean qscore.')
    group.add_argument('-a', '--accuracy', type=float,
        help='Filter reads by accuracy.')
    group.add_argument('-c', '--coverage', type=float,
        help='Filter reads by coverage (what fraction of the read aligns).')
    group.add_argument('-l', '--length', type=int, default=None,
        help='Filter reads by read length.')
    group.add_argument('--keep_unmapped', action='store_true',
        help='Include unmapped reads.')

    pgroup = group.add_mutually_exclusive_group()
    pgroup.add_argument('--primary_only', action=PrimSupAction, default=True, nargs=0,
                        help='Use only primary reads.')
    pgroup.add_argument('--keep_supplementary', action=PrimSupAction, default=False, nargs=0,
        help='Include supplementary alignments.')
    return parser


QSCORES_TO_PROBS = 10 ** (-0.1 * np.array(np.arange(100)))

def mean_q_from_qualities(qualities):
    return -10 * np.log10(np.mean(QSCORES_TO_PROBS[qualities]))


def filter_read(r, args, logger=None):
    """Decide whether a read should be filtered out, returning a bool"""

    if logger is None:
        logger = logging.getLogger('Filter')
    # filter secondary and unmapped reads
    if r.is_secondary:
        return True
    if r.is_unmapped and not args.keep_unmapped:
        return True
    if r.is_supplementary and not args.keep_supplementary:
        return True

    # Filters that apply on unaligned reads
    # filter quality
    if args.quality is not None:
        mean_q = mean_q_from_qualities(r.query_qualities)
        if mean_q < args.quality:
            logger.debug(f"Filtering {r.query_name} with quality {mean_q:.2f}")
            return True

    # filter length
    if args.length is not None:
        read_length = r.query_length
        if read_length < args.length:
            logger.info("Filtering {} by length ({:.2f}).".format(r.query_name, read_length))
            return True

    # Filters that only apply on aligned reads
    if not r.is_unmapped:
        # filter orientation
        if (r.is_reverse and args.orientation == 'fwd') or \
                (not r.is_reverse and args.orientation == 'rev'):
            return True

        if args.accuracy is not None or args.coverage is not None:
            stats, _ = stats_from_aligned_read(r)
            if args.accuracy is not None and stats['acc'] < args.accuracy:
                logger.info("Filtering {} by accuracy ({:.2f}).".format(r.query_name, stats['acc']))
                return True
            if args.coverage is not None and stats['coverage'] < args.coverage:
                logger.info("Filtering {} by coverage ({:.2f}).".format(r.query_name, stats['coverage']))
                return True
    # don't filter
    return False


def write_bam(bam, prefix, region, sequences, keep_supplementary=False):
    """Write subset bam.

    :param bam: str, path to input bam.
    :param prefix: str, prefix for output bam.
    :param bam: `Region` obj of region to subset to.
    :param sequences: set of query names to retain.
    :param keep_supplementary: bool, whether to retain supplementary alignments.
    """
    # filtered bam
    sequences = set(sequences)
    taken = set()
    output = '{}_{}.{}'.format(prefix, region.ref_name, os.path.basename(bam))
    src_bam = pysam.AlignmentFile(bam, "rb")
    out_bam = pysam.AlignmentFile(output, "wb", template=src_bam)
    for read in src_bam.fetch(region.ref_name, region.start, region.end):
        if read.is_secondary:
            continue
        if read.is_supplementary and not keep_supplementary:
            continue
        if read.query_name in sequences and (read.query_name, read.is_supplementary) not in taken:
            out_bam.write(read)
            taken.add((read.query_name, read.is_supplementary))
    src_bam.close()
    out_bam.close()
    pysam.index(output)
