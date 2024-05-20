import argparse
from collections import namedtuple, defaultdict
import itertools
import logging
import os
import shutil
import sys

from ncls import NCLS
import numpy as np
import pandas as pd
import pysam

Region = namedtuple('Region', 'ref_name start end')
AlignPos = namedtuple('AlignPos', ('qpos', 'qbase', 'rpos', 'rbase'))


class FastxWrite:
    def __init__(self, fname, mode='w', width=80, force_q=False, mock_q=10):
        self.fname = fname
        self.mode = mode
        self.width = width
        self.force_q = force_q
        self.mock_q = chr(33 + mock_q)
        self.aq = None

    def __enter__(self):
        if self.fname == '-':
            self.fh = sys.stdout
        else:
            self.fh = open(self.fname, self.mode)
        return self

    def __exit__(self, type, value, traceback):
        if self.fname != '-':
            self.fh.close()

    def write(self, name, seq, qual=None, comment=''):
        if comment != '':
            comment = ' {}'.format(comment)
        if (self.aq == 'a' and qual is not None) or \
                (self.aq == 'q' and qual is None):
            raise ValueError("Type of records changed whilst writing.")
        if qual or self.force_q:
            if qual is None:
                qual = self.mock_q * len(seq)
            self.aq = 'q'
            self.fh.write("@{}{}\n{}\n+\n{}\n".format(name, comment, seq, qual))
        else:
            self.aq = 'a'
            self.fh.write(">{}{}\n".format(name, comment))
            for _, chunk in enumerate(chunks(seq, self.width)):
                self.fh.write('{}\n'.format(''.join(chunk)))


def chunks(iterable, n):
    """Generate fixed length chunks of an interable.

    :param iterable: input sequence.
    :param n: chunk size.
    """
    it = iter(iterable)
    while True:
        chunk_it = itertools.islice(it, n)
        try:
            first_el = next(chunk_it)
        except StopIteration:
            return
        yield itertools.chain((first_el,), chunk_it)


def cat(files, output, chunks=1024*1024*10):
    """Concatenate a set of files.

    :param files: input filenames.
    :param output: output filenames.
    :param chunks: buffersize for filecopy.
    """
    with open(output, 'wb') as wfd:
        for f in files:
            with open(f, 'rb') as fd:
                shutil.copyfileobj(fd, wfd, chunks)


def split_fastx(fname, output, chunksize=10000):
    """Split records in a fasta/q into fixed lengths.

    :param fname: input filename.
    :param output: output filename.
    :param chunksize: (maximum) length of output records.
    """
    with FastxWrite(output, 'w') as fout:
        with pysam.FastxFile(fname, persist=False) as fin:
            for rec in fin:
                name = rec.name
                seq = rec.sequence
                qual = rec.quality
                if rec.comment is None:
                    comment = 'chunk_length={}'.format(chunksize)
                else:
                    comment = '{} chunk_length={}'.format(rec.comment, chunksize)
                if qual is None:
                    for i, s in enumerate(chunks(seq, chunksize)):
                        chunk_name = '{}_chunk{}'.format(name, i)
                        fout.write(chunk_name, s, comment=comment)
                else:
                    for i, (s, q) in enumerate(zip(chunks(seq, chunksize), chunks(qual, chunksize))):
                        chunk_name = '{}_chunk{}'.format(name, i)
                        fout.write(chunk_name, ''.join(s), ''.join(q), comment)


def split_fastx_cmdline():
    """Split records in a fasta/q file into chunks of a maximum size."""
    parser = argparse.ArgumentParser(
        prog='split_fastx',
        description='Split records in a fasta/q file into chunks of a maximum size.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input', help='Input fastax/q file.')
    parser.add_argument('output', help='Output fastax/q file.')
    parser.add_argument('chunksize', type=int, help='Maximum size of output sequences.')
    args = parser.parse_args()
    split_fastx(args.input, args.output, args.chunksize)


def fast_convert():
    """Convert between fasta<->fastq."""
    parser = argparse.ArgumentParser(
        prog='fast_convert',
        description='Convert between fasta<->fastq.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('convert', choices=['qq', 'aa', 'aq', 'qa'],
        help='Conversion code: from->to.')
    parser.add_argument('--discard_q', action='store_true',
        help='Discard quality information from fastq, use with --mock_q.')
    parser.add_argument('--mock_q', default=10, type=int,
        help='Mock quality value, valid for convert=aq|qq.')
    args = parser.parse_args()

    in_fmt = 'fastq'
    out_fmt= 'fasta'
    qflag = False # controls quality manipulation
    if args.convert == 'qq':
        out_fmt = 'fastq'
        if args.discard_q is not None:
            qflag = True
    elif args.convert == 'aa':
        in_fmt = 'fasta'
    elif args.convert == 'aq':
        in_fmt = 'fasta'
        out_fmt = 'fastq'
        qflag = True
    elif args.convert == 'qa':
        pass # default
    else:
        raise ValueError("convert must be 'qq', 'aq', 'qa', or 'aa'\n")

    with FastxWrite('-', force_q=out_fmt=='fastq', mock_q=args.mock_q) as fh_out:
        with pysam.FastxFile('-') as fh_in:
            for rec in fh_in:
                if out_fmt == 'fasta':
                    qual = None
                elif qflag:
                    qual = chr(33 + args.mock_q) * len(rec.sequence)
                else:
                    qual = rec.quality

                fh_out.write(rec.name, rec.sequence, qual, rec.comment)


def extract_long_reads():
    """Filter fastq to longest reads."""

    parser = argparse.ArgumentParser(description='Extract longest reads from a fastq.')
    parser.add_argument('input',
        help='Input .fasta/q file.')
    parser.add_argument('output',
        help='Output .fasta file.')
    filt = parser.add_mutually_exclusive_group(required=True)
    filt.add_argument('--longest', default=None, type=float,
        help='Percentage of longest reads to partition.')
    filt.add_argument('--bases', default=None, type=int,
        help='Maximum number of bases (subject to at least one read.)')
    parser.add_argument('--others', default=None,
        help='Write all other reads to file.')
    args = parser.parse_args()

    sys.stderr.write('Loading reads...\n')
    record_dict = dict()
    for rec in pysam.FastxFile(args.input):
        record_dict[rec.name] = (rec.sequence, rec.quality, rec.comment)

    ids = list(record_dict.keys())
    lengths = np.fromiter(
        (len(record_dict[i][0]) for i in ids),
        dtype=int, count=len(ids)
    )
    sys.stderr.write('Sorting reads...\n')
    if args.bases is None:
        # partial sort will do fine here
        max_reads = int(len(ids) * (args.longest / 100))
        longest = np.argpartition(lengths, -max_reads)[-max_reads:]
    else:
        # need a full sort
        order = np.argsort(lengths)[::-1]
        cumsum = 0
        last = 1
        for i, j in enumerate(order, 1):
            cumsum += lengths[j]
            if cumsum > args.bases:
                break
            last = i
        longest = order[:last]

    with FastxWrite(args.output, 'w') as out:
        for i in longest:
            name = ids[i]
            rec = record_dict[ids[i]]
            out.write(name, *rec)

    if args.others is not None:
        if args.bases is None:
            others = np.argpartition(lengths, -max_reads)[:-max_reads]
        else:
            others = order[last:]

        with FastxWrite(args.others, 'w') as out:
            for i in others:
                name = ids[i]
                rec = record_dict[ids[i]]
                out.write(name, *rec)


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


class SeqLen(argparse.Action):
    """Parse a sequence length from str such as 4.8mb or from fastx."""
    def __call__(self, parser, namespace, values, option_string=None):
        seq_len = None
        try:
            seq_len = int(values)
        except:
            suffixes = {
                'kb': 10 ** 3,
                'mb': 10 ** 6,
                'gb': 10 ** 9,
            }
            for suffix, multiplier in suffixes.items():
                if suffix in values.lower():
                    seq_len = int(multiplier * float(values.lower().replace(suffix, '')))
                    break
            if seq_len is None:
                try:
                    seq_len = sum(get_seq_lens(values))
                except:
                    raise ValueError('Could not get sequence length from {}.'.format(values))

        setattr(namespace, self.dest, seq_len)


def get_seq_lens(fastx):
    """Get sequence lengths from fastx file"""
    return [len(r.sequence) for r in pysam.FastxFile(fastx)]


def coverage_from_fastx():
    parser = argparse.ArgumentParser(
        prog='coverage_from_fastx',
        description='Estimate coverage from summed basecall and reference lengths',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('basecalls', help='input fastx file.')
    parser.add_argument('ref_len', action=SeqLen,
                        help='reference length (e.g. 4.8kb/mb/gb) or reference fastx from which to calculate length.')
    parser.add_argument('--coverage', type=int, help='Calculate fraction of reads required for this coverage.')
    parser.add_argument('--longest', action='store_true', default=False,
                        help='Use the longest reads when calculating fraction reads required for a given coverage.')

    args = parser.parse_args()

    basecall_lens = get_seq_lens(args.basecalls)
    total_basecalls_len = sum(basecall_lens)

    print('Total length of basecalls: {}'.format(total_basecalls_len))
    print('Total ref length: {}'.format(args.ref_len))
    print('Estimated coverage: {:.2f}'.format(total_basecalls_len / args.ref_len))

    if args.coverage is not None:
        if args.longest:
            basecall_lens.sort(reverse=True)
            cum_len = np.cumsum(basecall_lens)
            required_len = args.coverage * args.ref_len
            n_required = np.searchsorted(cum_len, required_len, side='right')
            msg = '% of longest reads required to achieve {}X coverage: {}'
            print(msg.format(args.coverage, 100 * n_required / len(cum_len)))
        else:
            msg = '% of randomly-selected reads required to achieve {}X coverage: {}'
            print(msg.format(args.coverage, 100 * args.coverage  * args.ref_len / total_basecalls_len))


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
    :returns: { str chrom: `ncls.NCLS` obj }.
    """
    tmptrees = defaultdict(list)
    trees = {}
    for chrom, start, stop in yield_from_bed(path_to_bed):
        tmptrees[chrom].append([start, stop])
    for chrom, d in tmptrees.items():
        a = np.array(d)
        trees[chrom] = NCLS(a[:, 0], a[:, 1], np.arange(a.shape[0]))
    return trees


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
        "ref": read.reference_name,
        "aligned_ref_len": read.reference_length,
        "ref_coverage": 100*float(read.reference_length) / read.header.lengths[read.reference_id],
        "mapq": read.mapping_quality,
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
        "ref": read.reference_name,
        "aligned_ref_len": aligned_ref_len,
        "ref_coverage": 100 * float(aligned_ref_len) / read.header.lengths[read.reference_id],
        "mapq": read.mapping_quality,
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
        mean_q = -10 * np.log10(np.mean(QSCORES_TO_PROBS[r.query_qualities]))
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
