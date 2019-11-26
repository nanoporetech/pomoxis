import argparse
from collections import namedtuple, defaultdict
import itertools
import shutil
import sys

from Bio import SeqIO
import intervaltree
import numpy as np
from pysam import FastxFile

Region = namedtuple('Region', 'ref_name start end')
AlignPos = namedtuple('AlignPos', ('qpos', 'qbase', 'rpos', 'rbase'))


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
    with open(output, 'w') as fout:
        with FastxFile(fname, persist=False) as fin:
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
                        fout.write(">{} {}\n{}\n".format(
                            chunk_name, comment, ''.join(s)))
                else:
                    for i, (s, q) in enumerate(zip(chunks(seq, chunksize), chunks(qual, chunksize))):
                        chunk_name = '{}_chunk{}'.format(name, i)
                        fout.write('@{} {}\n{}\n+\n{}\n'.format(
                            chunk_name, comment, ''.join(s), ''.join(q)))


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

    if qflag:
        def fq_gen(io):
            for rec in SeqIO.parse(io, in_fmt):
                rec.letter_annotations["phred_quality"] = [args.mock_q] * len(rec)
                yield rec
        sys.stderr.write('Creating/ignoring quality information in input.\n')
        SeqIO.write(fq_gen(sys.stdin), sys.stdout, out_fmt)
    else:
        SeqIO.convert(sys.stdin, in_fmt, sys.stdout, out_fmt)


def extract_long_reads():
    """Filter fastq to longest reads."""

    parser = argparse.ArgumentParser(description='Extract longest reads from a fastq.')
    parser.add_argument('input',
        help='Input .fastq file.')
    parser.add_argument('output',
        help='Output .fastq file.')
    filt = parser.add_mutually_exclusive_group(required=True)
    filt.add_argument('--longest', default=None, type=float,
        help='Percentage of longest reads to partition.')
    filt.add_argument('--bases', default=None, type=int,
        help='Maximum number of bases (subject to at least one read.)')
    parser.add_argument('--others', default=None,
        help='Write all other reads to file.')
    args = parser.parse_args()

    sys.stderr.write('Loading reads...\n')
    fmt = 'fastq'
    try:
        record_dict = SeqIO.index(args.input, fmt)
    except:
        fmt = 'fasta'
        record_dict = SeqIO.index(args.input, fmt)
    sys.stderr.write('Format is {}.\n'.format(fmt))

    ids = list(record_dict.keys())
    lengths = np.fromiter(
        (len(record_dict[i]) for i in ids),
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
        for i, j in enumerate(np.argsort(lengths), 1):
            cumsum += lengths[j]
            if cumsum > args.bases:
                break
            last = i
        longest = order[:last]

    SeqIO.write(
        (record_dict[ids[i]] for i in longest),
        args.output, fmt
    )

    if args.others is not None:
        longest = set(longest)
        SeqIO.write(
            (record_dict[ids[i]] for i in range(len(ids)) if i not in longest),
            args.others, fmt
        )


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
    return [len(r.sequence) for r in FastxFile(fastx)]


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
    :returns: { str chrom: `intervaltree.IntervalTree` obj }.
    """
    trees = defaultdict(intervaltree.IntervalTree)
    for chrom, start, stop in yield_from_bed(path_to_bed):
        trees[chrom].add(intervaltree.Interval(begin=start, end=stop))
    return trees
