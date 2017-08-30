import itertools
import os
import shutil
import sys

from Bio import SeqIO
from pysam import FastxFile


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
    with open(output,'wb') as wfd:
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
                        fout.write('@{} {}\n{}\n+{}\n'.format(
                            chunk_name, comment, ''.join(s), ''.join(q)))


def split_fastx_cmdline():
    """Split records in a fasta/q file into chunks of a maximum size."""
    fname, output, chunksize = sys.argv[1:]
    split_fastx(fname, output, int(chunksize))


def fast_convert(convert=None, mock_q=10):
    """Convert between fasta<->fastq.

    :param convert: conversion code: from->to, e.g. 'aq'.
    :param mockq: mock quality value, valid for convert=aq.
    """

    if convert is None:
        parser = argparse.ArgumentParser(description='fast_convert -- Convert between fasta<->fastq.')
        parser.add_argument('convert', choices=['qq', 'aa', 'aq', 'qa'],
            help='Conversion code: from->to.')
        parser.add_argument('--mock_q', default=mock_q, type=int,
            help='Mock quality value, valid for convert=aq.')
        args = parser.parse_args()
        convert, mock_q = args.convert, args.mock_q

    in_fmt = 'fastq'
    out_fmt= 'fasta'
    flag = False
    if convert == 'qq':
        out_fmt = 'fastq'
    elif convert == 'aa':
        in_fmt = 'fasta'
    elif convert == 'aq':
        in_fmt = 'fasta'
        out_fmt = 'fastq'
        flag = True
    elif convert == 'qa':
        pass # default
    else:
        raise ValueError("convert must be 'qq', 'aq', or aa'\n")
    
    if flag:
        def fq_gen(io):
            for rec in SeqIO.parse(io, in_fmt):
                rec.letter_annotations["phred_quality"] = [mock_q] * len(rec)
                yield rec
        sys.stderr.write('Mocking fastq from fasta\n')
        SeqIO.write(fq_gen(sys.stdin), sys.stdout, out_fmt)
    else:
        SeqIO.convert(sys.stdin, in_fmt, sys.stdout, out_fmt)

