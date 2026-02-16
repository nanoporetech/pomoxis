import argparse
import itertools
import sys

import pysam


class FastxWrite:
    def __init__(self, fname, width=80):
        self.fname = fname
        self.width = width

    def __enter__(self):
        if self.fname == "-":
            self.fh = sys.stdout
        else:
            self.fh = open(self.fname, "w")
        return self

    def __exit__(self, *args, **kwargs):
        if self.fname != "-":
            self.fh.close()

    def write(self, name, seq, qual=None, comment=""):
        if qual:
            self.fh.write("@{} {}\n{}\n+\n{}\n".format(name, comment, seq, qual))
        else:
            self.fh.write(">{} {}\n".format(name, comment))
            for _, chunk in enumerate(chunks(seq, self.width)):
                self.fh.write("{}\n".format("".join(chunk)))


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


def split_fastx(fname, output, chunksize=10000):
    """Split records in a fasta/q into fixed lengths.

    :param fname: input filename.
    :param output: output filename.
    :param chunksize: (maximum) length of output records.
    """
    with pysam.FastxFile(fname, persist=False) as f_in, FastxWrite(output) as f_out:
        for rec in f_in:
            name = rec.name
            seq = rec.sequence
            qual = rec.quality
            if rec.comment is None:
                comment = f"chunk_length={chunksize}"
            else:
                comment = f"{rec.comment} chunk_length={chunksize}"

            if qual is None:
                for i, s in enumerate(chunks(seq, chunksize)):
                    f_out.write(f"{name}_chunk{i}", s, comment=comment)
            else:
                for i, (s, q) in enumerate(
                    zip(chunks(seq, chunksize), chunks(qual, chunksize))
                ):
                    f_out.write(f"{name}_chunk{i}", "".join(s), "".join(q), comment)


def main():
    """Split records in a fasta/q file into chunks of a maximum size."""
    parser = argparse.ArgumentParser(
        prog="split_fastx",
        description="Split records in a fasta/q file into chunks of a maximum size.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("input", help="Input fastax/q file.")
    parser.add_argument("output", help="Output fastax/q file.")
    parser.add_argument("chunksize", type=int, help="Maximum size of output sequences.")
    args = parser.parse_args()
    split_fastx(args.input, args.output, args.chunksize)


if __name__ == "__main__":
    main()
