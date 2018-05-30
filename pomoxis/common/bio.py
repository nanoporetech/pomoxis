import numpy as np
from pysam import FastaFile

"""Bioinformatics helpers."""

comp = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'X': 'X', 'N': 'N',
    'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'x': 'x', 'n': 'n',
    #'-': '-'
}
comp_trans = str.maketrans(''.join(comp.keys()), ''.join(comp.values()))


def reverse_complement(seq):
    """Reverse complement sequence.

    :param: input sequence string.

    :returns: reverse-complemented string.
    """
    return seq.translate(comp_trans)[::-1]


def shotgun_library(fasta_file, mu, sigma, direction=(1,-1)):
    """Generate random fragment sequences of a given input sequence

    :param seq: input sequence.
    :param mu: mean fragment length.
    :param sigma: stdv of fragment length.
    :param direction: tuple represention direction of output sequences with
        respect to the input sequence.

    :yields: sequence fragments.

    .. note:: Could be made more efficient using buffers for random samples
        and handling cases separately.
    """
    fasta = FastaFile(fasta_file)
    seq_lens = [fasta.get_reference_length(x) for x in fasta.references]
    total_len = sum(seq_lens)
    seq_probs = [x / total_len for x in seq_lens]
    # FastaFile.fetch is proper slow, just read everything
    refs = fasta.references
    fasta = {k:fasta.fetch(k) for k in refs}

    def random_buffer(probs, size=10000):
        while True:
            buf = []
            for x, n in zip(range(len(probs)), np.random.multinomial(size, probs)):
                buf.extend([x]*n)
            np.random.shuffle(buf)
            for x in buf:
                yield x
    seq_chooser = random_buffer(seq_probs)
   
    # parameters for lognormal
    mean = np.log(mu / np.sqrt(1 + sigma**2 / mu**2))
    stdv = np.sqrt(np.log(1 + sigma**2 / mu**2))
    
    while True:
        # choose a seq based on length
        seq_i = next(seq_chooser)
        seq = fasta[refs[seq_i]]
        seq_len = seq_lens[seq_i]

        start = np.random.randint(0, seq_len)
        frag_length = int(np.random.lognormal(mean, stdv))
        move = np.random.choice(direction)
        end = max(0, start + move*frag_length)
        start, end = sorted([start, end])

        if end - start < 2:
            # Expand a bit to ensure we grab at least one base.
            start = max(0, start - 1)
            end += 1

        frag_seq = seq[start:end]
        if move == -1:
            frag_seq = reverse_complement(frag_seq)
        yield frag_seq, refs[seq_i], start, end, '+' if move == 1 else '-'

