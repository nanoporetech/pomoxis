import argparse
from concurrent.futures import ThreadPoolExecutor
import functools
import itertools
import logging

import numpy as np
from pysam import FastaFile

import scrappy

from pomoxis.common.bio import reverse_complement, shotgun_library


def worker(args, noise=None):
    model='rgrgr_r94'
    seq, ref, start, end, strand = args
    squiggle = scrappy.sequence_to_squiggle(
        seq, rescale=True).data(as_numpy=True, sloika=False)

    n = 1/np.sqrt(2)
    raw_data = np.concatenate([
        np.random.laplace(mean, n*stdv, int(dwell))
        for mean, stdv, dwell in squiggle
    ])
    if noise is not None:
        raw_data += np.random.normal(scale=noise, size=len(raw_data))

    raw = scrappy.RawTable(raw_data)
    raw.scale()
    post = scrappy.calc_post(raw, model, log=True)
    call, score, _ = scrappy.decode_post(post, model)
    return '>call_{}:{}-{}({}) seq_len={} call_len={} score={}\n{}'.format(
        ref, start, end, strand, len(seq), len(call), score, call)


def main():
    logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    parser = argparse.ArgumentParser(description='Simulate basecalls with scrappy.')

    parser.add_argument('fasta', help='Source sequence file.')
    parser.add_argument('ncalls', type=int, help='Number of basecalls to produce.')
    parser.add_argument('--mu', type=float, default=8000, help='mean fragment length.')
    parser.add_argument('--sigma', type=float, default=1000, help='stdv fragment length.')
    parser.add_argument('--noise', type=float, default=0.06, help='Additional Gaussian noise on signal.')
    parser.add_argument('--threads', type=int, default=None, help='number of worker threads.')

    args = parser.parse_args()

    regions = itertools.islice(shotgun_library(
        args.fasta, args.mu, args.sigma, direction=(1,-1)
    ), args.ncalls)

    _worker = functools.partial(worker, noise=args.noise)
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        for fasta in executor.map(_worker, regions):
            print(fasta)


if __name__ == '__main__':
    main()
