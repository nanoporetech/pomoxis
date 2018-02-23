import asyncio
from collections import namedtuple

from aiozmq import rpc

import mappy

import logging
logger = logging.getLogger(__name__)


class MiniMapServe(rpc.AttrHandler):

    def __init__(self, index, *args, map_opts={'preset':'map-ont'}, **kwargs):
        """bwa mem alignment server implementation using python binding.

        :param index: bwa index base path, or list thereof.
        :param map_opts: command line options for minimap

        """
        super().__init__(*args, **kwargs)
        self.logger = logging.getLogger('MiniMapServe')
        self.index = index
        self.map_opts = map_opts

        self.aligner = None
        self.aligner = mappy.Aligner(self.index, **map_opts)
        self.logger.info('Minimap service started.')

    def _clean_index(self):
        self.logger.info('Cleaning alignment proxy.')
        self.aligner = None

    def __del__(self):
            self._clean_index()

    @rpc.method
    def clean_index(self):
        """Destroy the aligner object, which will cleanup the index in memory."""
        self._clean_index()

    @rpc.method
    @asyncio.coroutine
    def align(self, sequence):
        """Align a base sequence.

        :param sequence: sequence to align.

        :returns: the output of bwa mem call.
        """
        if self.aligner is None:
            self.aligner = mappy.Aligner(self.index, **map_opts)
        self.logger.debug("Aligning sequence of length {}.".format(len(sequence)))
        # create named tuples from the results
        results = []
        for r in self.aligner.map(sequence):
             results.append((
                 r.q_st, r.q_en, '+' if r.strand > 0 else '-',
                 r.ctg, r.ctg_len, r.r_st, r.r_en,
                 r.mlen, r.blen, r.mapq, r.is_primary,
                 r.cigar_str
             ))
        self.logger.info("Aligned sequence of {} bases with {} hits.".format(
            len(sequence), len(results)
        ))
        return results

