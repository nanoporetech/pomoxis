import asyncio
import subprocess

from aiozmq import rpc

from pomoxis import get_prog_path, run_prog
from pomoxis.common import util

import logging
logger = logging.getLogger(__name__)


class BwaServe(rpc.AttrHandler):

    def __init__(self, index, *args, bwa_cmd=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.logger = logging.getLogger('BwaServe')
        if isinstance(index, (str, bytes)):
            index = [index]
        self.index = index
        self.bwa_cmd = 'mem -x ont2d'
        if bwa_cmd is not None:
            self.bwa_cmd  = bwa_cmd
        self.bwa_cmd = self.bwa_cmd.split()
        self.bwa = get_prog_path('bwa')
        self.logger.info("Found bwa at {}".format(self.bwa))

        for ind in self.index:
            try:
                run_prog(self.bwa, ['shm', ind])
            except Exception as e:
                logger.debug(e)
                raise RuntimeError('Cannot load bwa index "{}" into shared memory.'.format(ind))
        self.logger.info('BWA service started.')


    @rpc.method
    @asyncio.coroutine
    def align(self, sequence):
        self.logger.debug("Aligning sequence of length {}.".format(len(sequence)))
        
        if isinstance(sequence, bytes):
            sequence = sequence.decode('utf-8')
        sequence = '>seq\n{}\n'.format(sequence)
        
        #TODO move this to a process pool
        results = []
        returncode = 0
        for ind in self.index:
            proc = subprocess.Popen(
                [self.bwa, 'mem', '-x', 'ont2d', '-A1', '-B2', '-O2', '-E1', ind, '-'],
                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = proc.communicate(sequence.encode('utf-8'))
            if proc.returncode != 0:
                logger.warn(stderr)
            returncode = max(returncode, proc.returncode)
            results.append(stdout.decode('utf-8'))
        return ''.join(results), returncode

    @rpc.method
    def clean_shm(self):
        run_prog(self.bwa, ['shm', '-d'])


@asyncio.coroutine
def align_server(index, port):
    server = yield from rpc.serve_rpc(
        BwaServe(index), bind='tcp://127.0.0.1:{}'.format(port)
    )
    return server

@asyncio.coroutine
def align_client(port):
    client = yield from rpc.connect_rpc(
        connect='tcp://127.0.0.1:{}'.format(port),
    )
    return client
