import argparse
import asyncio
from multiprocessing import freeze_support
import subprocess

from aiozmq import rpc

from pomoxis import get_prog_path, run_prog, set_wakeup
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
        self.logger.info('Cleaning shared memory.')
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


class AlignClient(object):
    def __init__(self, port):
        """A synchronous alignment client.

        :param port: RPC server port.
        """
        self.port = port

    @asyncio.coroutine
    def _run_align(self, sequence):
        client = yield from align_client(self.port)
        results = yield from client.call.align(sequence)
        return results

    def align(self, sequence):
        """Align a given sequence."""
        loop = asyncio.get_event_loop()
        result = loop.run_until_complete(self._run_align(sequence))
        return result


def serve(args):
    loop = asyncio.get_event_loop()
    set_wakeup()
    server = loop.create_task(align_server(
        args.bwa_index, args.port
    ))

    @asyncio.coroutine
    def clean_up():
        client = yield from align_client(args.port)
        yield from client.call.clean_shm()
    
    logger.info('Alignment server running, awaiting requests...')
    try:
        loop.run_forever()
    except KeyboardInterrupt:
        logger.info('Shutting down server...')
        if args.clean:
            loop.run_until_complete(clean_up())
        logger.info('Server shut down.')

def send(args):
    client = AlignClient(args.port)
    for seq in args.sequences:
        print(client.align(seq))


def get_parser():
    parser = argparse.ArgumentParser('BWA alignment server/client.')
    subparsers = parser.add_subparsers(title='subcommands', description='valid commands', help='additional help', dest='command')
    subparsers.required = True

    sparser = subparsers.add_parser('server', help='Launch alignment server.')
    sparser.set_defaults(func=serve)
    sparser.add_argument('port', type=int, help='Port on which to serve.')
    sparser.add_argument('bwa_index', nargs='+', help='Filename path prefix for BWA index files.')
    sparser.add_argument('--clean', action='store_true', default=False, help='Clean bwa shm when done.')


    sparser = subparsers.add_parser('client', help='Test client.')
    sparser.set_defaults(func=send)
    sparser.add_argument('port', type=int, help='Port on which to serve.')
    sparser.add_argument('sequences', nargs='+', help='Base sequences to align.')

    return parser


def main():
    freeze_support()
    logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    args = get_parser().parse_args()
    args.func(args)


