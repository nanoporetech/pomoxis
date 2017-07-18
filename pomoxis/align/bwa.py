import argparse
import asyncio
from multiprocessing import freeze_support
import subprocess

from aiozmq import rpc

try:
    from bwapy import BwaAligner
except ImportError:
    BwaAligner = None


from pomoxis import get_prog_path, run_prog, set_wakeup
from pomoxis.common import util

import logging
logger = logging.getLogger(__name__)


class BwaServe(rpc.AttrHandler):

    def __init__(self, index, *args, clean=False, bwa_opts='-x ont2d', **kwargs):
        """bwa mem alignment server implementation using shared memory and
        subprocess calls.

        :param index: bwa index base path, or list thereof.
        :param clean: clean-up shared memory on exit.
        :param bwa_opts: command line options for bwa mem.

        """
        DeprecationWarning('The class BwapyServe should be used in preference.')
        super().__init__(*args, **kwargs)
        self.logger = logging.getLogger('BwaServe')
        if isinstance(index, (str, bytes)):
            index = [index]
        self.index = index
        self.clean = clean
        self.bwa_opts = bwa_opts.split()
        self.bwa = get_prog_path('bwa')
        self.logger.info("Found bwa at {}".format(self.bwa))

        self.loaded = False
        self._load_index()
        self.logger.info('bwa service started.')

    def _load_index(self):
        """Load indices into shared memory."""
        if not self.loaded:
            for ind in self.index:
                try:
                    run_prog(self.bwa, ['shm', ind])
                except Exception as e:
                    logger.debug(e)
                    raise RuntimeError('Cannot load bwa index "{}" into shared memory.'.format(ind))
        self.loaded = True

    def _clean_index(self):
        self.logger.info('Cleaning shared memory.')
        run_prog(self.bwa, ['shm', '-d'])
        self.loaded = False

    def __del__(self):
        if self.clean:
            self._clean_index()

    @rpc.method
    def clean_index(self):
        """Clean bwa shared memory, note that this clears any and all
        indices in memory not just the ones placed by this class.

        """
        self._clean_index()

    @rpc.method
    @asyncio.coroutine
    def align(self, sequence):
        """Align a base sequence.

        :param sequence: sequence to align.

        :returns: the output of bwa mem call.
        """
        self.logger.debug("Aligning sequence of length {}.".format(len(sequence)))
        if not self.loaded:
            self._load_index()       
 
        if isinstance(sequence, bytes):
            sequence = sequence.decode('utf-8')
        sequence = '>seq\n{}\n'.format(sequence)
        
        #TODO: move this to a thread pool
        results = []
        returncode = 0
        for ind in self.index:
            proc = subprocess.Popen(
                [self.bwa, 'mem'] + self.bwa_opts + [ind, '-'],
                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = proc.communicate(sequence.encode('utf-8'))
            if proc.returncode != 0:
                logger.warn(stderr)
            returncode = max(returncode, proc.returncode)
            results.append(stdout.decode('utf-8'))
        return ''.join(results), returncode


class BwapyServe(rpc.AttrHandler):

    def __init__(self, index, *args, clean=False, bwa_opts='-x ont2d', **kwargs):
        """bwa mem alignment server implementation using python binding.

        :param index: bwa index base path, or list thereof.
        :param clean: clean-up shared memory on exit.
        :param bwa_opts: command line options for bwa mem.

        """
        super().__init__(*args, **kwargs)
        self.logger = logging.getLogger('BwaServe')
        self.index = index
        self.bwa_opts = bwa_opts

        self.aligner = None
        if BwaAligner is None:
            raise ImportError(
                '{} requires BwaAligner which could not be imported.'.format(
                self.__class__.__name__
            ))
        self.aligner = BwaAligner(self.index, options=self.bwa_opts)
        self.logger.info('bwa service started.')

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
            self.aligner = BwaAligner(self.index, options=self.bwa_opts)
        self.logger.debug("Aligning sequence of length {}.".format(len(sequence)))
        return self.aligner.align_seq(sequence)


@asyncio.coroutine
def align_server(index, port, shm=False, clean=False):
    if shm:
        server = yield from rpc.serve_rpc(
            BwaServe(index, clean=clean), bind='tcp://127.0.0.1:{}'.format(port)
        )
    else:
        server = yield from rpc.serve_rpc(
            BwapyServe(index[0]), bind='tcp://127.0.0.1:{}'.format(port)
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

    def align(self, sequence, loop=None):
        """Align a given sequence."""
        if loop is None:
            loop = asyncio.get_event_loop()
        result = loop.run_until_complete(self._run_align(sequence))
        return result


def serve(args):
    @asyncio.coroutine
    def clean_up():
        # Doesn't seem right to have a client clean up the server
        logger.info('Creating cleanup client')
        client = yield from align_client(args.port)
        yield from client.call.clean_index()

    @asyncio.coroutine
    def _run():
        set_wakeup()
        server = yield from align_server(
            args.bwa_index, args.port, shm=args.shm, clean=not args.noclean
        )
        logger.info('Alignment server running, awaiting requests...')
        yield from server.wait_closed()

    loop = asyncio.get_event_loop()
    try:
        loop.run_until_complete(_run())
    except KeyboardInterrupt:
        logger.info('Shutting down server...')
        if not args.noclean:
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
    sparser.add_argument('--shm', action='store_true', help='Use deprecated shared memory implementation server.')
    sparser.add_argument('--noclean', action='store_true', default=False, help='Clean bwa shm when done.')

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


