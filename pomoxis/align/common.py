import argparse
import asyncio
from collections import namedtuple
from multiprocessing import freeze_support
import subprocess

from aiozmq import rpc

from pomoxis.align.minimap import MiniMapServe

from pomoxis import set_wakeup
from pomoxis.common import util

import logging
logger = logging.getLogger(__name__)


@asyncio.coroutine
def align_server(index, port, aligner, opts=''):
    bind = 'tcp://127.0.0.1:{}'.format(port)
    if aligner == 'minimap':
        server = yield from rpc.serve_rpc(
            MiniMapServe(index[0], map_opts=opts), bind=bind
        )
    else:
        raise ValueError("Unknown aligner '{}'.".format(aligner))
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
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        result = loop.run_until_complete(self._run_align(sequence))
        loop.close()
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
            args.index, args.port, args.aligner, opts=args.opts
        )
        logger.info('Alignment server running, awaiting requests...')
        yield from server.wait_closed()

    loop = asyncio.get_event_loop()
    try:
        loop.run_until_complete(_run())
    except KeyboardInterrupt:
        logger.info('Shutting down server...')
        loop.run_until_complete(clean_up())
        logger.info('Server shut down.')


def send(args):
    client = AlignClient(args.port)
    for seq in args.sequences:
        print(client.align(seq))


def to_dict(string):
   """Convert a string whitespace separate key:value pairs to a dict."""
   return {k:v for k,v in (x.split(":") for x in string.split())}


def get_parser():
    parser = argparse.ArgumentParser('Alignment server/client.')
    subparsers = parser.add_subparsers(title='subcommands', description='valid commands', help='additional help', dest='command')
    subparsers.required = True

    sparser = subparsers.add_parser('server', help='Launch alignment server.')
    sparser.set_defaults(func=serve)
    sparser.add_argument('port', type=int, help='Port on which to serve.')
    sparser.add_argument('index', nargs='+', help='Filename path prefix for index files.')
    sparser.add_argument('--aligner', choices=('minimap'), default='minimap', help='Choice of aligner.')
    sparser.add_argument('--opts', default=dict(), type=to_dict, help="Alignment options as 'key:value key:value ...'.")

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


