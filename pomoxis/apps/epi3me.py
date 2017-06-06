import argparse
import os
import sys
import shutil
import asyncio
import aiozmq
import zmq
import uuid
import struct
from functools import partial
import socket

from concurrent.futures import ProcessPoolExecutor as PoolExecutor
from multiprocessing import freeze_support, Process

from pomoxis import set_wakeup
from pomoxis.watcher import watch_path
from pomoxis.pyscrap import pyscrap

import logging


def basecall_file(fname, event_detect=True):
    score, basecall = 0, ''
    try:
        results = pyscrap.basecall_file(fname, event_detect=event_detect)
        if results is not None:
            score, basecall = results
    except:
        pass
    return score, basecall


@asyncio.coroutine
def main_router(addr, name=None):
    if name is None:
        name = uuid.uuid4()
    logging.info("Starting client {} connected to {}".format(name, addr))
    router = yield from aiozmq.create_zmq_stream(
        zmq.ROUTER
    )
    yield from router.transport.connect(addr)

    event_loop = asyncio.get_event_loop()
    executor = PoolExecutor(max_workers=1)
    event_loop.set_default_executor(executor)

    while True: 
        try:
            data = yield from router.read()
        except aiozmq.ZmqStreamClosed:
            logging.info('Client {} got stream closed'.format(name))
            break
        if data[1] == b'FINISHED':
            logging.info('Client {} got finished signal'.format(name))
            break

        #TODO: restructure so we can multiprocess
        fname = data[1].decode()
        score, basecall = yield from event_loop.run_in_executor(None, basecall_file, fname)

        data = data + [name.bytes, struct.pack('f', score), basecall.encode()]
        router.write(data)
    router.close()


@asyncio.coroutine
def main_dealer(path, outpath, output=None, port='*'):
    logging.info("Starting server...")
    ip = socket.gethostbyname(socket.gethostname())
    if ip == '127.0.0.1':
        logging.warn("ip found as 127.0.0.1, remote connections won't work.")
    dealer = yield from aiozmq.create_zmq_stream(
        zmq.DEALER,
        bind='tcp://{}:{}'.format(ip, port))
    addr = list(dealer.transport.bindings())[0]
    logging.info("Server started on {}.".format(addr))
    
    path = os.path.abspath(path)
    logging.info("Watching files at {}".format(path))

    @asyncio.coroutine
    def _callback(fname):
        dealer.write((fname.encode(),))
        
    event_loop = asyncio.get_event_loop()
    event_loop.create_task(watch_path(path, _callback, recursive=True))

    file_counter = 0
    call_counter = 0
    output_handle = None
    calls_per_file = 100000

    fail_folder = os.path.join(outpath, '_FAILED')
    if os.path.exists(fail_folder):
        os.rmdir(fail_folder)
    os.mkdir(fail_folder)
    pass_folder = os.path.join(outpath, '_PASSED')
    if os.path.exists(pass_folder):
        os.rmdir(pass_folder)
    os.mkdir(pass_folder)

    if output is not None:
        output, ext = os.path.splitext(output)
    def _start_new_file(fh, counter):
        if output is None:
            return sys.stdout
        if fh is not None:
            logging.info("Closing output file: {}".format(fh))
            fh.close()
        fname = '{}_{}{}'.format(output, counter, ext)
        logging.info("Starting new output file: {}".format(fname))
        return open(fname, 'w')

    try:
        while True:
            answer = yield from dealer.read()
            fname = answer[0].decode()
            client = uuid.UUID(bytes=answer[1])
            score = struct.unpack('f', answer[2])[0]
            basecall = answer[3].decode()
            if basecall != '':
                try:
                    shutil.move(fname, pass_folder)
                except:
                    logging.info("Could not move passed file {}, ...Penny Lane.".format(fname))
                if call_counter % calls_per_file == 0:
                    file_counter += 1
                    output_handle = _start_new_file(output_handle, file_counter)
                fasta_name = os.path.basename(fname)
                output_handle.write(">{} score={}\n{}\n".format(
                    fasta_name, score, basecall)
                )
                call_counter += 1
            else:
                logging.info("Basecall failed for {}.".format(fname))
                try:
                    shutil.move(fname, fail_folder)
                except:
                    logging.info("Could not move failed file {} to fail folder, ...Penny Lane.".format(fname))
                
    except KeyboardInterrupt:
        logging.info("Closing down server...")
        if output is not None:
            output_handle.close()

    dealer.close()


def run_dealer(args):
    set_wakeup()
    asyncio.get_event_loop().run_until_complete(main_dealer(args.path, args.outpath, output=args.output, port=args.port))

def run_router(args):
    set_wakeup()
    asyncio.get_event_loop().run_until_complete(main_router(args.addr))


def main():
    logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    parser = argparse.ArgumentParser('Epi3Me -- Distributed basecalling platform.')
    subparsers = parser.add_subparsers(title='subcommands', description='valid commands', help='additional help', dest='command')
    subparsers.required = True

    dealer = subparsers.add_parser('dealer', help='Create a ZMQ dealer (server).')
    dealer.set_defaults(func=run_dealer)
    dealer.add_argument('path', help='Path to watch for new files.')
    dealer.add_argument('outpath', help='Path to move finished files.')
    dealer.add_argument('--port', type=int, default='*', help='port on which to run.')
    dealer.add_argument('--output', type=str, default=None, help='output fasta file.')
    
    router = subparsers.add_parser('router', help='Create a ZMQ router (client).')    
    router.set_defaults(func=run_router)
    router.add_argument('--addr', help='Address to use, should include port.')


    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
