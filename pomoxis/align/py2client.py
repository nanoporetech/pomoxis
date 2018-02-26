from __future__ import print_function
import random
import struct
import sys
import time
import os

import zmq

from msgpack import ExtType, packb, unpackb

class AlignClient(object):
    # A synchronous Python2 alignment client 

    REQ_PREFIX = struct.Struct('=HH')
    REQ_SUFFIX = struct.Struct('=Ld')
    RESP = struct.Struct('=HHLd?')

    def __init__(self, port):
        self.port = port
        self.ctx = zmq.Context()
        self.sock = self.ctx.socket(zmq.DEALER)

        self.sock.connect('tcp://127.0.0.1:{}'.format(self.port))

        self.counter = 0
        self.prefix = self.REQ_PREFIX.pack(os.getpid() % 0x10000,
                                           random.randrange(0x10000))

    def packb(self, data):
        return packb(data, encoding='utf-8', use_bin_type=True)

    def unpackb(self, packed):
        return unpackb(packed, use_list=False, encoding='utf-8')

    def _new_id(self):
        self.counter += 1
        if self.counter > 0xffffffff:
            self.counter = 0
        return (self.prefix + self.REQ_SUFFIX.pack(self.counter, time.time()),
                self.counter)

    def align(self, sequence):
        header, req_id = self._new_id()
        bname = 'align'.encode('utf-8')
        bargs = self.packb([sequence.decode()])
        bkwargs = self.packb(dict())
        msg = [header, bname, bargs, bkwargs]
        self.sock.send_multipart(msg)
        while True:
            try:
                data = self.sock.recv_multipart(zmq.NOBLOCK)
            except zmq.ZMQError as e:
                time.sleep(0.1)
            else:
                header, banswer = data
                pid, rnd, res_req_id, timestamp, is_error = self.RESP.unpack(header)
                if res_req_id != req_id:
                    raise ValueError('Received response for request {}, but send {}.'.format(res_res_id, req_id))
                answer = self.unpackb(banswer)
                return answer


if __name__ == '__main__':
    port, seq = sys.argv[1:3]
    client = AlignClient(port)
    alignments = client.align(seq)
    print(alignments)
            
