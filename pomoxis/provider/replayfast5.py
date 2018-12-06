import argparse
import sys
import os
import shutil
import copy
import platform
from collections import Counter
from multiprocessing import freeze_support
from timeit import default_timer as now

import asyncio
from aiozmq import rpc
from concurrent.futures import ProcessPoolExecutor as PoolExecutor
from functools import partial

import numpy as np
import msgpack

from fast5_research.fast5_bulk import BulkFast5

from pomoxis import set_wakeup

import logging
logger = logging.getLogger(__name__)


class Fast5Data(np.ndarray):

    def __new__(cls, input_array, info=None, start=None, end=None, sample_rate=6000):
        """Helper class for transmitting raw and event data over RPC requests.

        :param input_array: numpy array to encode.
        :param info: metadata.
        :param start: start time of data.
        :param end: end time of data.
        :param sample_rate: sampling rate of source data.
        """
        obj = np.asarray(input_array).view(cls)
        obj.info = info
        obj.start = start
        obj.end = end
        obj.sample_rate = sample_rate
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        self.info = getattr(obj, 'info', None)
        self.start = getattr(obj, 'start', None)
        self.end = getattr(obj, 'end', None)
        self.sample_rate = getattr(obj, 'sample_rate', None)

    @staticmethod
    def encode(obj):
        """Serialize data to send over RPC."""
        data = {
            b'shape': obj.shape,
            b'info': obj.info, 
            b'start': obj.start,
            b'end': obj.end,
            b'sample_rate': obj.sample_rate, 
            b'data': obj.tobytes(),
            b'type': obj.dtype.descr,
            b'kind': obj.dtype.kind
        }
        return data

    @staticmethod
    def decode(obj):
        """Deserialize data from RPC request."""
        if obj[b'kind'] == b'V':
            descr = [(x.decode('utf-8'), y.decode('utf-8')) for x, y in obj[b'type']]
        else:
            descr = obj[b'type'][0][1].decode('utf-8')
        data = Fast5Data(
            np.fromstring(
                obj[b'data'],
                dtype=np.dtype(descr)
            ).reshape(obj[b'shape']),
            info=obj[b'info'].decode('utf-8'),
            start=obj[b'start'],
            end=obj[b'end'],
            sample_rate=obj[b'sample_rate'],
        )
        return data


translation_table = {
    0: (Fast5Data,
        lambda value: msgpack.packb(value, default=value.encode),
        lambda binary: msgpack.unpackb(binary, object_hook=Fast5Data.decode)
    ),
}


class ReplayChannel(rpc.AttrHandler):

    def __init__(self, fast5, channel, *args, good_class='strand', time_warp=1, **kwargs):
        """An RPC service for replaying a channel from a .fast5 file.

        :param fast5: input filename.
        :param channel: channel to simulate.
        :param good_class: read classification name of desirable reads.
        :param time_warp: time multiplier for playback speed.

        ..note:: `args` and `kwargs` are passed to `aiozmq.rpc.AttrHandler`.
        
        """
        super().__init__(*args, **kwargs)
        self.fast5 = fast5
        self.channel = channel
        self.good_class = good_class
        self.time_warp = time_warp
        self.logger = logging.getLogger('Replay Channel {}'.format(channel.zfill(4)))

        with BulkFast5(self.fast5) as fh:
            self.sample_rate = fh.sample_rate
            self.time_base = 1.0 / self.sample_rate
            self.reads = [x for x in fh.get_reads(self.channel)]
            self.read_starts = np.ascontiguousarray([x['read_start'] for x in self.reads])

        self.reset_time()

    @rpc.method
    def reset_time(self):
        """Reset internal time to zero."""
        self.start_time = now()
        self.sample_offset = 0
        self._current_event = 0
        self._current_read = 0

    @property
    def current_sample(self):
        """Current simulated sample (the simulation time)."""
        return int((now() - self.start_time) * self.sample_rate * self.time_warp + self.sample_offset)


    @property
    def experiment_time(self):
        """Current simulated time."""
        return self.current_sample * self.time_base 

    @property
    def time_saved(self):
        """Difference between current real time and simulated time. (How far
        we've fast-forwarded the channel.)

        """
        return self.sample_offset * self.time_base

    @property
    def current_event(self):
        """Index of the current event."""
        prev = self._current_event
        with BulkFast5(self.fast5) as fh:
            event_path = fh.__event_data__.format(self.channel)
            new = np.searchsorted(fh[event_path][prev:]['start'], self.current_sample)
        self._current_event = int(new + prev)
        return self._current_event


    @property
    def current_read(self):
        """Index of the current read."""
        search_start = max(0, self._current_read - 2)
        self._current_read = search_start - 1 + np.searchsorted(self.read_starts[search_start:], self.current_sample)
        return self._current_read


    @property
    def good_reads_to_now(self):
        """Reads with classification equal to `good_class` until current time."""
        return [x for x in self.reads[:self._current_read + 1] if x['classification'] == self.good_class]
        
    @property
    def total_good_reads(self):
        """Total number of good reads seen up to current simulation time."""
        return len(self.good_reads_to_now)

    @property
    def cumulative_good_read_time_to_now(self):
        """Total duration of reads with classification equal to `good_class` until current time."""
        return sum(read['read_length'] for read in self.good_reads_to_now) * self.time_base


    @rpc.method
    def get_events(self, n_events=400):
        """Return events from the start of the current read.

        :param nevents: maximum number of events to return.

        :returns: Serialized event data, see :class:`Fast5Data`.
        """
        self.logger.debug("Request for events at {}".format(self.current_sample))
        read = self.reads[self.current_read]
        if read['classification'] != self.good_class:
            return None
        else:
            start = int(read['event_index_start'])
            end = min(self.current_event, start + n_events)
            self.logger.debug("Fetching events [{}, {}] for read {} starting at {}. Current sample is {}.".format(
                start, end,
                self.current_read, int(self.sample_offset + self.reads[self.current_read]['read_start']),
                self.current_sample
            ))
            with BulkFast5(self.fast5) as fh:
                events = fh.get_events(self.channel, event_indices=[start, end])
                return Fast5Data(
                    events, info=read['read_id'].decode('utf-8'),
                    start=int(self.sample_offset + read['read_start']),
                    end=int(self.sample_offset + read['read_start'] + read['read_length'])
                )

    @rpc.method
    def get_raw(self, seconds=1, delay=0.5):
        """Return events from the start of the current read.

        :param nevents: maximum number of events to return

        :returns: Serialized raw data, see :class:`Fast5Data`.
        """
        self.logger.debug("Request for raw at {}".format(self.current_sample))
        read = self.reads[self.current_read]
        if read['classification'] != self.good_class:
            return None
        else:
            start = int(read['read_start'] + self.sample_rate * delay)
            end = int(min(self.current_sample, start + self.sample_rate * seconds))
            if end <= start:
                return None
            self.logger.debug("Fetching raw [{}, {}] for read {} starting at {}. Current sample is {}.".format(
                start, end,
                self.current_read, int(self.sample_offset + self.reads[self.current_read]['read_start']),
                self.current_sample
            ))
            with BulkFast5(self.fast5) as fh:
                raw = fh.get_raw(self.channel, raw_indices=[start, end])
                return Fast5Data(
                    raw, info=read['read_id'].decode('utf-8'),
                    start=int(self.sample_offset + read['read_start']),
                    end=int(self.sample_offset + read['read_start'] + read['read_length'])
                )

    @rpc.method
    def unblock(self, read_id, read_end):
        """Fast forward to next read, returns new current sample. If the
        given read has already ended, no action is taken.

        :param read_id: read_id of read to unblock.
        :param read_end: expiration time for request, used only for logging.
        """
        cur_read = self.reads[self.current_read]
        if cur_read['read_id'].decode('utf-8') == read_id:
            next_read = self.reads[self.current_read + 1]
            jump = int(next_read['read_start'] - self.current_sample)
            self.sample_offset += jump
            self.logger.debug("Unblocking read {}. Skipping ahead {} samples to sample {}.".format(
                read_id, jump, self.current_sample
            ))
            return self.current_sample, True
        else:
            self.logger.info("Unblock received too late for read {}. Received at {}, read ended at {}".format(
                read_id, self.current_sample, read_end
            ))
            return self.current_sample, False


class ReplayFast5(rpc.AttrHandler):
    def __init__(self, fast5, channels, *args, good_class='strand', time_warp=1, **kwargs):
        """Replay multiple channels from a .fast5 file.

        :param fast5: input filename.
        :param channel: list of channels to simulate.
        :param good_class: read classification name of desirable reads.
        :param time_warp: time multiplier for playback speed.

        ..note:: `args` and `kwargs` are passed to `aiozmq.rpc.AttrHandler`.

        """
        super().__init__(*args, **kwargs)
        self.fast5 = fast5
        self.channels = channels
        self.good_class = good_class
        self.time_warp = time_warp
        self.logger = logging.getLogger('Replay Fast5')
        if self.time_warp < 1:
            self.logger.warn("time_warp set to slow down time, behaviour may be undefined.")

        self.replay_channels = {
            channel:ReplayChannel(self.fast5, channel, good_class=self.good_class, time_warp=self.time_warp)
            for channel in self.channels
        }

        for chan in self.replay_channels.values():
            chan.reset_time()

    @rpc.method
    def time_saved(self):
        """Difference between current real time and simulated time. (How far
        we've fast-forwarded the channel.)

        """
        return sum(chan.time_saved for chan in self.replay_channels.values())

    @rpc.method
    def total_good_reads(self):
        """Total number of good reads seen up to current simulation time."""
        return sum(chan.total_good_reads for chan in self.replay_channels.values())

    @rpc.method
    def cumulative_good_read_time(self):
        """Total duration of reads with classification equal to `good_class` until current time."""
        return sum(chan.cumulative_good_read_time_to_now for chan in self.replay_channels.values())

    @rpc.method
    def get_events(self, channel, n_events=400):
        """Return events from the start of the current read.

        :param nevents: maximum number of events to return.
        :param channel: channel for which to get events.

        :returns: Serialized event data, see :class:`Fast5Data`.
        """
        return self.replay_channels[channel].get_events(n_events=n_events)

    @rpc.method
    def get_raw(self, channel, seconds=0.5):
        """Return raw data from the start of the current read.

        :param channel: channel for which to get raw data.

        :returns: Serialized raw data, see :class:`Fast5Data`.
        """
        return self.replay_channels[channel].get_raw(seconds=seconds)

    @rpc.method
    def unblock(self, channel, read_id, read_end):
        """Fast forward to next read, returns new current sample. If the
        given read has already ended, no action is taken.

        :param channel: channel to unblock.
        :param read_id: read_id of read to unblock.
        :param read_end: expiration time for request, used only for logging.
        """
        return self.replay_channels[channel].unblock(read_id, read_end)


@asyncio.coroutine
def replay_server(fast5, channels, port, good_class, time_warp=1):
    """Create an RPC server to replay data from a .fast5 file.

    :param fast5: input .fast5 file for simulation.
    :param channels: list of channels to simulate.
    :param port: port on which to listen for clients.
    :param good_class: read classification name of desirable reads.
    
    :returns: instance of :class:`ReplayFast5`.
    """
    server = yield from rpc.serve_rpc(
        ReplayFast5(fast5, channels, good_class=good_class, time_warp=time_warp),
        bind='tcp://127.0.0.1:{}'.format(port),
        translation_table=translation_table,
    )
    return server


@asyncio.coroutine
def replay_client(port):
    """Create an RPC client to request raw/event data and send unblock requests.

    :param port: server port.
    """
    client = yield from rpc.connect_rpc(
        connect='tcp://127.0.0.1:{}'.format(port),
        translation_table=translation_table
    )
    return client


def read_until_demo(fast5, channels, port=5555):
    """Simple demo of read until server and client application.

    :param fast5: input .fast5 file.
    :param channels: list of channels to simulate.
    :param port: port on which to run server and cliant.
    """
    logger = logging.getLogger('ReadUntil App')
    good_class = 'strand'
    time_warp=2
    event_loop = asyncio.get_event_loop()
    set_wakeup()
    # Setup replay service
    event_loop.create_task(replay_server(
        fast5, channels, port, good_class, time_warp=time_warp
    ))

    @asyncio.coroutine
    def read_until(port):
        client = yield replay_client(port)
        counter = 0
        while True:
            #This can be mostly rewritten for any real app
            for channel in channels:
                read_block = yield from client.call.get_events(channel)
                if read_block is None:
                    logger.debug("Channel not in '{}' classification".format(good_class))
                else:
                    logger.info("Analysing {} events".format(len(read_block)))
                    yield from asyncio.sleep(0.5/time_warp)
                    counter += 1
                    # We decided to unblock
                    if counter % 10:
                        # We took too a long time
                        yield from asyncio.sleep(5/time_warp)
                    client.call.unblock(channel, read_block.info, read_block.end)
            yield from asyncio.sleep(5/time_warp)
    event_loop.create_task(read_until(port))

    try:
        event_loop.run_forever()
    except KeyboardInterrupt:
        pass


class ExpandRanges(argparse.Action):
    """Translate a str like 1,2,3-5,40 to [1,2,3,4,5,40]"""
    def __call__(self, parser, namespace, values, option_string=None):
        import re
        elts = []
        for item in values.replace(' ', '').split(','):
            mo = re.search(r'(\d+)-(\d+)', item)
            if mo is not None:
                rng = [int(x) for x in mo.groups()]
                elts.extend(list(range(rng[0], rng[1] + 1)))
            else:
                elts.append(int(item))
        setattr(namespace, self.dest, elts)

if __name__ == '__main__':
    freeze_support()
    logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    parser = argparse.ArgumentParser('Simple zmq server and client demo.')
    parser.add_argument('fast5', help='Input fast5.')
    parser.add_argument('channels', action=ExpandRanges, help='Fast5 channel for source data.')
    args = parser.parse_args()
    read_until_demo(args.fast5, [str(x) for x in args.channels])
