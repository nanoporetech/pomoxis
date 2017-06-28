import argparse
import asyncio
from collections import defaultdict, Counter
from multiprocessing import freeze_support
from timeit import default_timer as now

from aiozmq import rpc
import numpy as np

from nanonet.eventdetection.filters import minknow_event_detect

from pomoxis import set_wakeup
from pomoxis.provider import replayfast5
from pomoxis.align import bwa
from pomoxis.pyscrap import pyscrap

import logging
logger = logging.getLogger(__name__)


def read_until_align_filter(fast5, bwa_index, channels, start_port=5555, targets=['Ecoli', 'yeast'], whitelist=False):
    """Demonstration read until application using scrappie and bwa to filter
    reads by identity.

    :param fast5: input bulk .fast5 file for read until simulation.
    :param bwa_index: list of bwa index files (just basename without additional
        suffixes).
    :param channels: list of channels to simulate.
    :param start_port: port on which to run .fast5 replay server, bwa alignment
        server will be run on the following numbered port.
    :param targets: list of reference names. If `whitelist` is `False`, reads
        aligning to these references will be ejected. If `whitelist` is `True`
        any read alinging to a reference other than those contained in this
        list will be ejected. Unidentified reads are never ejected.
    :param whitelist: see `target`.
    """

    logger = logging.getLogger('ReadUntil App')
    good_class = 'strand'
    time_warp=2
    event_loop = asyncio.get_event_loop()
    set_wakeup()

    port = start_port
    # Setup replay service
    replay_port = port
    replay_server = event_loop.create_task(replayfast5.replay_server(
        fast5, channels, replay_port, good_class, time_warp=time_warp
    ))
    port += 1
    # Setup alignment service
    align_port = port
    align_server = event_loop.create_task(bwa.align_server(
        bwa_index, align_port
    ))


    identified_reads = {}
    unidentified_reads = defaultdict(list)
    unblocks = Counter()

    ###
    # The read until app
    @asyncio.coroutine
    def poll_data(port):
        align_client = yield from bwa.align_client(align_port) 
        replay_client = yield from replayfast5.replay_client(replay_port)
        yield from asyncio.sleep(5)
        start_time = now()
        target_count = 0
        while True:
            time_saved = yield from replay_client.call.time_saved()
            total_pore_time = (now() - start_time) * len(channels)
            total_strand_time = yield from replay_client.call.cumulative_good_read_time()
            try:
                pore_time_saved = time_saved/total_pore_time
            except:
                pore_time_saved = 0
            try:
                strand_time_saved = time_saved/total_strand_time
            except:
                strand_time_saved = 0
            logger.info("Total pore time saved: {:.2f}% [{:.2f}/{:.2f}]".format(
                100 * pore_time_saved, time_saved, total_pore_time
            ))
            logger.info("Total strand time saved: {:.2f}% [{:.2f}/{:.2f}]".format(
                100 * strand_time_saved, time_saved, total_strand_time
            ))
            reads_analysed = set(identified_reads.keys()) | set(unidentified_reads.keys())
            all_reads = yield from replay_client.call.total_good_reads()
            ided = len(identified_reads)
            unided = len(unidentified_reads)
            missed = all_reads - len(reads_analysed)
            logger.info("identified/unidentified/missed reads: {}/{}/{}.".format(ided, unided, missed))
            logger.info("Unblocks (timely/late): {}/{}.".format(unblocks[True], unblocks[False]))
            logger.info("Total good reads: {}".format(target_count))

            for channel in channels:
                read_block = yield from replay_client.call.get_raw(channel)
                if read_block is None:
                    logger.debug("Channel not in '{}' classification".format(good_class))
                elif read_block.info in identified_reads:
                    logger.debug("Skipping because I've seen before.")
                    continue
                else:
                    logger.debug("Analysing {} samples".format(len(read_block)))
                    sample_rate = read_block.sample_rate
                    events = minknow_event_detect(
                        read_block, read_block.sample_rate, **{
                            'window_lengths':[3, 6], 'thresholds':[1.4, 1.1],
                            'peak_height':0.2
                        }
                    )
                    if len(events) < 100:
                        continue

                    #TODO: do this in a process pool
                    score, basecall = pyscrap.basecall_events(events)
                    #TODO: check sanity of basecall
                    if len(basecall) < 100:
                        continue

                    alignment, returncode = yield from align_client.call.align(basecall)
                    hits = []
                    if returncode != 0:
                        logger.warning('Alignment failed for {}'.format(read_block.info))
                    else:
                        recs = [x for x in alignment.split('\n') if len(x) > 0 and x[0] != '@']
                        for r in recs:
                            fields = r.split('\t')
                            if fields[2] != '*':
                                hits.append(fields[2])
                    logger.debug('{} aligns to {}'.format(read_block.info, hits))

                    if len(hits) == 1:
                        identified_reads[read_block.info] = hits[0]
                        # maybe got 0 or >1 previously
                        #TODO: there are some edges cases here
                        try:
                            del unidentified_reads[read_block.info]
                        except KeyError:
                            pass
                    else:
                        unidentified_reads[read_block.info].extend(hits)

                    if read_block.info in identified_reads:
                        good_read = whitelist
                        if identified_reads[read_block.info] not in targets:
                            good_read = not whitelist

                        if not good_read:
                            logger.info('Attempting to unblock channel {} due to contaminant.'.format(channel))
                            _, good_unblock = yield from replay_client.call.unblock(channel, read_block.info, read_block.end)
                            unblocks[good_unblock] += 1
                        else:
                            target_count += 1

    event_loop.create_task(poll_data(port))

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

def main():
    freeze_support()
    logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    parser = argparse.ArgumentParser(description="""
Read until with alignment filter. The read until simulator takes as input a
"bulk .fast5" file. These are an optional output of MinKnow which contains
the signal data for all channels of a Oxford Nanopore Technologies’ device
(including periods of time where a channel is not undergoing an strand
translocation). Outputting a bulk .fast5 can be configured when the user
starts and experiment in MinKnow.
""")
    parser.add_argument('fast5', help='Input fast5.')
    parser.add_argument('channels', action=ExpandRanges, help='Fast5 channel for source data.')
    parser.add_argument('bwa_index', nargs='+', help='Filename path prefix for BWA index files.')
    args = parser.parse_args()
    read_until_align_filter(args.fast5, args.bwa_index, [str(x) for x in args.channels])


if __name__ == '__main__':
    main()
