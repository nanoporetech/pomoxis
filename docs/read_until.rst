Read Until Simulator
--------------------

The Fast5 playback in pomoxis can act as a read until simulator. There is a
minimal API to fetch raw or event data from current reads and to unblock
channels.

Internally the playback classes construct an index across reads in a channel.
When any request is made, the current simulation sample is calculated from
the current (real) time and an accumulated offset. The simulation sample is used
to calculate the current read and associated raw and event data. When unblock
requests are made the simulation records the difference between the current
simulated sample and the start sample of the following read. This sample offset
is accumulated as subsequent unblock requests are made.

.. note::

    The read until simulator takes as input a "bulk .fast5" file. These are an
    optional output of MinKnow which contains the signal data for all channels
    of a Oxford Nanopore Technologies' device (including periods of time where
    a channel is not undergoing an strand translocation). Outputting a bulk
    .fast5 can be configured when the user starts and experiment in MinKnow.

To create an RPC server for an input .fast5 file:

.. code-block:: python

    from pomoxis.provider import replayfast5

    loop = asyncio.get_event_loop()
    loop.create_task(replayfast5.replay_server(
        fast5, channels, port, good_class, time_warp=time_warp
    )

The first three arguments should be self explanatory. `good_class` will set the
read classification of reads which are determined to be valid, usually this
should be set to `'strand'`. If a request to the server for data is made when
the current read is not of this class, the server will return `None`. The
optional `time_warp` parameter specifies a multiplier for the progression of
time, the default is `1`, but may be set to values larger than this (smaller
and things might break).

To create a client of the above server:

.. code-block:: python

    from pomoxis.provider import replayfast5

    client = yield replayfast5.replay_client(port)

On its own this is pretty useless, but one can quickly build a read until app

.. code-block:: python

    @asyncio.coroutine
    def read_until(port)
        client = yield replayfast5.replay_client(port)
        while True:
        for channel in channels:
            events = yield from client.call.get_events(channel)
            if read_block is None:
                print("Channel not in '{}' classification".format(good_class))
            else:
                client.call.unblock(channel, events.info, events.end)
    event_loop.create_task(read_until(port))

Here, the unblock call passes backsome meta data on the events. This is merely
for programming convenience, the API may be simplified in the future.

A more fullsome example of an concurrent alignment system is given in
:py:mod:`pomoxis.apps.read_until_filter`.
