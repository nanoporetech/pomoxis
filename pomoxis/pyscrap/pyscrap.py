import importlib
import imp
import os
import sys

from fast5_research.fast5 import Fast5
from nanonet.segment import segment
from nanonet.eventdetection.filters import minknow_event_detect


from cffi import FFI
ffi = FFI()

import logging
logger = logging.getLogger(__name__)

def get_shared_lib(name):
    """Cross-platform resolution of shared-object libraries, working
    around vagueries of setuptools
    """
    try:
        # after 'python setup.py install' we should be able to do this
        lib_file = importlib.import_module(name).__file__
    except Exception as e:
        try:
            # after 'python setup.py develop' this should work
            lib_file = imp.find_module(name)[1]
        except Exception as e:
            raise ImportError('Cannot locate C library "{}".'.format(name))
        else:
            lib_file = os.path.abspath(lib_file)
    finally:
        library = ffi.dlopen(lib_file)
    return library


library = get_shared_lib('pyscrap')

ffi.cdef("""
    typedef struct {
        double start;
        float length;
        float mean, stdv;
        int pos, state;
    } event_t;

    typedef struct {
        unsigned int n, start, end;
        event_t * event;
    } event_table;

    typedef struct {
        float score;
        char * bases;
    } basecall_result;

    void free_char(char* ptr);
    basecall_result basecall_events(event_t * events, const int n_events);
""")


def build_events(events):
    """Transform numpy events into scrappie event_t[].

    .. note::

        the memory for the events is managed by python and therefore
        should not be freed on the C side otherwise a doublefree will occur.

    """
    n_events = len(events)
    event = ffi.new('event_t[]', n_events)
    for i in range(n_events):
        event[i].start = events[i]['start']
        event[i].length = events[i]['length']
        event[i].mean = events[i]['mean']
        event[i].stdv = events[i]['stdv']
        event[i].pos = -1
        event[i].state = -1
    return event


def basecall_events(events):
    """Basecall numpy events using scrappie."""
    c_events = build_events(events)
    result = library.basecall_events(c_events, len(events))
    if result.bases == ffi.NULL:
        return None
    basecall = ffi.string(result.bases)
    library.free_char(result.bases)
    return -result.score/len(events), basecall.decode()


def basecall_file(fname=None, event_detect=True):
    """Read event data from file and print scrappie basecall.

    :param fname: filename to read data from (if not given assumed
        to be given on command line.
    :param event_detect: do event detection?

    """
    is_main = False
    if fname is None: #called as entrypoint
        fname = sys.argv[1]
        is_main = True

    # magic numbers
    ed_params = {
        'window_lengths':[4, 8],
        'thresholds': [1.5, 9.0],
        'peak_height': 0.2,
    }

    with Fast5(fname) as fh:
        if event_detect:
            events = minknow_event_detect(
                fh.get_read(raw=True), fh.sample_rate, **ed_params
            )
        else:
            events = fh.get_read()
    events, _ = segment(events, section='template') 

    results = basecall_events(events)
    if results is None:
        return None
    if is_main:
        print("{} score={}\n{}".format(fname, *results))
    else:
        return results
