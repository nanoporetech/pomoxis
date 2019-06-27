import asyncio
import os
from watchdog.observers import Observer
from watchdog.events import RegexMatchingEventHandler
from watchdog.utils import has_attribute, unicode_paths
import logging

logger = logging.getLogger(__name__)

EVENT_TYPE_MOVED = 'moved'
EVENT_TYPE_DELETED = 'deleted'
EVENT_TYPE_CREATED = 'created'
EVENT_TYPE_MODIFIED = 'modified'


def wait_for_file(fname):
    """Block until a filesize remains constant."""
    # watchdog cannot tell us when a newly created file has finished being
    #    written. Let's just wait until the filesize stops changing
    size = None
    while True:
        try:
            newsize = os.path.getsize(fname)
        except:
            newsize = None
        else:
            if newsize is not None and size == newsize:
                break
        print(size, newsize, fname)
        size = newsize

class AIORegexMatchingEventHandler(RegexMatchingEventHandler):
    def __init__(self, callback, loop=None, **kwargs):
        """asyncio compatible minimal regex matching event
        handler for watchdog.

        :param callback: function to apply to filenames.
        :param loop: ayncio-like event loop.
        """

        RegexMatchingEventHandler.__init__(self, **kwargs)
        self._loop = loop or asyncio.get_event_loop()
        self.callback = callback

    @asyncio.coroutine
    def _process_file(self, event):
        """Process an event when a file is created (or moved).

        :param event: watchdog event.

        :returns: result of applying `callback` to watched file.
        """
        if event.event_type == EVENT_TYPE_CREATED:
            fname = event.src_path
        else:
            fname = event.dest_path
        # need to wait for file to be closed
        wait_for_file(fname)
        result = yield from self.callback(fname)
        return result

    def dispatch(self, event):
        """Dispatch an event after filtering. We handle
        creation and move events only.
        
        :param event: watchdog event.

        :returns: None
        """
        if event.event_type not in (EVENT_TYPE_CREATED, EVENT_TYPE_MOVED):
            return
        if self.ignore_directories and event.is_directory:
            return

        paths = []
        if has_attribute(event, 'dest_path'):
            paths.append(unicode_paths.decode(event.dest_path))
        if event.src_path:
            paths.append(unicode_paths.decode(event.src_path))

        if any(r.match(p) for r in self.ignore_regexes for p in paths):
            return

        if any(r.match(p) for r in self.regexes for p in paths):
            self._loop.call_soon_threadsafe(asyncio.ensure_future, self._process_file(event))


class Watcher(object):
    def __init__(self, path, event_handler, recursive=False):
        """Wrapper around common watchdog idiom.

        :param path: path to watch for new files.
        :param event_handler: subclass of watchdog.events.FileSystemEventHandler.
        :param recursive: watch path recursively?
        """
        self.observer = Observer()
        self.observer.schedule(event_handler, path, recursive)

    def start(self):
        """Start observing path."""
        self.observer.start()

    def stop(self):
        """Stop observing path."""
        self.observer.stop()
        self.observer.join()


@asyncio.coroutine
def watch_path(path, callback, recursive=False, regexes=['.*\.fast5$']):
    """Watch a filepath indefinitely for new files, applying callback to files.

    :param path: path to watch.
    :param callback: callback to apply to newly found files.
    :param recursive: recursively watch path?
    :param regexes: filter files by applying regex.
    """
    handler = AIORegexMatchingEventHandler(callback=callback, regexes=regexes)
    watch = Watcher(path, event_handler=handler, recursive=recursive)
    logging.info('Starting to watch {} for new files matching {}.'.format(path, regexes))
    watch.start()
