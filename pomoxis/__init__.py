__version__ = '0.0.1'

import os
import sys
import asyncio
import platform
import subprocess

@asyncio.coroutine
def wakeup():
    while True:
        yield from asyncio.sleep(1)


def set_wakeup():
    if platform.system() == 'Windows':
        #SO/27480967
        asyncio.async(wakeup())


def get_prog_path(prog):
    prog_path = os.path.join(
        os.path.dirname(__file__), '..', 'exes', prog
    )
    if not (os.path.isfile(prog_path) and os.access(prog_path, os.X_OK)):
        raise RuntimeError('Cannot find executable "{}".'.format(prog))
    return prog_path


def show_prog_path():
    print(get_prog_path(sys.argv[1]))


def run_prog(prog, args, stdout=None):
    prog_path = get_prog_path(prog)
    if stdout is None:
        return subprocess.check_output([prog_path] + args)
    else:
        with open(stdout, 'w') as fh:
            return subprocess.call([prog_path] + args, stdout=fh, stderr=open(os.devnull, 'w'))

