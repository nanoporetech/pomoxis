__version__ = '0.1.0'

import os
import sys
import asyncio
import platform
import subprocess

@asyncio.coroutine
def wakeup():
    """No-op coroutine."""
    while True:
        yield from asyncio.sleep(1)


def set_wakeup():
    """Workaround suppression of `KeyboardInterrrupt` on Windows."""
    if platform.system() == 'Windows':
        #SO/27480967
        asyncio.async(wakeup())


def get_prog_path(prog):
    """Get the absolute path of bundled executables.

    :param prog: programme (file) name.

    :returns: absolute path to executable.
    """
    prog_path = os.path.join(
        os.path.dirname(__file__), '..', 'exes', prog
    )
    if not (os.path.isfile(prog_path) and os.access(prog_path, os.X_OK)):
        raise RuntimeError('Cannot find executable "{}".'.format(prog))
    return prog_path


def show_prog_path():
    """Print the path of bundled executables."""
    print(get_prog_path(sys.argv[1]))


def run_prog(prog, args, stdout=None):
    """Run one of the bundled executables.

    :param prog: programme name.
    :param args: commandline arguments for programme.
    :param stdout: filehandle/path for stdout of programme.

    :returns: result of `subprocess.call`.
    """
    prog_path = get_prog_path(prog)
    if stdout is None:
        return subprocess.check_output([prog_path] + args)
    else:
        with open(stdout, 'w') as fh:
            return subprocess.call([prog_path] + args, stdout=fh, stderr=open(os.devnull, 'w'))

