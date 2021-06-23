__version__ = '0.3.8'

import argparse
import os
import sys
import platform
import subprocess

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
    parser = argparse.ArgumentParser(
        prog='pomoxis_path',
        description='Print the path of bundled executables.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('program', help='name of program.')
    args = parser.parse_args()
    print(get_prog_path(args.program))


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

