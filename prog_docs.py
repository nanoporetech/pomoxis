#!/usr/bin/env python

import io
import os
import subprocess
import sys

header = """
Pomoxis Programs
================

Below you will find a listing with brief description of all the tools
within pomoxis. For more information concerning each tool simply run
it on the command line with the ``--help`` option.

More complex workflows (using tools from pomoxis) can be found within the
`katuali <https://github.com/nanoporetech/katuali>`_ software package.


"""

print(header)

def create_py_docs(prog, err=False):
    output = list()
    p = subprocess.run([prog, '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        p.check_returncode()
    except subprocess.CalledProcessError as e:
        sys.stderr.write("- Failed to run {} -h.\n".format(prog))
        sys.stderr.write("- STDOUT:\n {}".format(p.stdout.decode()))
        sys.stderr.write("- STDERR:\n {}".format(p.stderr.decode()))
        raise e

    found_spacer = False
    buffer = list()
    text = p.stderr if err else p.stdout
    for line in io.StringIO(text.decode()).readlines():
        line = line.rstrip()
        if found_spacer:
            buffer.append(line)
            break
        buffer.append(line)
        if line == '':
            found_spacer = True
    usage_lines = buffer[0:-2]
    output.append(prog)
    output.append('*'*len(prog))
    output.append('')
    output.append(buffer[-1])
    output.append('')
    output.append('.. code-block:: bash')
    output.append('')
    for b in usage_lines:
        b = b.replace('usage: ','')
        output.append('    {}'.format(b))
    output.append('')
    return output

program_blocks = dict()

# grab docs for bash scripts
sys.stderr.write("Gerenating docs for scripts.\n")
for r, d, f in os.walk('scripts'):
    for script in f:
        sys.stderr.write("* {}\n".format(script))
        program_blocks[script] = create_py_docs(script, err=True)

# create docs for all python entry points
sys.stderr.write("Gerenating docs for entry points.\n")
start = False
with open('setup.py', 'r') as fh:
    for line in fh.readlines():
        if start:
            if line.strip()[0] == ']':
                break
            else:
                prog_name = line.split()[0].replace("'",'')
                sys.stderr.write("* {}\n".format(prog_name))
                program_blocks[prog_name] = create_py_docs(prog_name)
        elif line.find("'console_scripts': [") != -1:
            start = True

for prog in sorted(program_blocks.keys()):
    for line in program_blocks[prog]:
        print(line)
