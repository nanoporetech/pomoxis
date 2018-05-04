import os
import sys
import shutil
import re
import shutil
import platform
from glob import glob
from setuptools import setup, find_packages, Extension
from setuptools import Distribution, Command
from setuptools.command.install import install
import pkg_resources

__path__ = os.path.dirname(__file__)
__pkg_name__ = 'pomoxis'
__pkg_path__ = os.path.join(os.path.join(__path__, __pkg_name__))

# Get the version number from __init__.py, and exe_path
verstrline = open(os.path.join(__pkg_name__, '__init__.py'), 'r').read()
vsre = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(vsre, verstrline, re.M)
if mo:
    __version__ = mo.group(1)
else:
    raise RuntimeError('Unable to find version string in "{}/__init__.py".'.format(__pkg_name__))

dir_path = os.path.dirname(__file__)
install_requires = []
with open(os.path.join(dir_path, 'requirements.txt')) as fh:
    reqs = (
        r.split('#')[0].strip()
        for r in fh.read().splitlines() if not r.startswith('#')
    )
    for req in reqs:
        if req == '':
            continue
        if req.startswith('git+https'):
            req = req.split('/')[-1].split('@')[0]
        install_requires.append(req)


# Read these from Makefile in the event someone modified list there
exes = []
with open(os.path.join(dir_path, 'Makefile')) as fh:
    for line in fh.readlines():
        tokens = line.split('=')
        if tokens[0] == 'BINARIES':
            exes = tokens[1].split()
            break

extensions = []
extra_requires={}

setup(
    name=__pkg_name__,
    version=__version__,
    url='https://github.com/nanoporetech/{}'.format(__pkg_name__),
    author='cwright',
    author_email='cwright@nanoporetech.com',
    description='Real-time analysis components.',
    dependency_links=[],
    ext_modules=extensions,
    install_requires=install_requires,
    tests_require=[].extend(install_requires),
    extras_require=extra_requires,
    packages=find_packages(exclude=['*.test', '*.test.*', 'test.*', 'test']),
    package_data={},
    zip_safe=False,
    test_suite='discover_tests',
    #place binaries as package data, below we'll copy them to standard path in dist
    data_files=[
        ('exes', [
            'bincache/{}'.format(x, x) for x in exes
        ])
    ],
    entry_points={
        'console_scripts': [
            'align_serve = {}.align.common:main'.format(__pkg_name__),
            'catalogue_errors = {}.common.catalogue_errors:main'.format(__pkg_name__),
            'coverage_from_bam = {}.common.coverage_from_bam:main'.format(__pkg_name__),
            'epi3me = {}.apps.epi3me:main'.format(__pkg_name__),
            'fast_convert = {}.common.util:fast_convert'.format(__pkg_name__),
            'long_fastx = {}.common.util:extract_long_reads'.format(__pkg_name__),
            'pomoxis_path = {}:show_prog_path'.format(__pkg_name__),
            'read_until_filter = {}.apps.read_until_filter:main'.format(__pkg_name__),
            'stats_from_bam = {}.common.stats_from_bam:main'.format(__pkg_name__),
            'summary_from_stats = {}.common.summary_from_stats:main'.format(__pkg_name__),
            'simulate_calls = {}.common.simulate_calls:main'.format(__pkg_name__),
            'split_fastx = {}.common.util:split_fastx_cmdline'.format(__pkg_name__),
            'subsample_bam = {}.common.subsample_bam:main'.format(__pkg_name__),
            'trim_alignments = {}.common.trim_alignments:main'.format(__pkg_name__),
            'common_errors_from_bam = {}.common.common_errors_from_bam:main'.format(__pkg_name__),
        ]
    },
    scripts=[
        'scripts/bwa_align',
        'scripts/mini_align',
        'scripts/mini_assemble',
        'scripts/assess_assembly',
        'scripts/intersect_assembly_errors',
    ]
)


# Nasty hack to get binaries into bin path
print("\nCopying utility binaries to your path.")
class GetPaths(install):
    def run(self):
        self.distribution.install_scripts = self.install_scripts
        self.distribution.install_libbase = self.install_libbase

def get_setuptools_script_dir():
    # Run the above class just to get paths
    dist = Distribution({'cmdclass': {'install': GetPaths}})
    dist.dry_run = True
    dist.parse_config_files()
    command = dist.get_command_obj('install')
    command.ensure_finalized()
    command.run()

    src_dir = glob(os.path.join(dist.install_libbase, 'pomoxis-*', 'exes'))[0]
    for exe in (os.path.join(src_dir, x) for x in os.listdir(src_dir)):
        print("Copying", os.path.basename(exe), '->', dist.install_scripts)
        shutil.copy(exe, dist.install_scripts)
    return dist.install_libbase, dist.install_scripts

get_setuptools_script_dir()

