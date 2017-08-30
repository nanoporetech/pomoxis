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

extra_requires = {
    'nanonet':'nanonet'
}

exes = ['minimap', 'minimap2', 'miniasm', 'racon', 'bwa', 'samtools']
scrappie_path = os.path.join('submodules', 'scrappie', 'src')

extensions = []

extensions.append(Extension(
    'pyscrap',
    sources=[os.path.join(__pkg_path__, 'pyscrap', 'pyscrap.c')] + [
        os.path.join(scrappie_path, x) for x in ('decode.c', 'layers.c', 'networks.c', 'nnfeatures.c', 'util.c', 'scrappie_matrix.c')
    ],
    include_dirs=[scrappie_path],
    extra_compile_args=['-pedantic', '-Wall', '-std=c99', '-march=native', '-ffast-math', '-DUSE_SSE2', '-DNDEBUG'],
    libraries=['blas']
))

model_cmd_name = 'model'
class CustomBasecallModel(Command):
    description = 'Use a custom scrappie basecall model'
    user_options = [
        ('scrappie-model=', None, 'path to scrappie model header file'),
    ]

    def initialize_options(self):
        self.scrappie_model = None

    def finalize_options(self):
        if self.scrappie_model is None:
            raise ValueError('Command "{}" requires scrappie-model=<value> option.'.format(model_cmd_name))
        else:
            if not os.path.exists(self.scrappie_model):
                raise ValueError('Model header {} does not exist.'.format(self.scrappie_model))

    def run(self):
        dst = os.path.join(__path__, 'submodules/scrappie/src/nanonet_events.h')
        if not os.path.exists(dst):
            raise RuntimeError("The model destination does not exists. Has scrappie renamed it's files (again)?")
        self.announce('Copying {} to {}.'.format(self.scrappie_model, dst))
        shutil.copyfile(self.scrappie_model, dst)


setup(
    name=__pkg_name__,
    version=__version__,
    cmdclass={
        model_cmd_name:CustomBasecallModel
    },
    url='https://git.oxfordnanolabs.local/research/{}'.format(__pkg_name__),
    author='cwright',
    author_email='cwright@nanoporetech.com',
    description='Real-time assembly.',
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
            'split_fastx = {}.common.util:split_fastx_cmdline'.format(__pkg_name__),
            'fast_convert = {}.common.util:fast_convert'.format(__pkg_name__),
            'pyscrap = {}.basecall.pyscrap.pyscrap:basecall_file'.format(__pkg_name__),
            'epi3me = {}.apps.epi3me:main'.format(__pkg_name__),
            'read_until_filter = {}.apps.read_until_filter:main'.format(__pkg_name__),
            'pomoxis_path = {}:show_prog_path'.format(__pkg_name__),
            'bwa_rpc = {}.align.bwa:main'.format(__pkg_name__),
        ]
    },
    scripts=[
        'scripts/bwa_align',
        'scripts/mini_align',
        'scripts/mini_assemble',
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

