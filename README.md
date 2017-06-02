Prototype fast assembly
=======================

Pomoxis contains a set of services to perform analysis of squiggles as they are
produced in real-time along with fast pipelines for generating draft assemblies.


Build
-----

Pomoxis should be installed inside a virtual environment. A Makefile is
provided to fetch, compile and install all direct dependencies into an
environment.

To setup the environment run:

    git clone --recursive https://git/research/pomoxis.git
    cd pomoxis
    make install
    . ./venv/bin/activate


Extras
------

The distribution bundles some common bioinformatics tools (some of which are not
currently used by pomoxis itself):

* miniasm
* minimap
* racon
* bwa
* samtools

The "raw" binaries are distributed under e.g.

    venv/lib/python{pyversion}/site-packages/pomoxis-{version}-py{pyversion}.egg/exes/
    
and wrapping python entry points are installed into the virtual environment.


Offline assembly
----------------

An offline assembly benchmarking script is available within the virtual
environment as:

    mini_assemble -i input.fastq -o output_directory -p my_assembly
