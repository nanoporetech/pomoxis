Prototype real-time analysis components
=======================================

[![Build Status](https://travis-ci.org/nanoporetech/pomoxis.svg?branch=master)](https://travis-ci.org/nanoporetech/pomoxis)

Pomoxis contains a set of services to perform analysis of squiggles as they are
produced in real-time along with fast pipelines for generating draft assemblies.

Documentation can be found at https://nanoporetech.github.io/pomoxis/.

Currently pomoxis supports only unix-like environments.
  

Build
-----

Pomoxis should be installed inside a virtual environment. A Makefile is
provided to create a fresh environment, and to fetch, compile and install
all direct dependencies into the environment.

To setup the environment run:

    git clone --recursive https://github.com/nanoporetech/pomoxis
    cd pomoxis
    make install
    . ./venv/bin/activate

Running the above within a pre-exisiting virtual environnment may well fail;
advanced may wish to simply run the `setup.py` file in the standard manner
after compiling the third party programs as in the `Makefile`.


Extras
------

The distribution bundles some common bioinformatics tools (some of which are not
currently used by pomoxis itself):

* miniasm
* minimap
* racon
* bwa
* samtools

These will be compiled and installed into the virtual environment created as above.
