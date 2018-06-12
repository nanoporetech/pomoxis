Prototype real-time analysis components
=======================================

[![Build Status](https://travis-ci.org/nanoporetech/pomoxis.svg?branch=master)](https://travis-ci.org/nanoporetech/pomoxis)

Pomoxis started life as a set of services to perform analysis of squiggles
as they are produced in real-time. This functionality still exists
although has been complemented with a number of offline analyses for
generating and analysing draft assemblies. Many of these tools are used by
the research data analysis group at Oxford Nanopore Technologies.

Documentation can be found at https://nanoporetech.github.io/pomoxis/.


Compatibility
-------------

Pomoxis is developed on Ubuntu 16.04, other recent linuxes should be
equally compatible (see build notes below. Pomoxis is known to work on
at least some MacOS High Sierra configurations, though some components
(notably scrappy) are known to not work on some MacOS configurations
(combinations of OS and xcode versions).


Build
-----

Pomoxis should be installed inside a virtual environment. A Makefile is
provided to create a fresh environment, and to fetch, compile and install
all direct dependencies into the environment.

> Before installing pomoxis is may be required to install some prerequisite
> packages, best installed by a package manager. On Ubuntu these are:
> * gcc-4.9
> * g++-4.9
> * zlib1g-dev
> * libncurses5-dev
> * python3-all-dev
> * libhdf5-dev
> * libatlas-base-dev
> * libopenblas-base
> * libopenblas-dev
> * libbz2-dev
> * liblzma-dev
> * libffi-dev

To setup the environment run:

    git clone --recursive https://github.com/nanoporetech/pomoxis
    cd pomoxis
    make install
    . ./venv/bin/activate
    

The installation of porechop (https://github.com/rrwick/Porechop)
requires a newer compiler than is a available on some systems. It may therefore
be necessary to install a newer compiler and set environment variables before
the `make install` step:

    # For porechop to be compiled on older systems set these, e.g.:
    export CXX="g++-4.9" CC="gcc-4.9"

Running the above within a pre-exisiting virtual environnment may well fail;
advanced may wish to simply run the `setup.py` file in the standard manner
after compiling the third party programs as in the `Makefile`.


Extras
------

The distribution bundles some common bioinformatics tools (some of which are not
currently used by pomoxis itself, but are used in the offline analysis scripts):

* miniasm
* minimap2
* racon
* bwa
* samtools
* porechop

These will be compiled and installed into the virtual environment created as above.
