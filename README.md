<img src="images/ONT_3Line_Black_RGB_1000px.png" alt="Oxford Nanopore Technologies logo" height="128">

Pomoxis - bioinformatics tools for nanopore research 
====================================================

Pomoxis comprises a set of basic bioinformatic tools tailored to nanopore
sequencing. Notably tools are included for generating and analysing draft
assemblies. Many of these tools are used by the research data analysis
group at Oxford Nanopore Technologies.

Documentation can be found at https://nanoporetech.github.io/pomoxis/.

© 2018 Oxford Nanopore Technologies Ltd.


Features
--------
 
 * Wraps third party tools with known good default parameters
   and methods of use.
 * Creates an isolated environment with all third-party tools.
 * Can be installed with conda.
 * Streamlines common short analysis chains.
 * Open source (Mozilla Public License 2.0).


Compatibility
-------------

Pomoxis is developed on Ubuntu 24.04, other recent linux distributions 
should be equally compatible (see build notes below).

Installation
------------

Much of pomoxis's functionality is dependent on thirty party tools. These
can be provided by the user, or can be installed from source with the 
help of the provided `Makefile`

**Installation with conda**

Pomoxis is available on [bioconda](https://bioconda.github.io/recipes/pomoxis/)
and so can be most easily installed with:

    conda install pomoxis

**Installation with pip**
  
Pomoxis can be installed without binary dependencies using pip:

    python3 -m venv venv
    . venv/bin/activate
    pip install --update pip
    pip install pomoxis

Using this method requires the user to add binaries on the `PATH`:

 * [minimap2](https://github.com/lh3/minimap2)
 * [samtools](https://github.com/samtools/samtools)
 * [bcftools](https://github.com/samtools/bcftools/)
 * [seqkit](https://github.com/shenwei356/seqkit)


**Installation from source**

With this method pomoxis will install itself into a an isolated virtual
environment. The installation will fetch, compile, and install all direct
dependencies into the environment. Use this method if you do not wish to
use conda, but will not be providing the third-party binaries.

> Before installing pomoxis is may be required to install some prerequisite
> packages, best installed by a package manager. On Ubuntu these are:
> * gcc
> * zlib1g-dev
> * libncurses5-dev
> * libhdf5-dev
> * libbz2-dev
> * liblzma-dev
> * make
> * wget (for fetching modules from github)
> * bzip2 (for extracting those modules)

To setup the environment run:

    git clone --recursive https://github.com/nanoporetech/pomoxis
    cd pomoxis
    make install
    . ./venv/bin/activate
    

**Installation from source without compiling third-party binaries**

Running the above within a pre-exisiting (virtual) environnment may well fail;
advanced may wish to simply run

    pip install ./pomoxis

in the standard manner after compiling the third party programs listed below
and ensuring they are present on the `PATH`. The `Makefile` provides a rule to 
copy the binaries into the python interpreter path if they are placed within a
directory named `bincache`. To make use of this facility run:

    make venv
    make copy_binaries


Third party binaries
--------------------

The distribution bundles some common bioinformatics tools:

* minimap2
* samtools
* bcftools
* seqkit


Help
----

**Licence and Copyright**

© 2018 Oxford Nanopore Technologies Ltd.

`pomoxis` is distributed under the terms of the Mozilla Public License 2.0.

**Research Release**

Research releases are provided as technology demonstrators to provide early
access to features or stimulate Community development of tools. Support for
this software will be minimal and is only provided directly by the developers.
Feature requests, improvements, and discussions are welcome and can be
implemented by forking and pull requests. However much as we would
like to rectify every issue and piece of feedback users may have, the 
developers may have limited resource for support of this software. Research
releases may be unstable and subject to rapid iteration by Oxford Nanopore
Technologies.
