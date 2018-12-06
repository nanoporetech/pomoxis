
Pomoxis - bioinformatics tools for nanopore research 
====================================================

Pomoxis comprises a set of basic bioinformatic tools tailored to nanopore
sequencing. Notably tools are included for generating and analysing draft
assemblies. Many of these tools are used by the research data analysis
group at Oxford Nanopore Technologies.

See :doc:`examples` for common simple tasks.

Features
--------

 * Wraps third party tools with known good default parameters
   and methods of use.
 * Creates an isolated environment with all third-party tools.
 * Streamlines common short analysis chains.
 * Includes a nanopore read simulator.
 * Server/client components for minimap2 and bwa.
 * Integrates into `katuali <https://github.com/nanoporetech/katuali>`_
   for performing more complex analysis pipelines.
 * Open source (Mozilla Public License 2.0).


Compatibility
-------------

Pomoxis is developed on Ubuntu 16.04, other recent linuxes should be
equally compatible (see build notes below). Pomoxis is known to work on
at least some MacOS High Sierra configurations, though some components,
notably scrappy, are known to not work on some MacOS configurations
(combinations of OS and xcode versions).


Installation
------------

Pomoxis will install itself into a an isolated virtual environment. The
installation will fetch, compile, and install all direct dependencies into the
environment.

.. note::

    Before installing pomoxis is may be required to install some prerequisite
    packages, best installed by a package manager. On Ubuntu these are:

    gcc-4.9 g++-4.9 zlib1g-dev libncurses5-dev python3-all-dev libhdf5-dev
    libatlas-base-dev libopenblas-base libopenblas-dev libbz2-dev liblzma-dev
    libffi-dev make python-virtualenv cmake wget bzip2

To setup the environment run:

.. code-block:: bash

    # export CXX="g++-4.9" CC="gcc-4.9"  # see below
    git clone --recursive https://github.com/nanoporetech/pomoxis
    cd pomoxis
    make install
    . ./venv/bin/activate
    

The installation of porechop (https://github.com/rrwick/Porechop)
requires a newer compiler than is a available on some systems. It may therefore
be necessary to install a newer compiler and set environment variables before
the `make install` step:

.. code-block:: bash

    # For porechop to be compiled on older systems set these, e.g.:
    export CXX="g++-4.9" CC="gcc-4.9"

Note also that racon requires at least `gcc>=4.8.5` to
`compile smoothly <https://github.com/isovic/racon/issues/57>`_.


**Alternative Installation Methods**

Running the above within a pre-exisiting virtual environnment may well fail;
advanced users may wish to simply run the `setup.py` file in the standard
manner after compiling the third party programs as in the `Makefile`.
`thomcuddihy <https://github.com/thomcuddihy>`_ has
`sketched <https://github.com/nanoporetech/pomoxis/issues/19#issuecomment-433255390>`_
how to use pomoxis within a conda enviroment.


Contents
--------

.. toctree::
   :maxdepth: 2

   examples

Full API reference
------------------

.. toctree::
   :maxdepth: 3
      
   pomoxis


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

