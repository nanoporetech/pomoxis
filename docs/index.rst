Welcome to Pomoxis' documentation!
==================================

Pomoxis comprises APIs and command line tools for interacting and analysing
Oxford Nanopore Technologies' data in real time.

Most tools expose an ZeroMQ RPC service for communication such that components
can run on disperate computing resources, though one could build apps with
components all running locally under a single event loop.

Installation
------------

Pomoxis should be installed inside a virtual environment. A Makefile is
provided to fetch, compile and install all direct dependencies into an
environment. (Some additional build dependencies are not installed via this
Makefile, see `.travis.yml` for additional requirements if things fail.

To setup the environment run:

.. code-block:: bash

    git clone --recursive https://github.com/nanoporetech/pomoxis.git
    cd pomoxis
    # For porechop to be compiled on older systems set these, e.g.:
    #    export CXX="g++-4.9" CC="gcc-4.9"
    make install
    . ./venv/bin/activate

The installation of porechop `porechop <https://github.com/rrwick/Porechop>`_
requires a newer compiler than is a available on some systems. It may therefore
be necessary to install a newer compiler and set variables as in the above.

See :doc:`examples` for common bundled workflows.

Contents
--------

.. toctree::
   :maxdepth: 2

   examples
   read_until

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

