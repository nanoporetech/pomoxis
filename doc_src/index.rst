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
Makefile, see `.gitlab-ci.yml` for additional requirements if things fail.

To setup the environment run:

.. code-block:: bash

    git clone --recursive https://git/research/pomoxis.git
    cd pomoxis
    make install
    . ./venv/bin/activate

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

