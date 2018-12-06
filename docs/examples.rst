Pomoxis Examples
================

Below you will find examples of some key tools within pomoxis. Not all programs
are listed, plenty more are available, and all have help information available
from their commandline interface.

More complex workflows (using tools from pomoxis) can be found within the
`katuali <https://github.com/nanoporetech/katuali>`_ software package.


Read Alignment with minimap2
----------------------------

A basic alignment workflow with commonl options is captured within the
`mini_align` program. This performs alignment of nanopore reads to a reference
sequence using appropriate parameters, performing the common sorting and indexing
post alignment steps, and exposes a minimal set of (optional) commonly useful
settings.

.. code-block:: bash

    mini_align [-h] -r <reference> -i <fastq>
    
    Align fastq/a formatted reads to a genome using minimap2.
    
        -h  show this help text.
        -r  reference, should be a fasta file. If correspondng minimap indices
            do not exist they will be created. (required).
        -i  fastq/a input reads (required).
        -I  split index every ~NUM input bases (default: 16G, this is larger
            than the usual minimap2 default).
        -a  aggressively extend gaps (sets -A1 -B2 -O2 -E1 for minimap2).
        -P  filter to only primary alignments (i.e. run samtools view -F 2308)
        -n  sort bam by read name.
        -c  chunk size. Input reads/contigs will be broken into chunks
            prior to alignment.
        -t  alignment threads (default: 1).
        -p  output file prefix (default: reads).
        -m  fill MD tag.


Fast de-novo assembly
---------------------

From a `.fastq` file containing basecalls one can perform an assembly and
consensus using the `mini_assemble` tool.

.. code-block:: bash

    mini_assemble [-h] -i <fastq>
    
    Assemble fastq/fasta formatted reads and perform POA consensus.
    
        -h  show this help text.
        -i  fastx input reads (required).
        -q  use qualities as is (default: false).
        -r  reference fasta for reference-guided consensus (instead of de novo assembly)
        -o  output folder (default: assm).
        -p  output file prefix (default: reads).
        -t  number of minimap and racon threads (default: 1).
        -m  number of racon rounds (default: 4).
        -n  number of racon shuffles (default: 1. If not 1, should be >10).
        -w  racon window length (default: 500).
        -c  trim adapters from reads prior to everything else.
        -e  error correct longest e% of reads prior to assembly, or an estimated coverage (e.g. 2x).
            For an estimated coverage, the -l option must be a fastx or a length (e.g. 4.8mb).
        -l  Reference length, either as a number (e.g. 4.8mb) or fastx from which length will be calculated.
        -x  log all commands before running.

For E.coli on a computer with 16 CPUs this should take around 5 minutes for
a 50-fold coverage dataset.

When a reference fasta is provided, a reference-guided consensus is generated
instead of a de novo assembly.

The options `-c` and `-e` can be used to improve the assembly quality at the 
expense of speed (particularly `-e`).


Consensus assessment
--------------------

To evaluate the accuracy of a consensus, the `assess_assembly` program will
perform the necessary alignments between an assembly and reference
sequence and provide post-processing on the aligner output to produce a variety
of summary information and graphs.

.. code-block:: bash

    assess_assembly [-h] -r <reference> -i <fastq>
    
    Calculate accuracy statistics for an assembly.
    
        -h  show this help text.
        -r  reference, should be a fasta file. If correspondng bwa indices
            do not exist they will be created. (required).
        -i  fastq/a input assembly (required).
        -c  chunk size. Input reads/contigs will be broken into chunks
            prior to alignment, 0 will not chunk (default 100000).
        -C  catalogue errors.
        -t  alignment threads (default: 1).
        -p  output file prefix (default: assm).
        -T  trim consensus to primary alignments of truth to assembly.

By default alignment is performed in chunks (of the assembly) for speed, though
this can be disabled. It can also be useful when faced with structural
variants.


Simulating nanopore reads
-------------------------

The program `simulate_calls` will produce a `.fasta` sequencing containing
sequences with similar error properties to those produced by nanopore
sequencing. This is achieved by using the
`scrappie <https://github.com/nanoporetech/scrappie>`_ python API to simulate
a nanopore ionic current trace and subsequently perform basecalling. The read
accuracy can be tuned by adding additional noise.

.. code-block:: bash

    simulate_calls [-h] fasta ncalls
    
    Simulate basecalls with scrappy.
    
    positional arguments:
      fasta              Source sequence file.
      ncalls             Number of basecalls to produce.
    
    optional arguments:
      -h, --help         show this help message and exit
      --mu MU            mean fragment length.
      --sigma SIGMA      stdv fragment length.
      --noise NOISE      Additional Gaussian noise on signal.
      --threads THREADS  number of worker threads.

Reads are produced as random fragments of the input sequence with a log-normal
length distribution.


Subsampling alignments
----------------------

As it is often useful to benchmark analysis techniques at various coverage
depths of a reference, the program `subsample_bam` implements a greedy method
to attempt uniform sampling of a `.bam` alignment file.


.. code-block:: bash

    subsample bam to uniform or proportional depth [-h] bam depth [depth ...]
    
    positional arguments:
      bam                   input bam file.
      depth                 Target depth.
    
    optional arguments:
      -h, --help            show this help message and exit
      -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                            Output prefix (default: sub_sampled)
      -r REGIONS [REGIONS ...], --regions REGIONS [REGIONS ...]
                            Only process given regions. (default: None)
      -p PROFILE, --profile PROFILE
                            Stride in genomic coordinates for depth profile.
                            (default: 1000)
      -O {fwd,rev}, --orientation {fwd,rev}
                            Sample only forward or reverse reads. (default: None)
      -t THREADS, --threads THREADS
                            Number of threads to use. (default: -1)
      -q QUALITY, --quality QUALITY
                            Filter reads by mean qscore. (default: None)
      -a ACCURACY, --accuracy ACCURACY
                            Filter reads by accuracy. (default: None)
      -c COVERAGE, --coverage COVERAGE
                            Filter reads by coverage (what fraction of the read
                            aligns). (default: None)
    
    Uniform sampling options:
      -x PATIENCE, --patience PATIENCE
                            Maximum iterations with no change in median coverage
                            before aborting. (default: 5)
      -s STRIDE, --stride STRIDE
                            Stride in genomic coordinates when searching for new
                            reads. Smaller can lead to more compact pileup.
                            (default: 1000)
    
    Proportional sampling options:
      -P, --proportional    Activate proportional sampling, thus keeping depth
                            variations of the pileup. (default: False)
      -S SEED, --seed SEED  Random seed for proportional downsampling of reads.
