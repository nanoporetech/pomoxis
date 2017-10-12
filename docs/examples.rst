Pomoxis Examples
================
Below you will find examples of some key tools and ways in which they may be
usefully combined.

Read Alignment with bwa
-----------------------

Pomoxis bundles bwa and samtools to accomplish fastq to bam conversion in a
streamlined manner with the `bwa_align` tool:

.. code-block:: bash

    bwa_align [-h] -r <reference> -i <fastx>
    
    Align fastq/a formatted reads to a genome using bwa.
    
        -h  show this help text.
        -r  reference, should be a fasta file. If correspondng bwa indices
            do not exist they will be created. (required).
        -i  fastq/a input reads (required).
        -a  aggresively extend gaps (sets -A1 -B2 -O2 -E1 for bwa mem).
        -c  chunk size. Input reads/contigs will be broken into chunks
            prior to alignment.
        -t  alignment threads (default: 1).
        -p  output file prefix (default: reads).
    -i and -r must be specified.

Read Alignment with minimap2
----------------------------

.. code-block:: bash

    mini_align [-h] -r <reference> -i <fastx>
    
    Align fastq/a formatted reads to a genome using minimap2.
    
        -h  show this help text.
        -r  reference, should be a fasta file. If correspondng minimap2 indices
            do not exist they will be created. (required).
        -i  fastq/a input reads (required).
        -a  aggresively extend gaps (sets -A1 -B2 -O2 -E1 for minimap2).
        -P  filter to only primary alignments (i.e. run samtools view -F 2308)
        -c  chunk size. Input reads/contigs will be broken into chunks
            prior to alignment.
        -t  alignment threads (default: 1).
        -p  output file prefix (default: reads).
    -i and -r must be specified.

Fast de-novo assembly (or reference-guided consensus)
-----------------------------------------------------

From a .fastq file containing basecalls one can perform an assembly and
consensus using the `mini_assemble` tool.

.. code-block:: bash

    mini_assemble [-h] -i <fastx>
    
    Assemble fastq/fasta formatted reads and perform POA consensus.
    
        -h  show this help text.
        -i  fastx input reads (required).
        -q  use qualities as is (default: false).
        -r  reference fasta for reference-guided consensus (instead of de novo assembly)
        -o  output folder (default: assm).
        -p  output file prefix (default: reads).
        -t  number of minimap and racon threads (default: 1).
    -i must be specified.

For E.coli on a computer with 16 CPUs this should take around 5 minutes for
a 50-fold coverage dataset.

When a reference fasta is provided, a reference-guided consensus is generated instead
of a de novo assembly.

Distributed basecalling
-----------------------

A distributed basecalling platform is included in the epi3me program. This
comprises a `dealer` and a `router` component. The former wil watch for .fast5
files appearing a a given path; numerous instances of the latter can be started
(on networked machines) to perform the basecalling. A current restriction is
that all machines used for basecalling must shared a filesystem.

To start watching a filesytem for new .fast5 files use the dealer component:

.. code-block:: bash

    usage: epi3me dealer [-h] [--port PORT] [--output OUTPUT] path outpath
    
    positional arguments:
      path             Path to watch for new files.
      outpath          Path to move finished files.
    
    optional arguments:
      -h, --help       show this help message and exit
      --port PORT      port on which to run.
      --output OUTPUT  output fasta file.

Once the dealer is running and picking up new files, one can start running
router jobs on the same or different machine:

.. code-block:: bash

    usage: epi3me router [-h] [--addr ADDR]
    
    optional arguments:
      -h, --help   show this help message and exit
      --addr ADDR  Address to use, should include port.

The platform was used to basecall the CliveOME2 dataset on an SGE cluster
within AWS EC2 with files being synced concurrently from AWS S3.
