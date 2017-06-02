Pomoxis Examples
================
Below you will find examples of some key tools and ways in which they may be
usefully combined.

Read Alignment with bwa
-----------------------

Pomoxis bundles bwa and samtools to accomplish fastq to bam conversion in a
streamlined manner with the `bwa_align` tool:

.. code-block:: bash

    bwa_align [-h] [-r reference] [-i fastq] [-p prefix]

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


Fast de-novo assembly
---------------------

From a .fastq file containing basecalls one can perform an assembly and
consensus using the `mini_assemble` tool.

.. code-block:: bash

    mini_assemble [-h] [-i fastq] [-o directory] [-p prefix]
    
    Assemble fastq formatted reads and perform POA consensus
    
        -h  show this help text
        -i  fastq input reads (required)
        -o  output folder (default: assm)
        -p  output file prefix (default: reads)
    -i must be specified.

For E.coli on a computer with 16 CPUs this should take around 5 minutes for
a 50-fold coverage dataset.
