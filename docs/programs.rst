
Pomoxis Programs
================

Below you will find a listing with brief description of all the tools
within pomoxis. For more information concerning each tool simply run
it on the command line with the ``--help`` option.

More complex workflows (using tools from pomoxis) can be found within the
`katuali <https://github.com/nanoporetech/katuali>`_ software package.



assess_assembly
***************

Calculate accuracy statistics for an assembly.

.. code-block:: bash

    assess_assembly [-h] -r <reference> -i <fastq>

catalogue_errors
****************

Create a catalogue of all query errors in a bam.

.. code-block:: bash

    catalogue_errors [-h] [--bed BED] [-t THREADS] [-o OUTDIR] bam

common_errors_from_bam
**********************

Get errors common to multiple assemblies aligned to ref.

.. code-block:: bash

    common_errors_from_bam [-h] [-o OUTPUT_PREFIX] bam ref_fasta

coverage_from_bam
*****************

Calculate read coverage depth from a bam.

.. code-block:: bash

    coverage_from_bam [-h] [-r REGIONS [REGIONS ...]] [-p PREFIX]
                             [-s STRIDE]
                             bam

coverage_from_fastx
*******************

Estimate coverage from summed basecall and reference lengths

.. code-block:: bash

    coverage_from_fastx [-h] [--coverage COVERAGE] [--longest]
                               basecalls ref_len

fast_convert
************

Convert between fasta<->fastq.

.. code-block:: bash

    fast_convert [-h] [--discard_q] [--mock_q MOCK_Q] {qq,aa,aq,qa}

intersect_assembly_errors
*************************

Assess errors which occur in the same reference position accross multiple assemblies.

.. code-block:: bash

    intersect_assembly_errors [-h] -r <reference> -i <fasta>

long_fastx
**********

Extract longest reads from a fastq.

.. code-block:: bash

    long_fastx [-h] (--longest LONGEST | --bases BASES) [--others OTHERS]
                      input output

mini_align
**********

Align fastq/a formatted reads to a genome using minimap2.

.. code-block:: bash

    mini_align [-h] -r <reference> -i <fastq>

mini_assemble
*************

Assemble fastq/fasta formatted reads and perform POA consensus.

.. code-block:: bash

    mini_assemble [-h] -i <fastq>

pomoxis_path
************

Print the path of bundled executables.

.. code-block:: bash

    pomoxis_path [-h] program

qscores_from_summary
********************

Extract Q scores from summary_from_stats output

.. code-block:: bash

    qscores_from_summary [-h] [--median] [--ref REF]
                                summaries [summaries ...]

ref_seqs_from_bam
*****************

Extract reference sequence that queries are aligned to

.. code-block:: bash

    ref_seqs_from_bam [-h] bam

split_fastx
***********

Split records in a fasta/q file into chunks of a maximum size.

.. code-block:: bash

    split_fastx [-h] input output chunksize

stats_from_bam
**************

Parse a bamfile (from a stream) and output summary stats for each read.

.. code-block:: bash

    stats_from_bam [-h] [--bed BED] [-m MIN_LENGTH] [-a] [-o OUTPUT]
                          [-s SUMMARY] [-t THREADS]
                          bam

subsample_bam
*************

Subsample a bam to uniform or proportional depth

.. code-block:: bash

    subsample_bam [-h] [-o OUTPUT_PREFIX] [-r REGIONS [REGIONS ...]]
                         [-p PROFILE] [-O {fwd,rev}] [-t THREADS] [-q QUALITY]
                         [-a ACCURACY] [-c COVERAGE] [--any_fail | --all_fail]
                         [-x PATIENCE] [-s STRIDE] [-P] [-S SEED]
                         bam depth [depth ...]

summary_from_stats
******************

Summarise output of `stats_from_bam`.

.. code-block:: bash

    summary_from_stats [-h] [-i INPUT] [-o OUTPUT]
                              [-p PERCENTILES [PERCENTILES ...]] [-pr]

trim_alignments
***************

Trim alignments in multiple bams to common overlap window.

.. code-block:: bash

    trim_alignments [-h] [-r REF_NAME] [-o OUTPUT_PREFIX]
                           [-f REFERENCE_FASTA]
                           bams [bams ...]

