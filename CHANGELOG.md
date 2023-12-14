# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.3.15] - 2023-12-2023
### Added
- `mini_align` can take a bam as input and optionally retain all or subset of bam tags. 

## [v0.3.14] - 2023-10-25
### Changed
- `subsample_bam` and `coverage_from_bam` now have unified read filtering options and logic.
### Added
- `filter_bam` to filter a bam with the same logic used in `subsample_bam`. 
### Fixed
- `subsample_bam` was previously subsampling proportionally before filtering resulting in lower than expected depth. 

## [v0.3.13] - 2023-06-23
### Changed
- `subsample_bam`: `--force_low_coverage` saves contigs with coverage below the target 
- `subsample_bam`: `--force_non_primary` saves multimapping for the subsampled reads
- `coverage_from_bam`: `--primary_only` considers only primary reads when computing the depth
- `bedtools`: upgraded to v2.31
- `porechop`: switched to using Artic version 
### Added
- Option `-C` for `mini_align` to copy fastx comments into bam tags
### Fixed
- Minor compatibility fixes to support `pandas>=2.0` 

## [v0.3.12] - 2023-02-09
### Changed
- `subsample_bam`: `--quality` filtering now uses mean error probability, not mean of quality scores as previously.
- `subsample_bam`: enable filtering for proportional subsampling.

## [v0.3.11] - 2022-11-16
### Fixed
- Fix crashes in `subsample_bam` with alignment filtering and `common_errors_from_bam`
- `assess_assembly -H` uses correct output directory.
- Handling of comments in bed files.
### Changed
- Added `Q(sub)` to summary output.
- Ported bed file handling from `intervaltrees` to [`ncls`](https://github.com/biocore-ntnu/ncls), speeding up assessment and multithreading efficiency.
### Added

## [v0.3.10] - 2022-02-22
### Fixed
- `stats_from_bam`: handle cigar strings using `=` and `X` instead of `M`.
### Changed
- Include mapping quality in `stats_from_bam` output.
### Added
- Handling of LRA bams in which NM tag is number of matches rather than edit distance. 
- Added an option (`-y`) to `assess_assembly` and `mini_align` to include supplementary alignments. 
- Added an option (`-d`) to `mini_align` and `assess_assembly` to select minimap2 alignment preset.
- Added accumulation of errors over a number of chunks (`-a` option in `summary_from_stats` and `assess_assembly`) to get better stats.
- Use `-L` option for `minimap2`.
- Updated versions of minimap2, samtools, bcftools, bedtools, seqkit in Makefile to the most recent ones.

## [v0.3.9] - 2021-08-18
### Fixed
- Reduced memory consumption of `catalogue_errors`.
- `fast_convert qa` now properly outputs a fasta file
- Fixed `long_fastx` `--others` option
- Fixed `split_fastx` fastq output
### Added
- `assess_homopolymers` can use multiple threads

## [v0.3.8] - 2021-06-22
### Changed
- Install `paftools.js` from minimap2 and `k8`

## [v0.3.7] - 2021-05-10
### Changed
- Speed improvements to several benchmarking and analysis scripts
### Fixed
- Quoted all variables in `mini_align` to handle spaces in inputs.

## [v0.3.6] - 2020-02-17
### Changed
 - `stats_from_bam` no longer throws exception when no alignments have been proceseed.
### Added
 - `coverage_from_bam` now has a `--one_file` option to better specify the output in the common usage.

## [v0.3.5] - 2020-02-01
### Removed
 - Python 3.5 support
### Added
 - Python >3.6 support
