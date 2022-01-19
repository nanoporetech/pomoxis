# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

### Unreleased
### Changed
- Include mapping quality in `stats_from_bam` output.
### Added
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
