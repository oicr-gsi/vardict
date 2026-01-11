# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.3] - 2026-01-10
### Fixed
- [GRD-1020](https://jira.oicr.on.ca/browse/GRD-1020) Reordered BAM input in `runVardict` task to `tumor|normal` as required by VarDict to ensure correct analysis (tumor first, normal second).

### Changed
- Updated module loading for `mergeVcfs` task to use `bcftools/1.9 tabix/1.9` instead of `vcftools`.
- Updated `mergeVcfs` task to use `bcftools concat` instead of `vcf-concat` and filter final VCF to PASS and Somatic variants only.
- Adjusted default VarDict parameters: `AF_THR=0.03` and `READ_POSITION_FILTER=8` to improve variant calling sensitivity.

### Added
- Added `-p` flag to `mkdir` in `splitBedByChromosome` task to avoid errors if the `split_beds` directory already exists.


## [1.0.2] - 2025-06-10
### Changed
- [GRD-833](https://jira.oicr.on.ca/browse/GRD-833) Changed the order of input files to ensure 
the correct order of data columns in vcf (NORMAL/TUMOR), satisfying the requirements of Neoantigen pipeline

## [1.0.1] - 2025-04-01
### Added
- Added bam index as input of workflow

## [1.0.0] - 2025-03-10
### Added
- [GRD-833](https://jira.oicr.on.ca/browse/GRD-833), first verion of the wdl along with README and vidarr files

