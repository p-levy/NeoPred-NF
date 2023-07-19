# Akazhiel/NeoPred-NF: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0

### Added

- `T1K` for DNA HLA typing aside from OptiType.
- Now proximal variants are reported on a per variant basis using `pysam`

### Changed

- `yaramap` now uses `samtools collate` and `samtools fastq` to generate the intermediate gzipped fastq to feed into yara mapper. Intermediate files are not removed.
- `VEP` now generates a bgzipped and tabix indexed VCF.

### Fixed

### Dependencies

### Deprecated

### Removed


## v0.9.0

### Added

- `Fastp` for trimming and quality control of input FASTQ.
- `arcasHLA` added for typing and expression of HLA alleles in RNA.

### Changed

- `trim_galore` and `FastQC` have been replaced with `fastp`
- `OptiType` has been removed for `neopredrna` and replaced with `arcasHLA`
- `VEP` and `pyensembl` bumped from version 104 to 109.

### Fixed

### Dependencies

### Deprecated

### Removed
