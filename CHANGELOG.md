# Changelog

The project follows [Semantic Versioning](https://semver.org/).

## 2.0.0 - 2026-07-12

### Added

- Selectable STAR, Bowtie2, and HISAT2 alignment.
- Shared SAMtools sorting, indexing, and statistics.
- Optional gene/exon featureCounts and prepared-matrix DESeq2 analysis.
- Self-contained regression fixtures and public CI.

### Changed

- Rebranded as the independent `mehdimerbah/nextflow-rnaseq` portfolio pipeline.
- Restricted repeated-lane samples to STAR; Bowtie2 and HISAT2 now reject them explicitly.

### Removed

- Unused class-era FastQC, Picard, StringTie, IGV, bedtools, local aggregation modules, profiles, and nf-core infrastructure workflows.
