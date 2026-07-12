# Local module provenance

These modules started from nf-core/modules components but contain pipeline-specific changes, so they are deliberately kept outside `modules/nf-core`.

- `hisat2/align` and `hisat2/build`: based on nf-core/modules `41dfa3f`; expose pipeline-specific read-group, unaligned-read, and index-memory behavior.
- `subread/featurecounts`: based on nf-core/modules `41dfa3f`; adds `--countReadPairs` for paired-end Subread 2.0.6 runs.
- `multiqc`: based on nf-core/modules `41dfa3f`; constrains Python below 3.13 for the pinned MultiQC 1.29 Conda environment.

Untouched upstream modules remain under `modules/nf-core` and are recorded in `modules.json`.
