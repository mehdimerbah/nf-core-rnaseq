#!/usr/bin/env bash
set -euo pipefail

profile="${NXF_PROFILE:-test,docker}"
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"

run_case() {
    local name="$1"
    shift
    nextflow run main.nf -profile "$profile" \
        --outdir "tests/results-$name" \
        -work-dir "tests/work-$name" \
        -ansi-log false \
        "$@"
    test -s "tests/results-$name/samtools/mapped.sorted.bam"
    test -s "tests/results-$name/multiqc/multiqc_report.html"
}

run_case star --aligner star

run_case bowtie2-deseq2 \
    --aligner bowtie2 \
    --run_deseq2 \
    --deseq2_samplesheet validation/small/deseq2_samplesheet.csv \
    --deseq2_counts validation/small/counts_for_deseq2.tsv \
    --deseq2_contrast_variable condition \
    --deseq2_reference WT \
    --deseq2_target KO \
    --deseq2_vst_nsub 10
test -s tests/results-bowtie2-deseq2/deseq2/KO_vs_WT.deseq2.results.tsv

run_case hisat2-featurecounts \
    --aligner hisat2 \
    --run_featurecounts \
    --run_featurecounts_exon
test -s tests/results-hisat2-featurecounts/featurecounts/mapped.gene.featureCounts.tsv
test -s tests/results-hisat2-featurecounts/featurecounts/mapped.exon.featureCounts.tsv
