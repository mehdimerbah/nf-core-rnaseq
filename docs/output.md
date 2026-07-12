# mehdimerbah/nextflow-rnaseq: Output

## Introduction

This document describes the outputs produced by the current mehdimerbah/nextflow-rnaseq workflow. The pipeline trims FASTQ reads, aligns them with one selected aligner, normalizes the alignment output through SAMtools, optionally runs featureCounts and DESeq2, and summarizes available logs with MultiQC.

All paths below are relative to the top-level output directory.

## Pipeline Overview

The active workflow runs these stages:

- [TrimGalore](#trimgalore) - adapter and quality trimming
- [STAR](#star) or [Bowtie2](#bowtie2) or [HISAT2](#hisat2) - selected with `--aligner`
- [SAMtools](#samtools) - sorted/indexed BAMs and alignment statistics
- [featureCounts](#featurecounts) - optional gene-level and opt-in exon-level counting
- [DESeq2](#deseq2) - optional differential expression analysis
- [MultiQC](#multiqc) - aggregate report
- [Pipeline Information](#pipeline-information) - execution metadata

## TrimGalore

<details><summary>Output files</summary>

- `trimgalore/*_val_1.fq.gz` and `trimgalore/*_val_2.fq.gz` - trimmed paired-end reads
- `trimgalore/*_trimmed.fq.gz` - trimmed single-end reads
- `trimgalore/*_trimming_report.txt` - trimming and adapter-removal summary

</details>

TrimGalore wraps Cutadapt and trims adapters and low-quality bases before alignment. The trimming reports are included in MultiQC.

## STAR

Produced when `--aligner star` is selected.

<details><summary>Output files</summary>

- `star/star/` - generated STAR genome index
- `star/*.Aligned.sortedByCoord.out.bam` - STAR alignment output
- `star/*Log.final.out`, `star/*Log.out`, `star/*Log.progress.out` - STAR logs
- `star/*.SJ.out.tab` - splice junction table

</details>

STAR is splice-aware and is the default aligner.

## Bowtie2

Produced when `--aligner bowtie2` is selected.

<details><summary>Output files</summary>

- `bowtie2/bowtie2/` - generated Bowtie2 genome index
- `bowtie2/*.bam` - Bowtie2 alignment output before the shared SAMtools sort step
- `bowtie2/*.bowtie2.log` - Bowtie2 alignment logs

</details>

Bowtie2 is not splice-aware, but is useful for quick mapping or comparison runs.

## HISAT2

Produced when `--aligner hisat2` is selected.

<details><summary>Output files</summary>

- `hisat2/hisat2/` - generated HISAT2 index when `--hisat2_index` is not supplied
- `hisat2/*.bam` - HISAT2 alignment output before the shared SAMtools sort step
- `hisat2/*.hisat2.summary.log` - HISAT2 alignment summary
- `hisat2/*.unmapped_1.fastq.gz` and `hisat2/*.unmapped_2.fastq.gz` - unmapped reads when `--save_unaligned` is enabled

</details>

HISAT2 is splice-aware. A pre-built index directory can be supplied with `--hisat2_index`.

## SAMtools

<details><summary>Output files</summary>

- `samtools/*.sorted.bam` - coordinate-sorted BAM files from the selected aligner
- `samtools/*.sorted.bam.bai` - BAM indexes
- `samtools/*.stats` - alignment statistics

</details>

All aligner branches feed into the same SAMtools sort, index, and stats steps so downstream outputs are consistent.

## featureCounts

Gene-level output is produced only when `--run_featurecounts` is enabled. Feature-level exon output is produced only when `--run_featurecounts_exon` is also enabled. Requires `--gtf`; convert GFF/GFF3 annotations to GTF before counting. Exon rows are identified by gene ID plus genomic coordinates and are not a complete alternative-splicing analysis.

<details><summary>Output files</summary>

- `featurecounts/*.gene.featureCounts.tsv` - gene-level counts
- `featurecounts/*.gene.featureCounts.tsv.summary` - gene-level assignment summary
- `featurecounts/*.exon.featureCounts.tsv` - exon-level counts, when `--run_featurecounts_exon` is enabled
- `featurecounts/*.exon.featureCounts.tsv.summary` - exon-level assignment summary, when `--run_featurecounts_exon` is enabled

</details>

By default, featureCounts uses `--primary`. With `--fc_count_multimappers`, it uses `-M -O` to count multi-mapping and overlapping reads.

## DESeq2

Produced only when `--run_deseq2` is enabled. Requires `--deseq2_samplesheet` and `--deseq2_counts`.

<details><summary>Output files</summary>

- `deseq2/*.deseq2.results.tsv` - differential expression results
- `deseq2/*.normalised_counts.tsv` - normalized counts
- `deseq2/*.deseq2.sizefactors.tsv` - size factors
- `deseq2/*.deseq2.dispersion.png` - dispersion plot
- `deseq2/*.deseq2.model.txt` - model information
- `deseq2/*.R_sessionInfo.log` - R session information

</details>

The current DESeq2 step consumes a prepared count matrix. It does not yet merge per-sample featureCounts outputs into a DESeq2 matrix inside the workflow.
The DESeq2 count matrix gene column defaults to `gene_id`, and the sample metadata ID column defaults to `experiment_accession`.

## MultiQC

<details><summary>Output files</summary>

- `multiqc/multiqc_report.html` - aggregate HTML report
- `multiqc/multiqc_data/` - parsed report data

</details>

MultiQC summarizes available TrimGalore, aligner, SAMtools, featureCounts, workflow, and software-version outputs.

## Pipeline Information

<details><summary>Output files</summary>

- `pipeline_info/execution_report_*.html` - Nextflow execution report
- `pipeline_info/execution_timeline_*.html` - task timeline
- `pipeline_info/execution_trace_*.txt` - task trace
- `pipeline_info/pipeline_dag_*.html` - workflow DAG
- `pipeline_info/nextflow_rnaseq_software_mqc_versions.yml` - software versions used by MultiQC
- `pipeline_info/params_*.json` - run parameters

</details>

These files are useful for debugging, reproducibility, and resource-usage review.
