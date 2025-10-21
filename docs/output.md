# nf-core/rnaseqpipeline: Output

## Introduction

This document describes the output produced by the nf-core/rnaseqpipeline. The pipeline performs comprehensive RNA-seq analysis including quality control, read trimming, alignment, duplicate removal, and transcript quantification. Most of the plots and summary statistics are consolidated in the MultiQC report, which provides an overview of results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes RNA-seq data using the following steps:

- [FastQC](#fastqc) - Raw read quality control
- [TrimGalore](#trimgalore) - Adapter trimming and quality filtering
- [FastQC (Trimmed)](#fastqc-trimmed) - Quality control on trimmed reads
- [STAR](#star-alignment) - Genome indexing and read alignment
- [SAMtools](#samtools) - BAM file processing and statistics
- [Picard MarkDuplicates](#picard-markduplicates) - PCR duplicate removal
- [StringTie](#stringtie-assembly--quantification) - Transcript assembly and quantification
- [StringTie Aggregation](#stringtie-aggregation) - Merge and aggregate expression matrices
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

The pipeline follows a standard RNA-seq analysis workflow: quality control → trimming → alignment → duplicate removal → quantification → reporting.

> **Note:** Exact directories depend on the modules you enable and the parameters used.

---

### FastQC

<details><summary>Output files</summary>

- `fastqc/*_fastqc.html` – interactive read QC per input FASTQ  
- `fastqc/*_fastqc.zip` – tabular metrics + images

</details>

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) performs quality control checks on raw sequence data coming from high throughput sequencing pipelines.

**What to check:** per-base quality (aim ≥ Q30), adapter content, GC distribution, overrepresented sequences. (Summaries appear in MultiQC.)

---

### TrimGalore

<details><summary>Output files</summary>

- `trimgalore/<sample>_trimmed.fq.gz` – trimmed FASTQ files
- `trimgalore/<sample>_trimming_report.txt` – trimming statistics and adapter removal summary

</details>

[TrimGalore](https://github.com/FelixKrueger/TrimGalore) is a wrapper around [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) and [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) that consistently applies quality and adapter trimming to FastQ files.

**Why it matters:** Removes low-quality bases and adapter sequences that can negatively impact downstream alignment and analysis. The trimming reports show how many reads were processed and what adapters were detected.

---

### FastQC (Trimmed)

<details><summary>Output files</summary>

- `fastqc/*_trimmed_fastqc.html` – interactive QC report for trimmed reads
- `fastqc/*_trimmed_fastqc.zip` – tabular metrics + images for trimmed reads

</details>

Quality control analysis performed on the trimmed reads to ensure the trimming process was effective.

**What to check:** Improvement in quality scores compared to raw reads, successful removal of adapter contamination, and overall read quality distribution.

---

### STAR alignment

<details><summary>Output files</summary>

- `star/<sample>.bam` (+ `.bai`) – coordinate-sorted alignments  
- `star/<sample>.SJ.out.tab` – detected splice junctions  
- `star/Log.final.out` and other `Log.*` – mapping summary & parameters

</details>

[STAR](https://github.com/alexdobin/STAR) (Spliced Transcripts Alignment to a Reference) is a fast RNA-seq read aligner designed to handle splice junctions in eukaryotic transcripts.

**Why it matters:** BAMs underpin expression estimates and allow IGV inspection; `Log.final.out` (parsed by MultiQC) shows % mapped/uniquely mapped, mismatch rates, chimeric fraction. STAR is particularly well-suited for RNA-seq as it can align reads across exon-exon junctions.

---

### SAMtools

<details><summary>Output files</summary>

- `samtools/<sample>.sorted.bam` (+ `.bai`) – coordinate-sorted BAM files
- `samtools/<sample>.stats.txt` – comprehensive alignment statistics
- `samtools/<sample>.flagstat.txt` – basic alignment statistics

</details>

[SAMtools](http://www.htslib.org/) is a suite of programs for interacting with high-throughput sequencing data in SAM/BAM format.

**Why it matters:** 
- **Sort**: Coordinates sorting enables efficient downstream processing and visualization
- **Index**: Creates index files required for fast random access to BAM files
- **Stats**: Provides detailed statistics about alignment quality, insert sizes, and coverage that are essential for QC

---

### Picard MarkDuplicates

<details><summary>Output files</summary>

- `picard/<sample>.markdup.bam` (+ `.bai`) – duplicates flagged  
- `picard/<sample>.markdup.metrics.txt` – duplication summary

</details>

[Picard MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard) identifies and flags PCR duplicate reads in BAM files.

**Why it matters:** PCR duplicates can artificially inflate expression levels and bias quantification results. Duplicate rates provide insight into library complexity and can influence differential expression analysis strategies. MultiQC collates these metrics for easy interpretation across samples.

---

### StringTie assembly & quantification

<details><summary>Output files</summary>

- `stringtie/<sample>.gtf` – per-sample assembled transcripts  
- `stringtie/<sample>_abundance.tsv` – transcript-level TPM/FPKM  

</details>

[StringTie](https://ccb.jhu.edu/software/stringtie/) assembles RNA-seq alignments into potential transcripts and estimates their abundances.

**Why it matters:** Enables gene/transcript quantification and potential novel isoform discovery. StringTie can identify new splice variants and provide accurate abundance estimates in TPM (Transcripts Per Million) and FPKM (Fragments Per Kilobase Million) units.

---

### StringTie Aggregation

<details><summary>Output files</summary>

- `gene_table/gene_table_TPM.tsv` – merged gene-level TPM expression matrix

</details>

Custom aggregation modules (`AGGREGATESTRINGTIE` and `MERGESTRINGTIE`) process individual StringTie outputs to create merged expression matrices.

**Why it matters:** 
- **Aggregation**: Removes duplicate gene entries and standardizes identifiers across samples
- **Merging**: Combines all samples into ready-to-use expression matrices for downstream analysis (e.g., DESeq2, edgeR)
- **Output formats**: Both gene-level and transcript-level quantifications are provided in standard formats

---

### MultiQC

<details><summary>Output files</summary>

- `multiqc/multiqc_report.html` – comprehensive HTML summary report  
- `multiqc/multiqc_data/` – parsed statistics in machine-readable formats
- `multiqc/multiqc_plots/` – individual plot files in multiple formats

</details>

[MultiQC](http://multiqc.info/) searches a given directory for analysis logs and compiles an HTML report with interactive plots and tables summarizing key metrics from the entire pipeline.

**What's included:** 
- **FastQC summaries**: Quality metrics for raw and trimmed reads
- **TrimGalore reports**: Adapter removal and trimming statistics  
- **STAR alignment**: Mapping rates, splice junction detection, and alignment metrics
- **SAMtools stats**: Detailed alignment statistics and quality metrics
- **Picard metrics**: PCR duplication rates and library complexity
- **Software versions table**: Complete list of tools and versions used
- **Workflow summary**: Pipeline parameters and execution details

**Why it's essential:** Provides a single, comprehensive view of your entire analysis, making it easy to identify potential issues, compare samples, and generate figures for publications and reports.

---


### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline.

**Key reports:**
- **Execution report**: Resource usage, task duration, and success/failure status for each process
- **Timeline**: Visual timeline showing when each task ran and how long it took
- **Trace file**: Detailed execution log with resource consumption metrics
- **Pipeline DAG**: Directed acyclic graph showing the workflow structure and dependencies

These reports allow you to:
- Troubleshoot errors and failed processes
- Optimize resource allocation for future runs  
- Track computational requirements and costs
- Document methods and software versions for reproducibility
- Monitor pipeline performance and identify bottlenecks
