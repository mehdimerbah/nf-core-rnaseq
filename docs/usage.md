
> Exact directories depend on the modules you enable.

---

## FastQC

<details><summary>Output files</summary>

- `fastqc/*_fastqc.html` – interactive read QC per input FASTQ  
- `fastqc/*_fastqc.zip` – tabular metrics + images

</details>

**What to check:** per-base quality (aim ≥ Q30), adapter content, GC distribution, overrepresented sequences. (Summaries appear in MultiQC.)

---

## STAR alignment *(if enabled)*

<details><summary>Output files</summary>

- `star/<sample>.bam` (+ `.bai`) – coordinate-sorted alignments  
- `star/<sample>.SJ.out.tab` – detected splice junctions  
- `star/Log.final.out` and other `Log.*` – mapping summary & parameters

</details>

**Why it matters:** BAMs underpin expression estimates and allow IGV inspection; `Log.final.out` (parsed by MultiQC) shows % mapped/uniquely mapped, mismatch rates, chimeric fraction.

---

## Picard MarkDuplicates 

<details><summary>Output files</summary>

- `picard/<sample>.markdup.bam` (+ `.bai`) – duplicates flagged  
- `picard/<sample>.markdup.metrics.txt` – duplication summary

</details>

**Why it matters:** duplicate rates contextualise library complexity and can influence DE analysis; MultiQC collates these metrics.

---

## StringTie assembly & quantification 

<details><summary>Output files</summary>

- `stringtie/<sample>.gtf` – per-sample assembled transcripts  
- `stringtie/<sample>_abundance.tsv` – transcript-level TPM/FPKM  
- `gene_table/` – merged matrices after merge/aggregation (e.g., `gene_table_TPM.tsv`)

</details>

**Why it matters:** enables gene/transcript quantification and potential novel isoform discovery; merged tables are ready for downstream analysis (e.g., DESeq2 / edgeR).

---

## MultiQC

<details><summary>Output files</summary>

- `multiqc/multiqc_report.html` – single summary report  
- `multiqc/multiqc_data/` – parsed stats for programmatic reuse  
- `multiqc/multiqc_plots/` – static plot assets

</details>

**What’s inside:** FastQC summaries, alignment & duplication barplots, and a **Software Versions** table derived from the pipeline — ideal for methods and reviews.

---

## Pipeline information (provenance & reproducibility)

<details><summary>Output files</summary>

- `pipeline_info/execution_report.html` – overview of all tasks (CPU, RAM, times)  
- `pipeline_info/execution_timeline.html` – Gantt chart of the run  
- `pipeline_info/execution_trace.txt` – per-task table (grep-able)  
- `pipeline_info/pipeline_dag.svg` – workflow graph  
- `pipeline_info/software_versions.yml` – exact tool versions  
- `pipeline_info/samplesheet.valid.csv` – schema-validated input sheet  
- `pipeline_info/params.json` – all runtime parameters (freeze for exact re-runs)

</details>

**Why it matters:** these artefacts make your run **auditable and reproducible**; include MultiQC, `software_versions.yml`, and `params.json` when sharing results.

