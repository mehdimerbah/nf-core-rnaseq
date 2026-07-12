/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/local/multiqc/main'
include { TRIMGALORE             } from '../modules/nf-core/trimgalore/main'
include { STAR_GENOMEGENERATE    } from '../modules/nf-core/star/genomegenerate/main'
include { STAR_ALIGN             } from '../modules/nf-core/star/align/main'
include { BOWTIE2_BUILD          } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN          } from '../modules/nf-core/bowtie2/align/main'
include { HISAT2_EXTRACTSPLICESITES } from '../modules/nf-core/hisat2/extractsplicesites/main'
include { HISAT2_BUILD           } from '../modules/local/hisat2/build/main'
include { HISAT2_ALIGN           } from '../modules/local/hisat2/align/main'
include { SAMTOOLS_SORT          } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX         } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS         } from '../modules/nf-core/samtools/stats/main'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_GENE  } from '../modules/local/subread/featurecounts/main'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_EXON  } from '../modules/local/subread/featurecounts/main'
include { DESEQ2_DIFFERENTIAL                          } from '../modules/nf-core/deseq2/differential/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_rnaseqpipeline_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASEQPIPELINE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:


    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Module: TRIMGALORE
    //
    TRIMGALORE (
        ch_samplesheet
    )

    TRIMGALORE.out.reads
            .set { ch_samplesheet_trimmed }

    ch_versions = ch_versions.mix(TRIMGALORE.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.log.collect{it[1]})

    //
    // Module: STAR_GENOMEGENERATE & STAR_ALIGN or BOWTIE2_BUILD & BOWTIE2_ALIGN or HISAT2_BUILD & HISAT2_ALIGN
    //
    // Choose either a fully custom reference or one configured through iGenomes.
    // All aligner branches emit BAMs that are normalized through SAMTOOLS_SORT.
    if ((params.fasta && !params.gtf) || (!params.fasta && params.gtf)) {
        error "Custom references must be provided as a pair. Set both --fasta and --gtf, or neither and use --genome/--igenomes_reference."
    }

    def reference_id = params.fasta ? 'custom_genome' : (params.genome ?: params.igenomes_reference)
    if (params.fasta && params.gtf) {
        ch_fasta = Channel.value( [ [ id: reference_id ], file(params.fasta) ] )
        ch_gtf   = Channel.value( [ [ id: reference_id ], file(params.gtf) ] )
    } else {
        if (!params.genomes || !reference_id || !params.genomes.containsKey(reference_id)) {
            error "No usable reference found. Provide --fasta and --gtf, or set --igenomes_reference/--genome to a key in conf/igenomes.config."
        }
        ch_fasta = Channel.value( [ [ id: reference_id ], file(params.genomes[reference_id].fasta) ] )
        ch_gtf   = Channel.value( [ [ id: reference_id ], file(params.genomes[reference_id].gtf) ] )
    }

    // Choose the alignment method based on params.aligner.
    if (params.aligner == 'star') {
        STAR_GENOMEGENERATE (
            ch_fasta,
            ch_gtf
        )

        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)

        STAR_ALIGN (
            ch_samplesheet_trimmed,
            STAR_GENOMEGENERATE.out.index.collect(),
            ch_gtf.collect(),
            false,
            'ILLUMINA',
            params.seq_center ?: ''
        )

        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)
        ch_alignment_bam = STAR_ALIGN.out.bam
        ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_final.collect{it[1]})

    } else if (params.aligner == 'bowtie2') {
        BOWTIE2_BUILD (
            ch_fasta
        )

        ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)

        BOWTIE2_ALIGN (
            ch_samplesheet_trimmed,
            BOWTIE2_BUILD.out.index.collect(),
            ch_fasta.collect(),
            false,  // This is to specify NOT saving unaligned reads
            false   // keep a consistent downstream sort/index/stats path for all aligners
        )

        ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)
        ch_alignment_bam = BOWTIE2_ALIGN.out.bam
        ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log.collect{it[1]})

    } else if (params.aligner == 'hisat2') {
        // Use pre-built index if provided, otherwise build it
        if (params.hisat2_index) {
            // Using pre-built HISAT2 index - use .collect() to make it a value channel
            // that can be reused for all samples
            ch_hisat2_index = Channel.value([[id: 'hisat2_index'], file(params.hisat2_index)])

            HISAT2_ALIGN (
                ch_samplesheet_trimmed,
                ch_hisat2_index,
                Channel.value([[id: 'none'], []])  // Empty splice sites channel as value
            )
        } else {
            // Build HISAT2 index from scratch
            HISAT2_EXTRACTSPLICESITES (
                ch_gtf
            )

            ch_versions = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)

            HISAT2_BUILD (
                ch_fasta,
                ch_gtf,
                HISAT2_EXTRACTSPLICESITES.out.txt
            )

            ch_versions = ch_versions.mix(HISAT2_BUILD.out.versions)

            HISAT2_ALIGN (
                ch_samplesheet_trimmed,
                HISAT2_BUILD.out.index.collect(),
                HISAT2_EXTRACTSPLICESITES.out.txt.collect()
            )
        }

        ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions)
        ch_alignment_bam = HISAT2_ALIGN.out.bam
        ch_multiqc_files = ch_multiqc_files.mix(HISAT2_ALIGN.out.summary.collect{it[1]})
    } else {
        error "Unsupported aligner '${params.aligner}'. Choose one of: star, bowtie2, hisat2."
    }

    //
    // Module: SAMTOOLS sort, index, stats
    //
    SAMTOOLS_SORT (
            ch_alignment_bam,
            ch_fasta.collect(),
            []
        )

    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    SAMTOOLS_SORT.out.bam
        .set { ch_bam }

    SAMTOOLS_INDEX (
        ch_bam
    )

    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    ch_bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }
        .set { ch_bam_bai }

    SAMTOOLS_STATS (
        ch_bam_bai,
        ch_fasta.collect()
    )

    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.collect{it[1]})

    // ========================================================================
    // FEATURECOUNTS - Read counting for gene expression quantification
    // ========================================================================

    def featurecounts_annotation = params.gtf

    if (params.run_featurecounts_exon && !params.run_featurecounts) {
        error "Exon-level featureCounts requested with --run_featurecounts_exon, but --run_featurecounts is not enabled."
    }

    // Only run featureCounts if requested and an annotation is available.
    if (params.run_featurecounts) {
        if (!featurecounts_annotation) {
            error "featureCounts requested with --run_featurecounts, but no --gtf annotation was provided. Convert GFF/GFF3 annotations to GTF before counting."
        }

        // Prepare input for featureCounts: combine BAM with annotation
        ch_bam
            .map { meta, bam -> [meta, bam, file(featurecounts_annotation)] }
            .set { ch_featurecounts_input }

        //
        // Module: FEATURECOUNTS - Gene-level counting
        // Uses -t exon -g gene_id.
        //
        FEATURECOUNTS_GENE (
            ch_featurecounts_input
        )

        ch_versions = ch_versions.mix(FEATURECOUNTS_GENE.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(FEATURECOUNTS_GENE.out.summary.collect{ _meta, summary -> summary })

        if (params.run_featurecounts_exon) {
            //
            // Module: FEATURECOUNTS - Feature-level exon counting
            // Uses feature-level exon counting.
            //
            FEATURECOUNTS_EXON (
                ch_featurecounts_input
            )

            ch_versions = ch_versions.mix(FEATURECOUNTS_EXON.out.versions)
            // Exon-level tables are published, but only gene-level summaries are sent to
            // MultiQC because both reports otherwise collapse to the same sample name.
        }
    }

    // ========================================================================
    // DESEQ2 - Differential Expression Analysis
    // ========================================================================

    if (params.run_deseq2) {
        if (!params.deseq2_samplesheet || !params.deseq2_counts) {
            error "DESeq2 requested with --run_deseq2, but --deseq2_samplesheet and --deseq2_counts were not both provided."
        }

        // DESeq2 requires:
        // 1. Contrast info: meta, contrast_variable, reference, target, formula, comparison
        // 2. Samplesheet and counts: meta2, samplesheet, counts
        // 3. Control genes (optional): meta3, control_genes_file
        // 4. Transcript lengths (optional): meta4, transcript_lengths_file

        ch_deseq2_contrast = Channel.of([
            [id: params.deseq2_contrast_id ?: 'KO_vs_WT'],  // meta
            params.deseq2_contrast_variable ?: 'condition', // contrast_variable
            params.deseq2_reference ?: 'WT',                // reference
            params.deseq2_target ?: 'KO',                   // target
            '',                                              // formula (empty = use contrast_variable)
            ''                                               // comparison (empty = auto-generate)
        ])

        ch_deseq2_input = Channel.of([
            [id: 'samplesheet'],                            // meta2
            file(params.deseq2_samplesheet),                // samplesheet
            file(params.deseq2_counts)                      // counts
        ])

        ch_deseq2_control_genes = Channel.of([
            [id: 'none'],                                   // meta3
            []                                               // control_genes_file (empty)
        ])

        ch_deseq2_lengths = Channel.of([
            [id: 'none'],                                   // meta4
            []                                               // transcript_lengths_file (empty)
        ])

        DESEQ2_DIFFERENTIAL (
            ch_deseq2_contrast,
            ch_deseq2_input,
            ch_deseq2_control_genes,
            ch_deseq2_lengths
        )

        ch_versions = ch_versions.mix(DESEQ2_DIFFERENTIAL.out.versions)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nextflow_rnaseq_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
    ch_multiqc_files.collect(),
    ch_multiqc_config.toList(),
    ch_multiqc_custom_config.toList(),
    ch_multiqc_logo.toList(),
    [],
    []
)

    emit:
        multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
        versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
