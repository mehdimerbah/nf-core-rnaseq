/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC as FASTQC; FASTQC as FASTQC_TRIMMED              } from '../modules/nf-core/fastqc/main'                // QC       // QC on trimmed reads
include { MULTIQC as MULTIQC; MULTIQC as MULTIQC_TRIMMED        } from '../modules/nf-core/multiqc/main'               // MULTIQC // MULTIQC on trimmed reads
include { TRIMGALORE             } from '../modules/nf-core/trimgalore/main'
include { STAR_GENOMEGENERATE    } from '../modules/nf-core/star/genomegenerate/main'
include { STAR_ALIGN             } from '../modules/nf-core/star/align/main'
include { BOWTIE2_BUILD          } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN          } from '../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_SORT          } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX         } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FAIDX         } from '../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_STATS         } from '../modules/nf-core/samtools/stats/main'
include { INDEXFASTA             } from '../modules/local/indexfasta/indexfasta.nf'
include { PICARD_MARKDUPLICATES  } from '../modules/nf-core/picard/markduplicates/main'
include { STRINGTIE_STRINGTIE    } from '../modules/nf-core/stringtie/stringtie/main'
include { AGGREGATESTRINGTIE     } from '../modules/local/aggregatestringtie/aggregatestringtie.nf'
include { MERGESTRINGTIE         } from '../modules/local/mergestringtie/mergestringtie.nf'
include { IGVREPORTS             } from '../modules/nf-core/igvreports/main'
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
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Module: TRIMGALORE
    //
    TRIMGALORE (
        ch_samplesheet
    )

    TRIMGALORE.out.reads
            .set { ch_samplesheet_trimmed }

    ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.log.collect{it[1]})

    //
    // MODULE: Run FastQC on TRIMMED
    //
    FASTQC_TRIMMED (
        ch_samplesheet_trimmed
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC_TRIMMED.out.versions.first())

    //
    // Module: STAR_GENOMEGENERATE & STAR_ALIGN or BOWTIE2_BUILD & BOWTIE2_ALIGN
    //
    ch_fasta = Channel.of( [ [ id: "${params.igenomes_reference}" ], [params.genomes[params.igenomes_reference].fasta] ] )
    ch_gtf = Channel.of( [ [ id: "${params.igenomes_reference}" ], [params.genomes[params.igenomes_reference].gtf] ] )

    // We could choose the alignment method based on params.aligner
    if (params.aligner == 'star') {
        STAR_GENOMEGENERATE (
            ch_fasta,
            ch_gtf
        )

        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions.first())

        STAR_ALIGN (
            ch_samplesheet_trimmed,
            STAR_GENOMEGENERATE.out.index.collect(),
            ch_gtf.collect(),
            [],
            [],
            []
        )

        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())
        ch_alignment_bam = STAR_ALIGN.out.bam
        ch_alignment_log = STAR_ALIGN.out.log_final

    } else if (params.aligner == 'bowtie2') {
        BOWTIE2_BUILD (
            ch_fasta
        )

        ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions.first())

        BOWTIE2_ALIGN (
            ch_samplesheet_trimmed,
            BOWTIE2_BUILD.out.index.collect(),
            ch_fasta.collect(),
            false,  // This is to specify NOT saving unaligned reads
            true    // sort_bam argment to specify outputing sorted BAM with .sorted prefix
        )

        ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())
        ch_alignment_bam = BOWTIE2_ALIGN.out.bam
        ch_alignment_log = BOWTIE2_ALIGN.out.log
        ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log.collect{it[1]})
    }

    //
    // Module: SAMTOOLS sort, index, stats
    //
    SAMTOOLS_SORT (
            ch_alignment_bam,
            ch_fasta.collect(),
            []
        )

    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_SORT.out.bam
        .set { ch_bam }

    SAMTOOLS_INDEX (
        ch_bam
    )

    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

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

    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.collect{it[1]})

    //
    // Module: IGV Reports - Generate IGV visualization reports
    //
    // Prepare FASTA with index for IGV
    SAMTOOLS_FAIDX (
        ch_fasta,
        [[], []],
        false  // get_sizes
    )

    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())

    //
    // Module: PICARD markduplicates
    //

    // Local Module IndexFasta: create reference .fai index file
    INDEXFASTA ( ch_fasta )
    INDEXFASTA.out.fai.set { ch_fai }

    PICARD_MARKDUPLICATES (
        ch_bam,
        ch_fasta.collect(),
        ch_fai.collect()
    )

    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_MARKDUPLICATES.out.metrics.collect{it[1]})

    //
    // Module: StringTie
    //
    // gtf_channel = Channel.of(params.genomes[params.igenomes_reference].gtf)

    STRINGTIE_STRINGTIE (
        PICARD_MARKDUPLICATES.out.bam,
        ch_gtf.map { meta, gtf -> gtf}.collect()
    )

    ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions.first())

    //
    // Module: AGGREGATESTRINGTIE - Remove duplicate gene entries
    //
    AGGREGATESTRINGTIE (
        STRINGTIE_STRINGTIE.out.abundance
    )

    ch_versions = ch_versions.mix(AGGREGATESTRINGTIE.out.versions.first())

    //
    // Module: AGGREGATESTRINGTIE & MERGESTRINGTIE
    //
    AGGREGATESTRINGTIE.out.abundance.map{meta, f -> f}.collect().set {merge_in}
    
    MERGESTRINGTIE ( merge_in )
    MERGESTRINGTIE.out.tsv.view { it -> "TPM table exported to $it" }

    ch_versions = ch_versions.mix(MERGESTRINGTIE.out.versions.first())

    ///////////////////////////////////
    ///////////////////////////////////


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'rnaseqpipeline_software_'  + 'mqc_'  + 'versions.yml',
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
