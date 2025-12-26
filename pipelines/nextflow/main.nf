#!/usr/bin/env nextflow

/*
 * yallHap Nextflow Pipeline
 *
 * Y-chromosome haplogroup inference for genomic samples.
 *
 * Author: yallHap Contributors
 * License: MIT
 */

nextflow.enable.dsl = 2

// Import modules
include { YALLHAP_CLASSIFY } from './modules/yallhap'

// Default parameters
params.input        = null
params.tree         = null
params.snp_db       = null
params.reference    = 'grch38'
params.ancient      = false
params.outdir       = './results'
params.help         = false

// Help message
def helpMessage() {
    log.info """
    =========================================
     yallHap Nextflow Pipeline v${workflow.manifest.version}
    =========================================

    Usage:
      nextflow run main.nf --input samples.csv --tree yfull_tree.json --snp_db ybrowse_snps.csv

    Required arguments:
      --input       Path to sample sheet CSV (sample_id,vcf_path)
      --tree        Path to YFull tree JSON file
      --snp_db      Path to SNP database CSV file

    Optional arguments:
      --reference   Reference genome: grch37, grch38, t2t [default: grch38]
      --ancient     Enable ancient DNA mode [default: false]
      --outdir      Output directory [default: ./results]
      --help        Show this help message

    Example:
      nextflow run main.nf \\
        --input samples.csv \\
        --tree data/yfull_tree.json \\
        --snp_db data/ybrowse_snps.csv \\
        --outdir results/

    Sample sheet format (CSV):
      sample_id,vcf_path
      SAMPLE1,/path/to/sample1.vcf.gz
      SAMPLE2,/path/to/sample2.vcf.gz
    """.stripIndent()
}

// Validate inputs
def validateParams() {
    if (params.help) {
        helpMessage()
        System.exit(0)
    }

    if (!params.input) {
        log.error "ERROR: --input is required"
        helpMessage()
        System.exit(1)
    }

    if (!params.tree) {
        log.error "ERROR: --tree is required"
        helpMessage()
        System.exit(1)
    }

    if (!params.snp_db) {
        log.error "ERROR: --snp_db is required"
        helpMessage()
        System.exit(1)
    }
}

// Main workflow
workflow {
    validateParams()

    // Read sample sheet
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            tuple(row.sample_id, file(row.vcf_path))
        }
        .set { samples_ch }

    // Reference files
    tree_file = file(params.tree)
    snp_db_file = file(params.snp_db)

    // Run yallHap classification
    YALLHAP_CLASSIFY(
        samples_ch,
        tree_file,
        snp_db_file,
        params.reference,
        params.ancient
    )

    // Collect all results into a single TSV
    YALLHAP_CLASSIFY.out.tsv
        .collectFile(
            name: 'yallhap_results.tsv',
            storeDir: params.outdir,
            keepHeader: true,
            skip: 1
        )
}

// Workflow completion
workflow.onComplete {
    log.info """
    =========================================
    Pipeline completed at: ${workflow.complete}
    Duration: ${workflow.duration}
    Success: ${workflow.success}
    Work directory: ${workflow.workDir}
    =========================================
    """.stripIndent()
}

