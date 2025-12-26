/*
 * yallHap Nextflow Module
 *
 * Process for Y-chromosome haplogroup classification using yallHap.
 */

process YALLHAP_CLASSIFY {
    tag "$sample_id"

    label 'process_low'

    container 'yallhap/yallhap:latest'

    publishDir "${params.outdir}/individual", mode: 'copy', pattern: '*.json'

    input:
    tuple val(sample_id), path(vcf)
    path tree
    path snp_db
    val reference
    val ancient

    output:
    tuple val(sample_id), path("${sample_id}.json"), emit: json
    tuple val(sample_id), path("${sample_id}.tsv"),  emit: tsv

    script:
    def ancient_flag = ancient ? '--ancient' : ''
    """
    yallhap classify \\
        ${vcf} \\
        --tree ${tree} \\
        --snp-db ${snp_db} \\
        --reference ${reference} \\
        --sample ${sample_id} \\
        ${ancient_flag} \\
        --format json \\
        --output ${sample_id}.json

    yallhap classify \\
        ${vcf} \\
        --tree ${tree} \\
        --snp-db ${snp_db} \\
        --reference ${reference} \\
        --sample ${sample_id} \\
        ${ancient_flag} \\
        --format tsv \\
        --output ${sample_id}.tsv
    """

    stub:
    """
    echo '{"sample": "${sample_id}", "haplogroup": "R-M343", "confidence": 0.95}' > ${sample_id}.json
    echo -e "sample\thaplogroup\tconfidence\tqc1\tqc2\tqc3\tqc4\tderived\tancestral\tmissing" > ${sample_id}.tsv
    echo -e "${sample_id}\tR-M343\t0.9500\t0.9500\t0.9500\t0.9500\t0.9500\t100\t50\t10" >> ${sample_id}.tsv
    """
}

/*
 * Process for batch classification of multiple samples from a multi-sample VCF
 */
process YALLHAP_BATCH {
    tag "batch"

    label 'process_medium'

    container 'yallhap/yallhap:latest'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path vcf_files
    path tree
    path snp_db
    val reference
    val ancient

    output:
    path "yallhap_batch_results.tsv", emit: tsv

    script:
    def ancient_flag = ancient ? '--ancient' : ''
    """
    yallhap batch \\
        ${vcf_files} \\
        --tree ${tree} \\
        --snp-db ${snp_db} \\
        --reference ${reference} \\
        ${ancient_flag} \\
        --output yallhap_batch_results.tsv
    """

    stub:
    """
    echo -e "sample\thaplogroup\tconfidence\tqc1\tqc2\tqc3\tqc4\tderived\tancestral\tmissing" > yallhap_batch_results.tsv
    echo -e "SAMPLE1\tR-M343\t0.9500\t0.9500\t0.9500\t0.9500\t0.9500\t100\t50\t10" >> yallhap_batch_results.tsv
    """
}

