// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process STRELKA_FILTER_INDEL {
    tag "Strelka_indel_filtering"
    label 'process_low'
    // publishDir "${params.outdir}",
    //     mode: params.publish_dir_mode,
    //     saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['patient']) }

    input:
    tuple val(meta), path(strelka_indel_vcf)

    output:
    tuple val(meta), path("*.vcf")      , emit: vcf

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${strelka_indel_vcf.getBaseName(times=2)}${options.suffix}" : "${strelka_indel_vcf.getBaseName(times=2)}"

    """
    strelka_indel_filter.py ${strelka_indel_vcf} ${prefix}.vcf
    """
}
