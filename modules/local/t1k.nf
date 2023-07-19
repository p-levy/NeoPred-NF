// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process T1K {
    tag "$meta.id"
    label 'process_medium'
    // publishDir "${params.outdir}",
    //     mode: params.publish_dir_mode,
    //     saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta.patient), path("*_allele.tsv"), emit: hla_allele
    tuple val(meta.patient), path("*_genotype.tsv"), emit: hla_genotype
    path  "*.version.txt"  , emit: version

    script:
    def software = getSoftwareName(task.process)
    def VERSION = "1.0.2"
    prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    run-t1k -f /t1k-ref/hlaidx_dna_seq.fa -1 ${reads[0]} -2 ${reads[1]} -t ${task.cpus} --preset hla -o ${prefix}

    echo $VERSION > ${software}.version.txt

    """
}
