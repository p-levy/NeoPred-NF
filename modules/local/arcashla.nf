// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCAS_HLA {
    tag "$meta.id"
    label 'process_high'
    // publishDir "${params.outdir}",
    //     mode: params.publish_dir_mode,
    //     saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta.patient), path("${prefix}/*.alleles.tsv"), emit: hla
    tuple val(meta.patient), path("${prefix}/*.loh.tsv"),    emit: loh
    path  "*.version.txt"  ,                                 emit: version

    script:
    def software = getSoftwareName(task.process)
    def VERSION = 0.5.0
    def base = "${bam.baseName}"
    prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    mkdir ${meta.id}

    arcasHLA extract -o ${meta.id}/extract -t ${task.cpus} $bam

    arcasHLA genotype ${meta.id}/${base}.extracted.1.fq.gz ${meta.id}/${base}.extracted.2.fq.gz -t ${task.cpus} -o ${meta.id}/genotype

    arcasHLA customize -G ${meta.id}/${meta.id}.genotype.json -o ${meta.id}/ref -t ${task.cpus}

    arcasHLA quant --ref ${meta.id}/${meta.id} -t ${task.cpus} -o ${meta.id}/quant ${meta.id}/${base}.extracted.1.fq.gz ${meta.id}/${base}.extracted.2.fq.gz

    cat $VERSION > ${software}.version.txt

    """
}
