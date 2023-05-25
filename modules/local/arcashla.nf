// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARCASHLA {
    tag "$meta.id"
    label 'process_high'
    // publishDir "${params.outdir}",
    //     mode: params.publish_dir_mode,
    //     saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta.patient), path("*.alleles.tsv"), emit: hla
    tuple val(meta.patient), path("*.loh.tsv"),    emit: loh
    path  "*.version.txt"  ,                                 emit: version

    script:
    def software = getSoftwareName(task.process)
    def VERSION = "0.5.0"
    def base = "${bam.baseName}"
    prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    mkdir ${meta.id}

    arcasHLA extract -o . -t ${task.cpus} $bam

    arcasHLA genotype ${base}.extracted.1.fq.gz ${base}.extracted.2.fq.gz -t ${task.cpus} -o .

    arcasHLA customize -G ${meta.id}.genotype.json -o . -t ${task.cpus}

    arcasHLA quant --LOH --ref ${meta.id} -t ${task.cpus} -o . ${base}.extracted.1.fq.gz ${base}.extracted.2.fq.gz

    echo $VERSION > ${software}.version.txt

    """
}
