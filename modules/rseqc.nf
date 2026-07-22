/************************************************************************
* RSeQC
************************************************************************/
process rseqc_bam_stat {
    label 'rseqc'
    tag "$meta.sample"

    conda "${projectDir}/envs/rseqc.yaml"

    if ( params.softlink_results ) {
        publishDir "${params.output}/${params.rseqc_dir}/bam_stat", pattern: "*.txt"
    } else {
        publishDir "${params.output}/${params.rseqc_dir}/bam_stat", mode: 'copy', pattern: "*.txt"
    }

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.sample}.bam_stat.txt"), emit: bam_stat

    script:
    """
    bam_stat.py -i ${bam} > ${meta.sample}.bam_stat.txt
    """
}

process rseqc_read_duplication {
    label 'rseqc'
    tag "$meta.sample"

    if ( params.softlink_results ) {
        publishDir "${params.output}/${params.rseqc_dir}/read_duplication", pattern: "*.DupRate*"
    } else {
        publishDir "${params.output}/${params.rseqc_dir}/read_duplication", mode: 'copy', pattern: "*.DupRate*"
    }

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.DupRate*"), emit: read_duplication_files

    script:
    """
    read_duplication.py -i ${bam} -o ${meta.sample}.DupRate
    """
}

process rseqc_gene_body_coverage {
    label 'rseqc'
    tag "$meta.sample"

    if ( params.softlink_results ) {
        publishDir "${params.output}/${params.rseqc_dir}/gene_body_coverage", pattern: "*geneBodyCoverage*"
    } else {
        publishDir "${params.output}/${params.rseqc_dir}/gene_body_coverage", mode: 'copy', pattern: "*geneBodyCoverage*"
    }

    input:
    tuple val(meta), path(bam), path(bai)
    path bed

    output:
    tuple val(meta), path("*geneBodyCoverage*"), emit: gene_body_coverage_files

    script:
    """
    geneBody_coverage.py \
        -i ${bam} \
        -r ${bed} \
        -o ${meta.sample}.geneBodyCoverage
    """
}

process rseqc_infer_experiment {
    label 'rseqc'
    tag "$meta.sample"

    publishDir "${params.output}/${params.rseqc_dir}/infer_experiment",
        mode: params.softlink_results ? 'symlink' : 'copy',
        pattern: "*.txt"

    input:
    tuple val(meta), path(bam), path(bai)
    path bed

    output:
    tuple val(meta), path("${meta.sample}.infer_experiment.txt"), emit: infer_experiment

    script:
    """
    infer_experiment.py -i ${bam} -r ${bed} > ${meta.sample}.infer_experiment.txt
    """
}

process rseqc_read_distribution {
    label 'rseqc'
    tag "$meta.sample"

    publishDir "${params.output}/${params.rseqc_dir}/read_distribution",
        mode: params.softlink_results ? 'symlink' : 'copy',
        pattern: "*.txt"

    input:
    tuple val(meta), path(bam), path(bai)
    path bed

    output:
    tuple val(meta), path("${meta.sample}.read_distribution.txt"), emit: read_distribution

    script:
    """
    read_distribution.py -i ${bam} -r ${bed} > ${meta.sample}.read_distribution.txt
    """
}

process rseqc_junction_annotation {
    label 'rseqc'
    tag "$meta.sample"

    publishDir "${params.output}/${params.rseqc_dir}/junction_annotation",
        mode: params.softlink_results ? 'symlink' : 'copy',
        pattern: "*junction*"

    input:
    tuple val(meta), path(bam), path(bai)
    path bed

    output:
    tuple val(meta), path("*junctionAnnotation*"), emit: junction_annotation_files

    script:
    """
    junction_annotation.py \
        -i ${bam} \
        -r ${bed} \
        -o ${meta.sample}.junctionAnnotation \
        2> ${meta.sample}.junctionAnnotation.log
    """
}

process rseqc_junction_saturation {
    label 'rseqc'
    tag "$meta.sample"

    publishDir "${params.output}/${params.rseqc_dir}/junction_saturation",
        mode: params.softlink_results ? 'symlink' : 'copy',
        pattern: "*junctionSaturation*"

    input:
    tuple val(meta), path(bam), path(bai)
    path bed

    output:
    tuple val(meta), path("*junctionSaturation*"), emit: junction_saturation_files

    script:
    """
    junction_saturation.py \
        -i ${bam} \
        -r ${bed} \
        -o ${meta.sample}.junctionSaturation \
        2> ${meta.sample}.junctionSaturation.log
    """
}

process rseqc_inner_distance {
    label 'rseqc'
    tag "$meta.sample"

    publishDir "${params.output}/${params.rseqc_dir}/inner_distance",
        mode: params.softlink_results ? 'symlink' : 'copy',
        pattern: "*inner_distance*"

    input:
    tuple val(meta), path(bam), path(bai)
    path bed

    output:
    tuple val(meta), path("*inner_distance*"), emit: inner_distance_files

    script:
    """
    inner_distance.py \
        -i ${bam} \
        -r ${bed} \
        -o ${meta.sample}.inner_distance
    """
}

