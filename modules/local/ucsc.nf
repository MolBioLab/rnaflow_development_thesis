process ucsc_gtf_to_genepred {
    label 'ucsc_gtftogenepred'
    tag "annotation"

    input:
    path gtf

    output:
    path "annotation.genePred", emit: genepred

    script:
    """
    gtfToGenePred \
        ${gtf} \
        annotation.genePred
    """
}

process ucsc_genepred_to_bed {
    label 'ucsc_genepredtobed'
    tag "annotation"

    publishDir "${params.output}/${params.rseqc_dir}/reference",
        mode: params.softlink_results ? 'symlink' : 'copy',
        pattern: "*.bed"

    input:
    path genepred

    output:
    path "annotation.bed", emit: bed

    script:
    """
    genePredToBed \
        ${genepred} \
        annotation.bed
    """
}

workflow ucsc_gtf_to_bed {

    take:
    gtf

    main:
    ucsc_gtf_to_genepred(gtf)
    ucsc_genepred_to_bed(ucsc_gtf_to_genepred.out.genepred)

    emit:
    bed = ucsc_genepred_to_bed.out.bed
}
