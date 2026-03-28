/************************************************************************
* STAR INDEX (genomeGenerate)
************************************************************************/

process starindex {

    label 'star'

    input:
    path(reference)
    path(annotation)

    output:
    tuple path(reference), path("star_index")

    script:
    """
    rm -rf star_index

    STAR \
        --runMode genomeGenerate \
        --runThreadN ${task.cpus} \
        --genomeDir star_index \
        --genomeFastaFiles ${reference} \
        --sjdbGTFfile ${annotation} \
        --sjdbOverhang ${params.star_sjdbOverhang ?: 100}
    """
}


/************************************************************************
* STAR ALIGN
************************************************************************/

process star {

    label 'star'
    tag "$meta.sample"

    if ( params.softlink_results ) {
        publishDir "${params.output}/${params.star_dir}", pattern: "*.sorted.bam"
    } else {
        publishDir "${params.output}/${params.star_dir}", mode: 'copy', pattern: "*.sorted.bam"
    }

    input:
    tuple val(meta), path(reads)
    tuple path(reference), path(star_index)
    path(annotation)
    val(additionalParams)

    output:
    tuple val(meta), path("${meta.sample}.sorted.bam"), emit: sample_bam
    path "${meta.sample}.Log.final.out", emit: log

    script:

    // detect gzipped fastq
    def gz = reads[0].name.endsWith('.gz')
    def read_cmd = gz ? '--readFilesCommand zcat' : ''

    // only use tmp directory when not using conda/mamba
    def tmp_opt = (workflow.profile.contains('conda') || workflow.profile.contains('mamba')) \
        ? '' : "--outTmpDir tmp-star-${meta.sample}"

    if ( !meta.paired_end ) {

        """
        rm -rf tmp-star-${meta.sample}

        STAR \
            --genomeDir ${star_index} \
            --readFilesIn ${reads[0]} \
            ${read_cmd} \
            --sjdbGTFfile ${annotation} \
            --runThreadN ${task.cpus} \
            --outFileNamePrefix ${meta.sample}. \
            --outSAMtype BAM SortedByCoordinate \
            ${tmp_opt} \
            ${additionalParams}

        mv ${meta.sample}.Aligned.sortedByCoord.out.bam ${meta.sample}.sorted.bam

        rm -rf tmp-star-${meta.sample}
        """

    } else {

        """
        rm -rf tmp-star-${meta.sample}

        STAR \
            --genomeDir ${star_index} \
            --readFilesIn ${reads[0]} ${reads[1]} \
            ${read_cmd} \
            --sjdbGTFfile ${annotation} \
            --runThreadN ${task.cpus} \
            --outFileNamePrefix ${meta.sample}. \
            --outSAMtype BAM SortedByCoordinate \
            ${tmp_opt} \
            ${additionalParams}

        mv ${meta.sample}.Aligned.sortedByCoord.out.bam ${meta.sample}.sorted.bam

        rm -rf tmp-star-${meta.sample}
        """
    }
}

