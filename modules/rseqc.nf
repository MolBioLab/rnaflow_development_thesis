/************************************************************************
* RSeQC
************************************************************************/
process rseqc_gtf_to_bed {
    label 'rseqc'
    tag "annotation"

    if ( params.softlink_results ) {
        publishDir "${params.output}/${params.rseqc_dir}/reference", pattern: "*.bed"
    } else {
        publishDir "${params.output}/${params.rseqc_dir}/reference", mode: 'copy', pattern: "*.bed"
    }

    input:
    path gtf

    output:
    path("annotation.bed"), emit: bed

    script:
    """
    python - <<'PY'
    import re
    from collections import defaultdict

    gtf_file = "${gtf}"
    out_file = "annotation.bed"

    def parse_attr(attr):
        d = {}
        for m in re.finditer(r'(\\S+) "([^"]+)"', attr):
            d[m.group(1)] = m.group(2)
        return d

    transcripts = defaultdict(lambda: {
        "chrom": None,
        "strand": None,
        "exons": []
    })

    with open(gtf_file) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\\n").split("\\t")
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attr = fields
            if feature != "exon":
                continue

            info = parse_attr(attr)
            transcript_id = info.get("transcript_id")
            if transcript_id is None:
                continue

            start = int(start) - 1
            end = int(end)

            rec = transcripts[transcript_id]
            rec["chrom"] = chrom
            rec["strand"] = strand
            rec["exons"].append((start, end))

    with open(out_file, "w") as out:
        for transcript_id, rec in transcripts.items():
            exons = sorted(rec["exons"])
            if not exons:
                continue

            tx_start = exons[0][0]
            tx_end = exons[-1][1]
            chrom = rec["chrom"]
            strand = rec["strand"]

            block_count = len(exons)
            block_sizes = ",".join(str(e - s) for s, e in exons) + ","
            block_starts = ",".join(str(s - tx_start) for s, e in exons) + ","

            out.write(
                f"{chrom}\\t{tx_start}\\t{tx_end}\\t{transcript_id}\\t0\\t{strand}\\t"
                f"{tx_start}\\t{tx_end}\\t0\\t{block_count}\\t{block_sizes}\\t{block_starts}\\n"
            )
    PY
    """
}

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

