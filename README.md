# Extended RNAflow: A Nextflow Pipeline for Bulk RNA-seq Differential Gene Expression Analysis

An extended version of [RNAflow](https://github.com/hoelzer-science/rnaflow), developed to provide **alternative splice-aware read alignment with STAR or HISAT2** and some **additional post-alignment RNA-seq quality assessment with RSeQC** (not full). The workflow retains the core RNAflow approach to preprocessing, gene-level quantification, differential expression analysis, and reporting.

> [!IMPORTANT]
> This extended pipeline was developed and run using **Docker**. **Conda and Singularity/Apptainer are not yet supported or validated for the newly developed STAR and RSeQC functionality.** Please use the Docker profile for workflows that enable these extensions. This limitation concerns the extended implementation, not necessarily the upstream RNAflow project.

## Origin and attribution

This repository is a derivative of **RNAflow** by the Hoelzer laboratory. For the original project's features, documentation, citation instructions, and history, please see **[the original RNAflow repository](https://github.com/hoelzer-science/rnaflow)**.

- Upstream RNAflow baseline used for development: `f2d26a3418adc10d24b5042ad653081a9500245b`.

This README describes the extended differential-expression workflow. Other features inherited from upstream RNAflow, including transcriptome assembly, are not re-documented here; refer to the upstream project.

## Pipeline overview

![Overview of the extended RNAflow pipeline](figures/developed_pipeline_overview.jpg)

*Fig 1. Pipeline overview shows the STAR/HISAT2 alignment options and RSeQC's modules.*

The analysis proceeds through:

1. **Preprocessing:** FastQC on raw reads → fastp trimming/filtering → FastQC on processed reads; optional SortMeRNA rRNA removal.
2. **Alignment:** choose HISAT2 (original/default) or STAR (added), using a genome FASTA and gene annotation; produce sorted BAM and index files.
3. **Post-alignment QC (optional):** invoke RSeQC analyses on alignment files, using suitable reference/annotation inputs where required.
4. **Quantification and expression:** featureCounts → optional TPM-based filtering → DESeq2 differential expression analysis.
5. **Downstream analysis and reporting:** inherited piano/WebGestalt pathway analysis options and MultiQC summaries, where supported by the chosen workflow.

## What's new in this extension?

### Optional STAR alignment

The original HISAT2 alignment branch is retained, while `--aligner star` selects an alternative STAR alignment branch. The genome index is generated for the selected aligner, and aligned data are passed to the existing downstream quantification and DEG analysis stages.

| Parameter | Meaning | Default / behavior |
| --- | --- | --- |
| `--aligner hisat2` | Use the original HISAT2 alignment branch | **Default** if omitted |
| `--aligner star` | Use the added STAR alignment branch | Opt-in |

### Optional RSeQC modules (original paper: https://academic.oup.com/bioinformatics/article/28/16/2184/325191)

Some RSeQC analyses (https://rseqc.sourceforge.net/) are controlled by the pipeline-level flag `--run_rseqc true`. When the flag is `false` or omitted, the added RSeQC workflow is skipped.

| RSeQC module / tool | Purpose | Key inputs and applicability |
| --- | --- | --- |
| `bam_stat.py` | Summary of alignment composition and mapping characteristics | Alignment BAM |
| `read_duplication.py` | Estimates read duplication characteristics | Alignment BAM |
| `geneBody_coverage.py` | Normalized 5′–3′ gene-body coverage; detects potential positional bias | Alignment BAM and transcript/gene-model BED |
| `infer_experiment.py` | Infers library strandedness/orientation | Alignment BAM and reference gene-model BED |
| `read_distribution.py` | Summarizes read assignment across annotated genomic features | Alignment BAM and compatible BED annotation |
| `junction_saturation.py` | Assesses how splice-junction discovery changes with read subsampling | Spliced alignment BAM and compatible gene-model annotation; dependent on installed RSeQC version |
| `inner_distance.py` | Estimates the inner distance between paired reads | **Paired-end** alignment BAM and compatible BED annotation; not applicable to single-end libraries |

## Requirements

- Linux environment with a working **Docker daemon** and permission to run containers.
- **Nextflow** and a compatible Java runtime installed on the host.
- Sufficient disk space for FASTQ files, indices, intermediate BAMs, Nextflow `work/`, and results.
- Reference genome FASTA and matching gene annotation GTF/GFF for the organism being analyzed, unless using a supported upstream RNAflow download RNAflow option.

## Quick start (Docker)

Run these commands **from the root directory of the cloned extended repository**, not from the upstream RNAflow repository:

```bash
# Verify the execution environment
nextflow -version
docker --version
docker info

# Inspect the available arguments in this checkout
nextflow run . --help

# Run the extended workflow on prepared input files
nextflow run . -profile docker,local \
  --reads samples.csv \
  --genome genome_files.csv \
  --aligner star \
  --run_rseqc true \
  --max_cores 6 \
  --cores 6 \
  -resume
```
To use the HISAT2 branch, substitute `--aligner hisat2` or omit `--aligner`. To skip RSeQC, use `--run_rseqc false` or omit the flag.

### Input samplesheet

RNAflow's original `--reads` format is CSV with this exact header:

```csv
Sample,R1,R2,Condition,Source,Strandedness
control_1,/path/control_1_R1.fastq.gz,/path/control_1_R2.fastq.gz,control,,0
control_2,/path/control_2_R1.fastq.gz,/path/control_2_R2.fastq.gz,control,,0
treated_1,/path/treated_1_R1.fastq.gz,/path/treated_1_R2.fastq.gz,treated,,0
treated_2,/path/treated_2_R1.fastq.gz,/path/treated_2_R2.fastq.gz,treated,,0
```

- `Sample`: unique sample label.
- `R1`, `R2`: FASTQ paths. For single-end samples, leave `R2` blank.
- `Condition`: experimental group for DEG comparisons. The upstream workflow expects at least two replicates per condition.
- `Source`: optional sample/patient or paired-design label; leave blank if not applicable.
- `Strandedness`: `0` = unstranded, `1` = forward stranded, `2` = reverse stranded, according to upstream RNAflow's convention. Set this to the actual library type.

### Optional configuration examples

```bash
# Default/original alignment, with the added RSeQC steps disabled
nextflow run . -profile docker,local \
  --reads samples.csv --genome genome_files.csv \
  --aligner hisat2 --run_rseqc false

# STAR alignment with the optional RSeQC workflow enabled
nextflow run . -profile docker,local \
  --reads samples.csv --genome genome_files.csv \
  --aligner star --run_rseqc true

# Avoid repeating completed tasks when resuming the same run
nextflow run . -profile docker,local \
  --reads samples.csv --genome genome_files.csv \
  --aligner star --run_rseqc true -resume
```

Other inherited upstream RNAflow controls include `--skip_sortmerna`, `--skip_read_preprocessing`, `--deg`, `--tpm`, and `--featurecounts_additional_params`. Refer to the [original RNAflow documentation](https://github.com/hoelzer-science/rnaflow) for their full usage.

## Expected outputs and quality-control guidance

Output paths depend on the current workflow implementation and configured results directory. Look for the following **types** of results rather than treating the illustrative labels below as guaranteed exact folder names:

| Stage | Outputs / checks |
| --- | --- |
| FastQC / fastp | Raw and trimmed read quality, adapters, trimming statistics |
| SortMeRNA, if enabled | rRNA filtering logs / processed reads |
| HISAT2 or STAR | Alignment log/summary, coordinate-sorted BAM and BAI |
| RSeQC `bam_stat` | Alignment summary statistics |
| RSeQC `read_duplication` | Duplicate-read metrics and associated plots/data |
| RSeQC `geneBody_coverage` | Gene-body coverage tables and/or plots from 5′ to 3′ |
| RSeQC `infer_experiment` | Estimates of read orientation relative to annotated genes |
| RSeQC `read_distribution` | Distribution across exonic, intronic, and other annotated regions |
| RSeQC `junction_saturation` | Splice-junction saturation output/plots |
| RSeQC `inner_distance` | Paired-read insert/inner-distance distribution |
| featureCounts | Gene-level count matrix and assignment summary |
| DESeq2 | Differential gene-expression tables and diagnostic plots |
| MultiQC | Combined QC report for supported tools and formats |

## Troubleshooting

**Docker process fails or cannot access files.** Confirm Docker is running (`docker info`) and that all FASTQ, genome, and annotation paths are readable from the host and accessible to Docker. Ensure sufficient storage for `work/` and image pulls.

**STAR index fails or the process is killed.** Examine RAM use and Nextflow process logs. Whole-genome STAR indices can demand substantially more memory than HISAT2 indices; use appropriately sized hardware and check process resource limits.

**Conda/Singularity execution fails.** For extended STAR/RSeQC operations, switch to the tested `-profile docker,local` configuration. The extended modules have not yet been adapted and validated for these alternate environments.

## Citation and license

This project builds on the original [RNAflow repository](https://github.com/hoelzer-science/rnaflow). Please cite RNAflow and the relevant third-party tools according to the upstream project's [citation guidance](https://github.com/hoelzer-science/rnaflow/blob/master/citations.md).

The upstream RNAflow repository is licensed under **GPL-3.0**. Retain the upstream license file and applicable attribution/notices in the derivative repository; see [`LICENSE`](LICENSE) for the repository's applicable license text.
