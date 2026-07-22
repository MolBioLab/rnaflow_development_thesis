#!/usr/bin/env python3

import argparse
import csv
import os
import re
from collections import defaultdict


def parse_attrs(attr_string):
    return dict(re.findall(r'(\S+)\s+"([^"]+)"', attr_string))


def parse_gtf_gene_biotypes(gtf_file, biotype_attr):
    gene_to_biotype = {}

    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            attrs = parse_attrs(fields[8])
            gene_id = attrs.get("gene_id")
            biotype = attrs.get(biotype_attr)

            if gene_id and biotype:
                gene_to_biotype[gene_id] = biotype

    return gene_to_biotype


def sample_name_from_counts_file(path):
    name = os.path.basename(path)

    for suffix in [".counts.tsv", ".featureCounts.tsv", ".tsv"]:
        if name.endswith(suffix):
            return name[: -len(suffix)]

    return os.path.splitext(name)[0]


def read_featurecounts_one_file(counts_file, gene_to_biotype):
    sample = sample_name_from_counts_file(counts_file)
    biotype_counts = defaultdict(int)

    with open(counts_file) as f:
        reader = csv.reader(f, delimiter="\t")

        header = None
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            header = row
            break

        if header is None or "Geneid" not in header:
            raise RuntimeError(f"Could not find featureCounts header in {counts_file}")

        gene_idx = header.index("Geneid")
        count_idx = len(header) - 1

        for row in reader:
            if not row:
                continue

            gene_id = row[gene_idx]
            count = int(float(row[count_idx]))
            biotype = gene_to_biotype.get(gene_id, "Unknown")
            biotype_counts[biotype] += count

    return sample, biotype_counts


def write_multiqc_custom_content(out_file, sample_to_counts):
    all_biotypes = sorted(
        {bt for counts in sample_to_counts.values() for bt in counts.keys()}
    )

    with open(out_file, "w") as out:
        out.write("# plot_type: bargraph\n")
        out.write("# id: featurecounts_biotype\n")
        out.write("# section_name: featureCounts: Biotypes\n")
        out.write("# description: Read counts summed by gene biotype from featureCounts gene-level counts.\n")
        out.write("# pconfig:\n")
        out.write("#     id: featurecounts_biotype_plot\n")
        out.write("#     title: featureCounts: Biotypes\n")
        out.write("#     ylab: Read counts\n")

        out.write("Sample\t" + "\t".join(all_biotypes) + "\n")

        for sample in sorted(sample_to_counts):
            counts = sample_to_counts[sample]
            values = [str(counts.get(bt, 0)) for bt in all_biotypes]
            out.write(sample + "\t" + "\t".join(values) + "\n")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--counts",
        required=True,
        nargs="+",
        help="One or more featureCounts gene-level TSV files"
    )
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--biotype", default="gene_biotype")

    args = parser.parse_args()

    gene_to_biotype = parse_gtf_gene_biotypes(args.gtf, args.biotype)

    sample_to_counts = {}

    for counts_file in args.counts:
        sample, biotype_counts = read_featurecounts_one_file(
            counts_file,
            gene_to_biotype
        )
        sample_to_counts[sample] = biotype_counts

    write_multiqc_custom_content(args.out, sample_to_counts)


if __name__ == "__main__":
    main()
