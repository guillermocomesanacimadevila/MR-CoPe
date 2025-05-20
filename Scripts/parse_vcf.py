#!/usr/bin/env python3

import gzip
import csv
import sys

vcf_path = sys.argv[1]
csv_path = sys.argv[2]
is_log10 = sys.argv[3].lower() == "true"

def parse_sample_field(field):
    es, se, lp, af, rsid = field.split(":")
    try:
        pval = 10 ** (-float(lp)) if is_log10 else float(lp)
    except:
        pval = None
    return {
        "BETA": float(es),
        "SE": float(se),
        "PVALUE": pval,
        "EAF": float(af),
        "SNP": rsid
    }

with gzip.open(vcf_path, 'rt') if vcf_path.endswith('.gz') else open(vcf_path, 'r') as vcf_in, \
     open(csv_path, 'w', newline='') as csv_out:

    writer = None
    for line in vcf_in:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            header = line.strip().lstrip("#").split("\t")
            format_col = header.index("FORMAT")
            sample_col = format_col + 1
            writer = csv.DictWriter(csv_out, fieldnames=[
                "SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "PVALUE", "EAF"
            ])
            writer.writeheader()
            continue

        fields = line.strip().split("\t")
        if len(fields) <= sample_col:
            continue
        chrom, pos, snp_id, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
        try:
            stats = parse_sample_field(fields[sample_col])
            writer.writerow({
                "SNP": stats["SNP"],
                "CHR": chrom,
                "BP": pos,
                "A1": alt,
                "A2": ref,
                "BETA": stats["BETA"],
                "SE": stats["SE"],
                "PVALUE": stats["PVALUE"],
                "EAF": stats["EAF"]
            })
        except Exception as e:
            print(f"⚠️ Skipping line at position {pos}: {e}")

