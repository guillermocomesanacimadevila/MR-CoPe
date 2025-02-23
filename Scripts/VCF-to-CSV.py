# ==== Call libraries ==== #

import os
import csv
import pandas as pd

# ==== Function accumulator ==== #

def vcf_to_csv(vcf_file_path, csv_file_path):
    with open(vcf_file_path, "r") as vcf_file:
        with open(csv_file_path, "w", newline="") as csv_file:
            csv_writer = csv.writer(csv_file)
            for line in vcf_file:
                if line.startswith("##"):
                    continue
                elif line.startswith("#"):
                    header = line.strip("#").strip().split("\t")
                    csv_writer.writerow(header)
                else:
                    data = line.strip().split("\t")
                    csv_writer.writerow(data)

            return data

# Exposure VCF and CSVº paths
exposure_csv_path = "~/cpep_MR/Data/exposure.vcf"
exposure_csv_path = "~/cpep_MR/Data/exposure.csv"

# Outcome VCF AND CSV paths
outcome_vcf_path = "~/cpep_MR/Data/outcome.vcf"
outcome_csv_path = "~/cpep_MR/Data/outcome.csv"

vcf_to_csv(os.path.expanduser(exposure_csv_path), os.path.expanduser(exposure_csv_path))
vcf_to_csv(os.path.expanduser(outcome_vcf_path), os.path.expanduser(outcome_csv_path))