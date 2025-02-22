#!bin/bash

# Download exposure GWAS data
wget https://gwas.mrcieu.ac.uk/files/ieu-b-5067/ieu-b-5067.vcf.gz >> exposure.vcf.gz

# Download outcome GWAS data
wget https://gwas.mrcieu.ac.uk/files/ieu-b-110/ieu-b-110.vcf.gz >> outcome.vcf.gz

echo "VCFs downloaded!"

# gzunip VCFs
gunzip exposure.vcf.gz
gunzip outcome.vcf.gz

echo "Process completed!"
