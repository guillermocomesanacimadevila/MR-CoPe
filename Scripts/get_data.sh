#!bin/bash

# Download exposure GWAS data
curl -O https://gwas.mrcieu.ac.uk/files/ieu-b-5067/ieu-b-5067.vcf.gz 

# Download outcome GWAS data
curl -O https://gwas.mrcieu.ac.uk/files/ieu-b-110/ieu-b-110.vcf.gz 

echo "VCFs downloaded!"

# gzunip VCFs
gunzip ieu-b-110.vcf.gz
gunzip ieu-b-5067.vcf.gz 

echo "Unzipping done!"

# Rename VCFs for clarity
mv ieu-b-110.vcf exposure.vcf
mv ieu-b-5067.vcf outcome.vcf

echo "Process completed!"
