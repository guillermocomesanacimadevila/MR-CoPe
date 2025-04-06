## MR-CoPe Description 🧬🧬

Guillermo Comesaña - MSc Bioinformatics @University of Bath

Christian Pepler - MRes Bioscience @Cardiff University

## Tools & Languages
<p align="left">
  <img src="https://github.com/user-attachments/assets/5e678fc0-9597-4252-98dd-eb9aaccc823e" alt="Python" width="60" style="margin: 0 10px;"/>
  <img src="https://github.com/user-attachments/assets/a49b35ad-c2f7-4cbe-b755-47ebe3330866" alt="R" width="72" style="margin: 0 10px; position: relative; top: -2px;"/>
  <img src="https://github.com/user-attachments/assets/4bbcf45e-d572-45e9-a16c-3ff379e72390" alt="Bash" width="65" style="margin: 0 10px;"/>
  <img src="https://github.com/user-attachments/assets/805532d9-fc8b-446f-aac6-933cc4aa6185" alt="Git" width="65" style="margin: 0 10px;"/>
  <img src="https://github.com/user-attachments/assets/bfc30e37-cb64-4d59-8cec-52ab5c12fab7" alt="Docker" width="65" style="margin: 0 10px;"/>
  <img src="https://github.com/user-attachments/assets/0427f54d-9e05-4969-91d1-13af16c3fb42" alt="SQL" width="100" style="margin: 0 10px;"/>
</p>

Bioinformatics 💻 / Data Science 📈 / Genetic Epidemiology 🧬

## What is Mendelian Randomisation

![1-s2 0-S0753332222003419-gr1](https://github.com/user-attachments/assets/b51c516e-c858-4d13-8529-8683abdf1e09)

## Run Pipeline -> (UNIX-based)
To run the pipeline, give the script execute permission and run it:

```bash
chmod +x run_pipeline.sh && ./run_pipeline.sh
```
- ** System Requirements ** 
** UNIX-based terminal ** 
** Python3 > 3.10 **
** R > 4.31 ** 
** Java > 17.00 **
** Nextflow > 24.11.0 ** 

## Pipeline Workflow

### Data Validation (Synthetic GWAS)

### Data Validation (Real GWAS)

#### MR(1)
SNPs -> LDl-c -> Alzheimer´s Disease Risk
SNPs -> Alzheimer´s Disease Risk -> LDL-c

#### MR(2)
SNPs -> Iodine-c -> Alzheimer´s Disease Risk
SNPs -> Alzheimer´s Disease Risk -> Iodine-c


### Mendelian Randomisation -> Statistical Methods

#### Inverse Variance Weighted (IVW)

β_IVW = ( ∑ wᵢ · β_Y,ᵢ · β_X,ᵢ ) / ( ∑ wᵢ · β_X,ᵢ² )

Where:
- β_Y,ᵢ = SNP–outcome association
- β_X,ᵢ = SNP–exposure association
- wᵢ = 1 / SE(β_Y,ᵢ)² = inverse variance of outcome association

#### Weighted Median Estimate (WME)

WME = median( β_Y,ᵢ / β_X,ᵢ ), weighted by 1 / SE(β_Y,ᵢ)²

Where:
- β_Y,ᵢ = SNP–outcome association
- β_X,ᵢ = SNP–exposure association
- SE(β_Y,ᵢ) = standard error of the SNP–outcome association

#### MR-Egger

β_Y,ᵢ = α + β_MR · β_X,ᵢ + εᵢ

Where:
- β_Y,ᵢ = SNP–outcome association
- β_X,ᵢ = SNP–exposure association
- α = intercept term (captures directional pleiotropy)
- β_MR = estimated causal effect (the slope)
- εᵢ = error term
