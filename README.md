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
To run the pipeline, clone the repo, and give the script execute permission and run it:

```bash
git clone https://github.com/guillermocomesanacimadevila/MR-CoPe.git
```

```bash
cd MR-CoPe
```

```bash
chmod +x run_mrcope.sh && ./run_mrcope.sh
```

## System Requirements ->
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

#### MR(2)
SNPs -> Iodine-c -> Alzheimer´s Disease Risk

### Mendelian Randomisation -> Statistical Methods

#### Inverse Variance Weighted (IVW) Estimator

The IVW method estimates the causal effect (**β<sub>IVW</sub>**) by performing a weighted regression of SNP-outcome effects (**β̂<sub>Yi</sub>**) on SNP-exposure effects (**β̂<sub>Xi</sub>**), **without an intercept**:

> β̂<sub>Yi</sub> = β<sub>IVW</sub> × β̂<sub>Xi</sub> + ε<sub>i</sub>

To obtain **β<sub>IVW</sub>**, we minimize the weighted sum of squared residuals:

> β<sub>IVW</sub> = ∑ w<sub>i</sub> × β̂<sub>Xi</sub> × β̂<sub>Yi</sub> / ∑ w<sub>i</sub> × β̂<sub>Xi</sub><sup>2</sup>

where weights are defined as the inverse of the variance of the outcome effect estimates:

> w<sub>i</sub> = 1 / SE<sub>Yi</sub><sup>2</sup>

#### Weighted Median Estimate (WME)

#### MR-Egger

### Genetic Instrument Strength

#### F-statistic

### Heterogeneity and Horizontal Pleiotropy

#### Cochran Q Statistic + p-Val -> IVW/WME

#### I2 Statistic -> MR-Egger

#### Egger Intercept 
$$
Z_i = \beta_0 + \beta_1 \cdot SE_i^{-1} + \epsilon_i
$$
