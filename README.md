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

#### 📐 Inverse-Variance Weighted (IVW) Estimator

The IVW method estimates the causal effect (β_IVW) by performing a weighted regression of 
SNP-outcome effects (β_Yi) on SNP-exposure effects (β_Xi), **without an intercept**:

    β_Yi = β_IVW · β_Xi + ε_i

To obtain β_IVW, we minimize the weighted sum of squared residuals:

             ⎡  ∑ (w_i · β_Xi · β_Yi)  ⎤
    β_IVW = ⎢ ------------------------ ⎥
             ⎣   ∑ (w_i · β_Xi²)      ⎦

Where the weights are defined as the inverse variance of the outcome effect estimates:

    w_i = 1 / SE_Yi²

#### 🧮 Weighted Median Estimator (WME)

The Weighted Median Estimator provides a consistent causal effect estimate 
even when up to 50% of the instruments are invalid.

Given a set of ratio estimates:

    β_i = β_Yi / β_Xi

Each estimate is weighted by the inverse variance of β_Yi:

    w_i = 1 / SE_Yi²

The WME is the **median** of the β_i values, ordered and weighted by w_i.

This method is robust to violations of the exclusion restriction, assuming that
at least 50% of the total weight comes from valid instruments.

#### 📐 MR-Egger Regression

MR-Egger allows for directional pleiotropy by adding an intercept term to the regression model:

$$
\hat{\beta}_{Yi} = \beta_0 + \beta_{\text{Egger}} \cdot \hat{\beta}_{Xi} + \epsilon_i
$$

To estimate \( \beta_{\text{Egger}} \), we perform a **weighted linear regression** of SNP-outcome effects on SNP-exposure effects **with an intercept**. The causal effect is given by the slope:

$$
\beta_{\text{Egger}} = \frac{\sum_{i} w_i (\hat{\beta}_{Xi} - \bar{\beta}_X)(\hat{\beta}_{Yi} - \bar{\beta}_Y)}{\sum_{i} w_i (\hat{\beta}_{Xi} - \bar{\beta}_X)^2}
$$

Where the weights are defined as:

$$
w_i = \frac{1}{SE_{Yi}^2}
$$

The intercept \( \beta_0 \) captures **average directional pleiotropy**.

### Genetic Instrument Strength

#### F-statistic

### Heterogeneity and Horizontal Pleiotropy

#### Cochran Q Statistic + p-Val -> IVW/WME

#### I2 Statistic -> MR-Egger

#### Egger Intercept 
$$
Z_i = \beta_0 + \beta_1 \cdot SE_i^{-1} + \epsilon_i
$$
