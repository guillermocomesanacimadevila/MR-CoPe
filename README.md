## 🧬 MR-CoPe: Mendelian Randomisation Automated Causal Pipeline Engine

A reproducible and scalable pipeline for Mendelian Randomisation analysis using summary-level GWAS data.

---

### 👥 Authors

- **Guillermo Comesaña Cimadevila**  
  MSc Bioinformatics, *University of Bath*  

- **Christian Pepler**  
  MRes Bioscience, *Cardiff University*

---

### 🛠️ Tools & Languages

<p align="left">
  <img src="https://github.com/user-attachments/assets/5e678fc0-9597-4252-98dd-eb9aaccc823e" alt="Python" width="60" style="margin: 0 10px;"/>
  <img src="https://github.com/user-attachments/assets/a49b35ad-c2f7-4cbe-b755-47ebe3330866" alt="R" width="72" style="margin: 0 10px; position: relative; top: -2px;"/>
  <img src="https://github.com/user-attachments/assets/4bbcf45e-d572-45e9-a16c-3ff379e72390" alt="Bash" width="65" style="margin: 0 10px;"/>
  <img src="https://github.com/user-attachments/assets/805532d9-fc8b-446f-aac6-933cc4aa6185" alt="Git" width="65" style="margin: 0 10px;"/>
  <img src="https://github.com/user-attachments/assets/bfc30e37-cb64-4d59-8cec-52ab5c12fab7" alt="Docker" width="65" style="margin: 0 10px;"/>
  <img src="https://github.com/user-attachments/assets/0427f54d-9e05-4969-91d1-13af16c3fb42" alt="SQL" width="100" style="margin: 0 10px;"/>
</p>

### 🔬 Focus Areas

- **Bioinformatics** 💻  
- **Data Science & Statistical Inference** 📈  
- **Genetic Epidemiology & Causal Inference** 🧬  

---

## 🧬 What is Mendelian Randomisation?

<img src="https://github.com/user-attachments/assets/8f2e8edc-e66d-44e6-86c0-1da112b27dd1" width="600" alt="Mendelian Randomisation Diagram" />

**Mendelian Randomisation (MR)** is a method that uses genetic variants as natural experiments (instrumental variables) to estimate the causal effect of a modifiable exposure (e.g., LDL cholesterol, iodine, fatty acids) on an outcome (e.g., Alzheimer’s disease). It leverages the random allocation of alleles at conception, mimicking the properties of a randomised controlled trial.

---

## 🚀 Run the Pipeline 

Follow the steps below to set up and run the MR-CoPe pipeline.

### 1. Clone the Repository

```bash
git clone https://github.com/guillermocomesanacimadevila/MR-CoPe.git
cd MR-CoPe
```

### 2. Ensure Nextflow is Available

```bash
$($(find / -name nextflow -type f 2>/dev/null | head -n 1))
```

### 3. Run the Pipeline

```bash
chmod +x run_mrcope.sh && ./run_mrcope.sh
```

---

## 🔧 System Requirements

| Tool       | Minimum Version |
|------------|-----------------|
| Python     | ≥ 3.10          |
| R          | ≥ 4.3.1         |
| Java       | ≥ 17.0          |
| Nextflow   | ≥ 24.11.0       |
| Docker     | ≥ 20.10         |
| Terminal   | UNIX-based      |

---

## 🧪 Pipeline Workflow


## 🧪 Data Validation

To assess the robustness and reliability of the causal inference pipeline, validation was performed using both synthetic and real genome-wide association study (GWAS) datasets.

---

### 📊 Synthetic GWAS Validation

Synthetic GWAS datasets were generated to simulate known causal effects. These datasets serve as a testbed for verifying the performance and bias characteristics of the Mendelian Randomisation (MR) estimators under controlled conditions (e.g., varying instrument strength, sample size, and pleiotropy).

---

### 📈 Real GWAS Validation

#### MR(1): Fatty Acids → Alzheimer’s Disease Risk

This MR analysis evaluates the potential **causal effect of genetically predicted circulating fatty acid levels** (e.g., LDL-C, omega-3) on the **risk of Alzheimer’s Disease (AD)**.

This analysis helps validate the pipeline in a real-world, biologically relevant context.

---

## 🔍 Mendelian Randomisation – Statistical Methods

Mendelian Randomisation (MR) uses genetic variants as instrumental variables (IVs) to estimate the causal effect of an exposure on an outcome. Below are the main estimators used in two-sample MR, including assumptions and derivations.

---

### 📐 Inverse-Variance Weighted (IVW) Estimator

The IVW estimator assumes all instruments are valid (i.e., no pleiotropy). It regresses the SNP-outcome associations (β_Yi) on SNP-exposure associations (β_Xi), **without an intercept**:

    β_Yi = β_IVW · β_Xi + ε_i

To estimate β_IVW, we minimize the weighted sum of squared residuals:

    ∑ w_i · (β_Yi - β_IVW · β_Xi)²

The closed-form solution is:

             ⎡  ∑ (w_i · β_Xi · β_Yi) ⎤
    β_IVW = ⎢ ------------------------ ⎥
             ⎣   ∑ (w_i · β_Xi²)      ⎦

Where the weights are the inverse variance of the outcome effects:

    w_i = 1 / SE_Yi²

---

### 🧮 Weighted Median Estimator (WME)

The Weighted Median Estimator gives a consistent estimate even when up to 50% of the instruments are invalid.

For each SNP, the ratio estimate is:

    β_i = β_Yi / β_Xi

And the weights are:

    w_i = 1 / SE_Yi²

The estimator computes the **weighted median** of the β_i values — the value where 50% of the total weight lies on either side. This is robust to violations of the exclusion restriction, as long as over 50% of the total weight comes from valid instruments.

---

### 📐 MR-Egger Regression

MR-Egger extends the IVW approach by including an intercept to account for directional (unbalanced) pleiotropy:

    β_Yi = α + β_Egger · β_Xi + ε_i

Where:
- α is the **intercept**, which captures the average pleiotropic effect
- β_Egger is the **causal effect estimate**, corrected for pleiotropy

This method relies on the **InSIDE** assumption: the Instrument Strength is Independent of the Direct Effect.

---

## 📊 Instrument Strength

### F-statistic (per SNP)

Used to assess the strength of the genetic instruments:

    F = (β_Xi²) / (SE_Xi²)

A rule of thumb: F > 10 suggests the instrument is strong.

---

## 🔬 Heterogeneity & Pleiotropy Tests

### Cochran’s Q Statistic

Assesses heterogeneity among SNP-specific causal estimates:

    Q = ∑ w_i · (β_Yi - β_IVW · β_Xi)²

Large Q suggests potential invalid instruments or pleiotropy.

---

### I² Statistic (MR-Egger)

Measures heterogeneity in instrument strength:

    I² = 1 - (1 / mean(F))

Low I² may indicate regression dilution bias in MR-Egger.

---

### MR-Egger Intercept Test

Used to detect **directional pleiotropy**:

    Z_i = α + β_Egger · (1 / SE_i) + ε_i

Where:
- Z_i is the standardized SNP-outcome effect
- α is the Egger intercept (should be ≈ 0 under no pleiotropy)

A significant intercept indicates directional pleiotropy bias.

---

## 📊 Causal Effect Estimation: IVW vs. WME vs. MR-Egger

This plot visualizes all three Mendelian Randomisation (MR) methods estimating the causal effect.

![IVW vs WME vs Egger](https://github.com/user-attachments/assets/91b3ba23-4560-407a-ae52-39545bda4dbb)

### 📌 Interpretation

- **Red line** – *Inverse-Variance Weighted (IVW)*
- **Green line** – *Weighted Median Estimator (WME)*
- **Blue line** – *MR-Egger Regression*
- **Orange dashed line** – *MR-Egger intercept*
- **Data points**: Each dot represents a single SNP's ratio estimate (β_Yi / β_Xi), with **bubble size** reflecting inverse-variance weight (1/SE²).

---
