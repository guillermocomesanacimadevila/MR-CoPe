# Install required packages
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR@0.4.26")

# Load package
library(TwoSampleMR)

# Set working directory
setwd("/Users/guillermocomesanacimadevila/Desktop/CRyPTIC_cleaning/MR R")

# Load CSV
filtered_snps <- read.csv("filtered_SNPs_for_MR.csv", header=TRUE)
filtered_snps

# Check is has been imported correctly
dim(filtered_snps)
print(colnames(filtered_snps))

# Ensure column names are clean
colnames(filtered_snps) <- gsub("\\s+", "", colnames(filtered_snps))
print(colnames(filtered_snps))

# Rename columns to match TwoSampleMR requirements
colnames(filtered_snps)[colnames(filtered_snps) == "A1_exp"] <- "effect_allele"
colnames(filtered_snps)[colnames(filtered_snps) == "A2_exp"] <- "other_allele"
colnames(filtered_snps)[colnames(filtered_snps) == "EAF_exp"] <- "eaf"

colnames(filtered_snps)[colnames(filtered_snps) == "A1_out"] <- "effect_allele.outcome"
colnames(filtered_snps)[colnames(filtered_snps) == "A2_out"] <- "other_allele.outcome"
colnames(filtered_snps)[colnames(filtered_snps) == "EAF_out"] <- "eaf.outcome"

# Ensure required columns exist
required_cols <- c("SNP", "BETA_exp", "SE_exp", "PVALUE_exp", "effect_allele", "other_allele", "eaf")
missing_cols <- setdiff(required_cols, colnames(filtered_snps))

if (length(missing_cols) > 0) {
  stop(paste("ERROR: Missing required columns:", paste(missing_cols, collapse = ", ")))
}

# Define Exposure and Outcome datasets
exposure_dat <- na.omit(
  data.frame(
    SNP = filtered_snps$SNP,
    beta = as.numeric(filtered_snps$BETA_exp),
    se = as.numeric(filtered_snps$SE_exp),
    pval = as.numeric(filtered_snps$PVALUE_exp),
    eaf = as.numeric(filtered_snps$eaf),
    effect_allele = as.character(filtered_snps$effect_allele),
    other_allele = as.character(filtered_snps$other_allele)
  )
)

outcome_dat <- na.omit(
  data.frame(
    SNP = filtered_snps$SNP,
    beta = as.numeric(filtered_snps$BETA_out),
    se = as.numeric(filtered_snps$SE_out),
    pval = as.numeric(filtered_snps$PVALUE_out),
    eaf = as.numeric(filtered_snps$eaf.outcome),
    effect_allele = as.character(filtered_snps$effect_allele.outcome),
    other_allele = as.character(filtered_snps$other_allele.outcome)
  )
)

# Save formatted exposure and outcome data
write.csv(exposure_dat, "exposure_dat.csv", row.names = FALSE)
write.csv(outcome_dat, "outcome_dat.csv", row.names = FALSE)

# Format data for MR analysis
exposure_dat <- format_data(exposure_dat, type = "exposure")
outcome_dat <- format_data(outcome_dat, type = "outcome")

# Merge datasets manually (skip harmonise_data)
harmonised_data <- merge(exposure_dat, outcome_dat, by = "SNP")
harmonised_data$mr_keep <- TRUE 

# Run MR methods
res_ivw <- mr(harmonised_data, method = "mr_ivw")  
res_egger <- mr(harmonised_data, method = "mr_egger_regression")  
res_wmedian <- mr(harmonised_data, method = "mr_weighted_median")  

# Sensitivity analyses
heterogeneity <- mr_heterogeneity(harmonised_data)
pleiotropy <- mr_pleiotropy_test(harmonised_data)

# Combine results into a table
results <- data.frame(matrix(ncol = 19, nrow = 1))

# Assign column names
colnames(results) <- c(
  "Exposure", "Outcome", "N_SNPs", 
  "IVW_OR", "IVW_Lower_95%", "IVW_Upper_95%", "IVW_Pval", 
  "WM_OR", "WM_Lower_95%", "WM_Upper_95%", "WM_Pval", 
  "Egger_OR", "Egger_Lower_95%", "Egger_Upper_95%", "Egger_Pval",
  "Egger_Intercept_Pval", "IVW_Q_Statistic_P", "I2_Statistic", "Mean_F_Statistic"
)

# Assign values to the table
results[1,1] <- "LDL-C"
results[1,2] <- "Alzheimer's Disease"
results[1,3] <- nrow(harmonised_data)

# Assign IVW results
if (!is.null(res_ivw$b)) {
  results[1,4:7] <- c(exp(res_ivw$b), exp(res_ivw$b + qnorm(.025) * res_ivw$se), exp(res_ivw$b + qnorm(.975) * res_ivw$se), res_ivw$pval)
} else {
  results[1,4:7] <- NA
}

# Assign Weighted Median results
if (!is.null(res_wmedian$b)) {
  results[1,8:11] <- c(exp(res_wmedian$b), exp(res_wmedian$b + qnorm(.025) * res_wmedian$se), exp(res_wmedian$b + qnorm(.975) * res_wmedian$se), res_wmedian$pval)
} else {
  results[1,8:11] <- NA
}

# Assign Egger results
if (!is.null(res_egger$b)) {
  results[1,12:15] <- c(exp(res_egger$b), exp(res_egger$b + qnorm(.025) * res_egger$se), exp(res_egger$b + qnorm(.975) * res_egger$se), res_egger$pval)
} else {
  results[1,12:15] <- NA
}

results[1,16] <- pleiotropy$pval 
results[1,17] <- ifelse("Inverse variance weighted" %in% heterogeneity$method, heterogeneity$Q_pval[heterogeneity$method == "Inverse variance weighted"], NA)

if ("Inverse variance weighted" %in% heterogeneity$method && !is.null(heterogeneity$I2[heterogeneity$method == "Inverse variance weighted"])) {
  results[1,18] <- heterogeneity$I2[heterogeneity$method == "Inverse variance weighted"]
} else {
  results[1,18] <- NA
}

results[1,19] <- mean(harmonised_data$F_stat, na.rm = TRUE)

# Save results
write.csv(results, "MR_Formatted_Results.csv", row.names = FALSE)

print("MR Analysis completed successfully. Results saved to 'MR_Results_Per_SNP.csv'.")