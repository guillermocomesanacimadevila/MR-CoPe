# === Install packages === #
install.packages("ggplot2")
install.packages("dplyr")

# === Call libraries === #
library(ggplot2)
library(dplyr)

# Set working directory
setwd("/Users/guillermocomesanacimadevila/cpep_MR/Cleaned Data")

# Load exposure GWAS data
exposure <- read.csv("exp_stats.csv", header=TRUE)

# Ensure P_VALUE is numeric
exposure$P_VALUE <- as.numeric(exposure$P_VALUE)

# Convert P_VALUE to -log10 scale
exposure$logP <- -log10(exposure$P_VALUE)

# Create scatter plot
plot <- ggplot(exposure, aes(x = POS, y = logP)) +
  geom_point(alpha = 0.6, color = "#0072B2", size = 1.2) +  
  theme_minimal(base_size = 14) +  
  labs(
    title = "Manhattan-style Scatter Plot",
    x = "Genomic Position (POS)",
    y = expression(-log[10](P)),
    caption = "Source: exp_stats.csv"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.major = element_line(color = "grey80", linetype = "dashed"),
    panel.grid.minor = element_blank()
  ) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red")  

# Display plot
print(plot)

# Outcome dataset
outcome <- read.csv("outcome_stats.csv", header=TRUE)

# Ensure P_VAL = NUMERIC
outcome$P_VALUE <- as.numeric(outcome$P_VALUE)

# Convert P_VALUE to -log10 scale
outcome$logP <- -log10(outcome$P_VALUE)

# Create scatter plot
plot2 <- ggplot(outcome, aes(x = POS, y = logP)) +
  geom_point(alpha = 0.6, color = "#0072B2", size = 1.2) +  
  theme_minimal(base_size = 14) +  
  labs(
    title = "Manhattan-style Scatter Plot",
    x = "Genomic Position (POS)",
    y = expression(-log[10](P)),
    caption = "Source: outcome_stats.csv"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.major = element_line(color = "grey80", linetype = "dashed"),
    panel.grid.minor = element_blank()
  ) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red")  

# Display plot
print(plot2)