rm(list = ls())
# Load the ggplot2 package
library(ggplot2)

# Read the data into a data frame
df <- read.csv("out.snpden_chr5", sep = "\t", header = TRUE)

# Ensure the BIN_START and SNP_COUNT columns are numeric
df$BIN_START <- as.numeric(df$BIN_START)
df$SNP_COUNT <- as.numeric(df$SNP_COUNT)

# Create a bar plot of SNP count
ggplot(df, aes(x = BIN_START, y = SNP_COUNT)) +
  geom_bar(stat = "identity") +
  labs(x = "Bin Start Position", y = "SNP Count", title = "chr5 SNP Density")

# Read the data into a data frame
df <- read.csv("out.snpden_chrZ", sep = "\t", header = TRUE)

# Ensure the BIN_START and SNP_COUNT columns are numeric
df$BIN_START <- as.numeric(df$BIN_START)
df$SNP_COUNT <- as.numeric(df$SNP_COUNT)

# Create a bar plot of SNP count
ggplot(df, aes(x = BIN_START, y = SNP_COUNT)) +
  geom_bar(stat = "identity") +
  labs(x = "Bin Start Position", y = "SNP Count", title = "chrZ SNP Density")

